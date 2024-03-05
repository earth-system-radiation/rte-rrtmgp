
#include "mo_gas_optics_kernels.h"
#include "rrtmgp_conversion.h"
#include <limits>


void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp, int2d const &flavor,
                   real1d const &press_ref_log, real1d const &temp_ref, real press_ref_log_delta, real temp_ref_min,
                   real temp_ref_delta, real press_ref_trop_log, real3d const &vmr_ref, real2d const &play,
                   real2d const &tlay, real3d const &col_gas, int2d const &jtemp, real6d const &fmajor, real5d const &fminor,
                   real4d const &col_mix, bool2d const &tropo, int4d const &jeta, int2d const &jpress) {
  using yakl::SB;
  using yakl::intrinsics::merge;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real2d ftemp ("ftemp" ,ncol,nlay);
  real2d fpress("fpress",ncol,nlay);

  real tiny = std::numeric_limits<real>::min();

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
    // index and factor for temperature interpolation
    jtemp(icol,ilay) = (int) ((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta);
    jtemp(icol,ilay) = std::min(ntemp - 1, std::max(1, jtemp(icol,ilay))); // limit the index range
    ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta;

    // index and factor for pressure interpolation
    real locpress = 1. + (log(play(icol,ilay)) - press_ref_log(1)) / press_ref_log_delta;
    jpress(icol,ilay) = std::min(npres-1, std::max(1, (int)(locpress)));
    fpress(icol,ilay) = locpress - (real)(jpress(icol,ilay));

    // determine if in lower or upper part of atmosphere
    tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log;
  });

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int iflav=1; iflav<=nflav; iflav++) {   // loop over implemented combinations of major species
  //       for (int itemp=1; itemp<=2; itemp++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nlay,ncol,nflav,2) , YAKL_LAMBDA (int ilay, int icol, int iflav , int itemp) {
    yakl::FSArray<int,1,SB<2>> igases;

    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));
    igases(1) = flavor(1,iflav);
    igases(2) = flavor(2,iflav);

    // compute interpolation fractions needed for lower, then upper reference temperature level
    // compute binary species parameter (eta) for flavor and temperature and
    //  associated interpolation index and factors
    real ratio_eta_half = vmr_ref(itropo,igases(1),(jtemp(icol,ilay)+itemp-1)) /
                          vmr_ref(itropo,igases(2),(jtemp(icol,ilay)+itemp-1));
    col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases(1)) + ratio_eta_half * col_gas(icol,ilay,igases(2));
    real eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(itemp,iflav,icol,ilay), 0.5,
                     col_mix(itemp,iflav,icol,ilay) > 2. * tiny);
    real loceta = eta * (neta-1.0);
    jeta(itemp,iflav,icol,ilay) = std::min((int)(loceta)+1, neta-1);
    real feta = fmod(loceta, 1.0);
    // compute interpolation fractions needed for minor species
    real ftemp_term = ((2.0 - itemp) + (2.0 * itemp - 3.0 ) * ftemp(icol,ilay));
    fminor(1,itemp,iflav,icol,ilay) = (1. - feta) * ftemp_term;
    fminor(2,itemp,iflav,icol,ilay) =       feta  * ftemp_term;
    // compute interpolation fractions needed for major species
    fmajor(1,1,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay);
    fmajor(2,1,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(2,itemp,iflav,icol,ilay);
    fmajor(1,2,itemp,iflav,icol,ilay) =       fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay);
    fmajor(2,2,itemp,iflav,icol,ilay) =       fpress(icol,ilay)  * fminor(2,itemp,iflav,icol,ilay);
  });
}



void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real3d const &tau_abs, real3d const &tau_rayleigh,
                              real3d const &tau, real3d const &ssa, real3d const &g) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real tiny = std::numeric_limits<real>::min();

  int constexpr TILE_SIZE=8;
  int colTiles = ncol / TILE_SIZE + 1;
  int gptTiles = ngpt / TILE_SIZE + 1;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int tcol=1; tcol<=colTiles; tcol++) {
  //     for (int tgpt=1; tgpt<=gptTiles; tgpt++) {
  //       for (int itcol=1; itcol<=TILE_SIZE; itcol++) {
  //         for (int itgpt=1; itgpt<=TILE_SIZE; itgpt++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<5>(nlay,colTiles,gptTiles,TILE_SIZE,TILE_SIZE) , YAKL_LAMBDA (int ilay, int tcol, int tgpt, int itcol, int itgpt) {
    int icol = (tcol-1)*TILE_SIZE + itcol;
    int igpt = (tgpt-1)*TILE_SIZE + itgpt;

    if ( icol <= ncol && igpt <= ngpt ) {
      real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
      tau(icol,ilay,igpt) = t;
      g  (icol,ilay,igpt) = 0.;
      if(t > 2. * tiny) {
        ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
      } else {
        ssa(icol,ilay,igpt) = 0.;
      }
    }
  });
}



void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
                           real2d const &tlay, real2d const &tlev, real1d const &tsfc, int sfc_lay, real6d const &fmajor,
                           int4d const &jeta, bool2d const &tropo, int2d const &jtemp, int2d const &jpress,
                           int1d const &gpoint_bands, int2d const &band_lims_gpt, real4d const &pfracin, real temp_ref_min,
                           real totplnk_delta, real2d const &totplnk, int2d const &gpoint_flavor, real2d const &sfc_src,
                           real3d const &lay_src, real3d const &lev_src_inc, real3d const &lev_src_dec) {
  using yakl::COLON;
  using yakl::intrinsics::merge;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;
  using yakl::fortran::Bounds;

  real3d pfrac          ("pfrac"          ,ngpt,nlay,ncol);
  real3d planck_function("planck_function",nbnd,nlay+1,ncol);
  real1d one            ("one"            ,2);
  one = 1;

  // Calculation of fraction of band's Planck irradiance associated with each g-point
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nlay,ncol,ngpt) , YAKL_LAMBDA (int ilay, int icol, int igpt) {
    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));  //WS moved itropo inside loop for GPU
    int iflav = gpoint_flavor(itropo, igpt); //eta interpolation depends on band's flavor
    // interpolation in temperature, pressure, and eta
    int jpress_loc = jpress(icol,ilay)+itropo;
    int jtemp_loc  = jtemp(icol,ilay);

    // inlining interpolate3D
    pfrac(igpt,ilay,icol) = ( fmajor(1,1,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc  ) +
                              fmajor(2,1,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc  ) +
                              fmajor(1,2,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc  , jtemp_loc  ) +
                              fmajor(2,2,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc  ) ) +
                            ( fmajor(1,1,2,iflav,icol,ilay) * pfracin(igpt, jeta(2,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc+1) +
                              fmajor(2,1,2,iflav,icol,ilay) * pfracin(igpt, jeta(2,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc+1) +
                              fmajor(1,2,2,iflav,icol,ilay) * pfracin(igpt, jeta(2,iflav,icol,ilay)  , jpress_loc  , jtemp_loc+1) +
                              fmajor(2,2,2,iflav,icol,ilay) * pfracin(igpt, jeta(2,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc+1) );
  });

  //
  // Planck function by band for the surface
  // Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  //
  // for (int icol=1; icol<=ncol; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , ncol , YAKL_LAMBDA (int icol) {
    interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk, planck_function.slice<1>(COLON,1,icol),nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // for (int igpt=1; igpt<=ngpt; igpt++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
    sfc_src(igpt,icol) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 1, icol);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ncol,nlay) , YAKL_LAMBDA (int icol, int ilay) {
    // Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
    interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk, planck_function.slice<1>(COLON,ilay,icol),nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // Explicitly unroll a time-consuming loop here to increase instruction-level parallelism on a GPU
  // Helps to achieve higher bandwidth
  //
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    lay_src(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,icol);
  });

  // compute level source irradiances for each g-point, one each for upward and downward paths
  // for (int icol=1; icol<=ncol; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , ncol , YAKL_LAMBDA (int icol) {
    interpolate1D(tlev(icol,1), temp_ref_min, totplnk_delta, totplnk, planck_function.slice<1>(COLON,1,icol),nPlanckTemp,nbnd);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=2; ilay<=nlay+1; ilay++) {
  parallel_for( YAKL_AUTO_LABEL() , Bounds<2>(ncol,{2,nlay+1}) , YAKL_LAMBDA (int icol, int ilay) {
    interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk, planck_function.slice<1>(COLON,ilay,icol),nPlanckTemp,nbnd);
  });

  //
  // Map to g-points
  //
  // Same unrolling as mentioned before
  //
  // for (int icol=1; icol<=ncol; icol+=2) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  );
    lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  );
    if (icol < ncol) {
      lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1);
      lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1);
    }
  });
}



// compute Rayleigh scattering optical depths
void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp,
                          int2d const &gpoint_flavor, int2d const &band_lims_gpt, real4d const &krayl, int idx_h2o,
                          real2d const &col_dry, real3d const &col_gas, real5d const &fminor, int4d const &jeta,
                          bool2d const &tropo, int2d const &jtemp, real3d const &tau_rayleigh) {
  using yakl::intrinsics::merge;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nlay,ncol,ngpt) , YAKL_LAMBDA (int ilay, int icol, int igpt) {
    int itropo = merge(1,2,tropo(icol,ilay)); // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int iflav = gpoint_flavor(itropo, igpt);
    // Inlining interpolate2D
    real k = fminor(1,1,iflav,icol,ilay) * krayl(igpt, jeta(1,iflav,icol,ilay)  , jtemp(icol,ilay)  ,itropo) +
             fminor(2,1,iflav,icol,ilay) * krayl(igpt, jeta(1,iflav,icol,ilay)+1, jtemp(icol,ilay)  ,itropo) +
             fminor(1,2,iflav,icol,ilay) * krayl(igpt, jeta(2,iflav,icol,ilay)  , jtemp(icol,ilay)+1,itropo) +
             fminor(2,2,iflav,icol,ilay) * krayl(igpt, jeta(2,iflav,icol,ilay)+1, jtemp(icol,ilay)+1,itropo);
    tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay));
  });
}



// compute minor species optical depths
#ifdef RRTMGP_CPU_KERNELS

  void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                                int nminor, int nminork, int idx_h2o, int idx_tropo, int2d const &gpt_flv,
                                real3d const &kminor, int2d const &minor_limits_gpt, bool1d const &minor_scales_with_density,
                                bool1d const &scale_by_complement, int1d const &idx_minor, int1d const &idx_minor_scaling,
                                int1d const &kminor_start, real2d const &play, real2d const &tlay, real3d const &col_gas,
                                real5d const &fminor, int4d const &jeta, int2d const &layer_limits, int2d const &jtemp, real3d const &tau) {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;
    using yakl::fortran::Bounds;

    #ifdef YAKL_AUTO_PROFILE
      auto timername = std::string(YAKL_AUTO_LABEL());
      yakl::timer_start(timername.c_str());
    #endif

    real constexpr PaTohPa = 0.01;

    int extent = size(scale_by_complement,1);

    #ifdef YAKL_ARCH_OPENMP
      #pragma omp parallel for collapse(2)
    #endif
    for (int icol=1; icol<=ncol; icol++) {
      for (int ilay=1; ilay<=nlay; ilay++) {
        // This check skips individual columns with no pressures in range
        if ( layer_limits(icol,1) <= 0 || ilay < layer_limits(icol,1) || ilay > layer_limits(icol,2) ) continue;

        int myjtemp = jtemp(icol,ilay);

        for (int imnr=1; imnr<=extent; imnr++) {
          real scaling = col_gas(icol,ilay,idx_minor(imnr));
          if (minor_scales_with_density(imnr)) {
            // NOTE: P needed in hPa to properly handle density scaling.
            scaling = scaling * (PaTohPa * play(icol,ilay)/tlay(icol,ilay));

            if (idx_minor_scaling(imnr) > 0) {  // there is a second gas that affects this gas's absorption
              real mycol_gas_imnr = col_gas(icol,ilay,idx_minor_scaling(imnr));
              real vmr_fact = 1. / col_gas(icol,ilay,0);
              real dry_fact = 1. / (1. + col_gas(icol,ilay,idx_h2o) * vmr_fact);
              // scale by density of special gas
              if (scale_by_complement(imnr)) { // scale by densities of all gases but the special one
                scaling = scaling * (1. - mycol_gas_imnr * vmr_fact * dry_fact);
              } else {
                scaling = scaling *          mycol_gas_imnr * vmr_fact * dry_fact;
              }
            }
          } // minor_scalse_with_density(imnr)

          // What is the starting point in the stored array of minor absorption coefficients?
          int minor_start = kminor_start(imnr);
          // Which gpoint range does this minor gas affect?
          int gptS = minor_limits_gpt(1,imnr);
          int gptE = minor_limits_gpt(2,imnr);
          for (int igpt=gptS; igpt<=gptE; igpt++) {
            // Interpolation of absorption coefficient and calculation of optical depth
            int iflav = gpt_flv(idx_tropo,igpt); // eta interpolation depends on flavor
            int minor_loc = minor_start + (igpt - gptS); // add offset to starting point
            // Inlined interpolate2D
            real kminor_loc = fminor(1,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)  , myjtemp  ) +
                              fminor(2,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)+1, myjtemp  ) +
                              fminor(1,2,iflav,icol,ilay) * kminor(minor_loc, jeta(2,iflav,icol,ilay)  , myjtemp+1) +
                              fminor(2,2,iflav,icol,ilay) * kminor(minor_loc, jeta(2,iflav,icol,ilay)+1, myjtemp+1);
            tau(igpt,ilay,icol) += kminor_loc * scaling;
          }
        }
      }
    }
    #ifdef YAKL_AUTO_PROFILE
      yakl::timer_stop(timername.c_str());
    #endif
  }

#else

  void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                                int nminor, int nminork, int idx_h2o, int idx_tropo, int2d const &gpt_flv,
                                real3d const &kminor, int2d const &minor_limits_gpt, bool1d const &minor_scales_with_density,
                                bool1d const &scale_by_complement, int1d const &idx_minor, int1d const &idx_minor_scaling,
                                int1d const &kminor_start, real2d const &play, real2d const &tlay, real3d const &col_gas,
                                real5d const &fminor, int4d const &jeta, int2d const &layer_limits, int2d const &jtemp, real3d const &tau) {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;
    using yakl::fortran::Bounds;

    real constexpr PaTohPa = 0.01;

    int extent = size(scale_by_complement,1);

    // for (int ilay=1; ilay<=nlay; ilay++) {
    //   for (int icol=1; icol<=ncol; icol++) {
    //     for (int igpt0=0; igpt0<=max_gpt_diff; igpt0++) {
    parallel_for( YAKL_AUTO_LABEL() , Bounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      // This check skips individual columns with no pressures in range
      if ( layer_limits(icol,1) <= 0 || ilay < layer_limits(icol,1) || ilay > layer_limits(icol,2) ) {
      } else {
        real myplay  = play (icol,ilay);
        real mytlay  = tlay (icol,ilay);
        int  myjtemp = jtemp(icol,ilay);
        real mycol_gas_h2o = col_gas(icol,ilay,idx_h2o);
        real mycol_gas_0   = col_gas(icol,ilay,0);

        for (int imnr=1; imnr<=extent; imnr++) {
          real scaling = col_gas(icol,ilay,idx_minor(imnr));
          if (minor_scales_with_density(imnr)) {
            // NOTE: P needed in hPa to properly handle density scaling.
            scaling = scaling * (PaTohPa * myplay/mytlay);

            if (idx_minor_scaling(imnr) > 0) {  // there is a second gas that affects this gas's absorption
              real mycol_gas_imnr = col_gas(icol,ilay,idx_minor_scaling(imnr));
              real vmr_fact = 1. / mycol_gas_0;
              real dry_fact = 1. / (1. + mycol_gas_h2o * vmr_fact);
              // scale by density of special gas
              if (scale_by_complement(imnr)) { // scale by densities of all gases but the special one
                scaling = scaling * (1. - mycol_gas_imnr * vmr_fact * dry_fact);
              } else {
                scaling = scaling *          mycol_gas_imnr * vmr_fact * dry_fact;
              }
            }
          } // minor_scalse_with_density(imnr)

          // What is the starting point in the stored array of minor absorption coefficients?
          int minor_start = kminor_start(imnr);
          // Which gpoint range does this minor gas affect?
          int gptS = minor_limits_gpt(1,imnr);
          int gptE = minor_limits_gpt(2,imnr);
          for (int igpt=gptS; igpt<=gptE; igpt++) {
            // Interpolation of absorption coefficient and calculation of optical depth
            int iflav = gpt_flv(idx_tropo,igpt); // eta interpolation depends on flavor
            int minor_loc = minor_start + (igpt - gptS); // add offset to starting point
            // Inlined interpolate2D
            real kminor_loc = fminor(1,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)  , myjtemp  ) +
                              fminor(2,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)+1, myjtemp  ) +
                              fminor(1,2,iflav,icol,ilay) * kminor(minor_loc, jeta(2,iflav,icol,ilay)  , myjtemp+1) +
                              fminor(2,2,iflav,icol,ilay) * kminor(minor_loc, jeta(2,iflav,icol,ilay)+1, myjtemp+1);
            real tau_minor = kminor_loc * scaling;
            tau(igpt,ilay,icol) += tau_minor;
          }
        }
      }
    });
  }

#endif


// compute minor species optical depths
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2d const &gpoint_flavor, int2d const &band_lims_gpt, real4d const &kmajor,
                              real4d const &col_mix, real6d const &fmajor, int4d const &jeta, bool2d const &tropo,
                              int2d const &jtemp, int2d const &jpress, real3d const &tau) {
  using yakl::intrinsics::merge;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  // optical depth calculation for major species
  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     // optical depth calculation for major species
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nlay,ncol,ngpt) , YAKL_LAMBDA (int ilay, int icol, int igpt) {
    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));  // WS: moved inside innermost loop

    // binary species parameter (eta) and col_mix depend on band flavor
    int iflav = gpoint_flavor(itropo, igpt);
    // interpolation in temperature, pressure, and eta
    int jpress_loc = jpress(icol,ilay)+itropo;
    int jtemp_loc  = jtemp (icol,ilay);

    // inlined interpolate3D
    real tau_major = col_mix(1,iflav,icol,ilay) * ( fmajor(1,1,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc  ) +
                                                    fmajor(2,1,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc  ) +
                                                    fmajor(1,2,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc  , jtemp_loc  ) +
                                                    fmajor(2,2,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc  ) ) +
                     col_mix(2,iflav,icol,ilay) * ( fmajor(1,1,2,iflav,icol,ilay) * kmajor(igpt, jeta(2,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc+1) +
                                                    fmajor(2,1,2,iflav,icol,ilay) * kmajor(igpt, jeta(2,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc+1) +
                                                    fmajor(1,2,2,iflav,icol,ilay) * kmajor(igpt, jeta(2,iflav,icol,ilay)  , jpress_loc  , jtemp_loc+1) +
                                                    fmajor(2,2,2,iflav,icol,ilay) * kmajor(igpt, jeta(2,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc+1) );

    tau(igpt,ilay,icol) += tau_major;
   });
}



// Compute minor and major species opitcal depth from pre-computed interpolation coefficients
//   (jeta,jtemp,jpress)
void compute_tau_absorption(int max_gpt_diff_lower, int max_gpt_diff_upper, int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp, int nminorlower,
                            int nminorklower, int nminorupper, int nminorkupper, int idx_h2o, int2d const &gpoint_flavor,
                            int2d const &band_lims_gpt, real4d const &kmajor, real3d const &kminor_lower, real3d const &kminor_upper,
                            int2d const &minor_limits_gpt_lower, int2d const &minor_limits_gpt_upper, bool1d const &minor_scales_with_density_lower,
                            bool1d const &minor_scales_with_density_upper, bool1d const &scale_by_complement_lower,
                            bool1d const &scale_by_complement_upper, int1d const &idx_minor_lower, int1d const &idx_minor_upper,
                            int1d const &idx_minor_scaling_lower, int1d const &idx_minor_scaling_upper, int1d const &kminor_start_lower,
                            int1d const &kminor_start_upper, bool2d const &tropo, real4d const &col_mix, real6d const &fmajor,
                            real5d const &fminor, real2d const &play, real2d const &tlay, real3d const &col_gas, int4d const &jeta,
                            int2d const &jtemp, int2d const &jpress, real3d const &tau, bool top_at_1) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  int2d itropo_lower("itropo_lower",ncol,2);
  int2d itropo_upper("itropo_upper",ncol,2);

  int huge  = std::numeric_limits<int>::max();
  int small = std::numeric_limits<int>::min();

  if (top_at_1) {

    // for (int icol=1; icol<=ncol; icol++){
    parallel_for( YAKL_AUTO_LABEL() , ncol , YAKL_LAMBDA (int icol) {
      itropo_lower(icol,2) = nlay;
      // itropo_lower(icol,1) = minloc(play(icol,:), dim=1, mask=tropo(icol,:))
      {
        int minloc = huge;
        real mn = (real) huge;
        for (int i=1; i<=nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = std::min(mn,play(icol,i));
          }
        }
        itropo_lower(icol,1) = minloc;
      }

      itropo_upper(icol,1) = 1;
      // itropo_upper(icol,2) = maxloc(play(icol,:), dim=1, mask=(.not. tropo(icol,:)))
      {
        int maxloc = small;
        real mx = (real) small;
        for (int i=1; i<=nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = std::max(mx,play(icol,i));
          }
        }
        itropo_upper(icol,2) = maxloc;
      }
    });

  } else {  // top_at_1

    // for (int icol=1; icol<=ncol; icol++){
    parallel_for( YAKL_AUTO_LABEL() , ncol , YAKL_LAMBDA ( int icol ) {
      itropo_lower(icol,1) = 1;
      // itropo_lower(icol,2) = minloc(play(icol,:), dim=1, mask=tropo(icol,:))
      {
        int minloc = huge;
        real mn = (real) huge;
        for (int i=1; i<=nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = std::min(mn,play(icol,i));
          }
        }
        itropo_lower(icol,2) = minloc;
      }

      itropo_upper(icol,2) = nlay;
      // itropo_upper(icol,1) = maxloc(play(icol,:), dim=1, mask=(.not.tropo(icol,:)))
      {
        int maxloc = small;
        real mx = (real) small;
        for (int i=1; i<=nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = std::max(mx,play(icol,i));
          }
        }
        itropo_upper(icol,1) = maxloc;
      }

    });

  }  // top_at_1

  // ---------------------
  // Major Species
  // ---------------------
  gas_optical_depths_major(ncol, nlay, nbnd, ngpt, nflav, neta, npres, ntemp, gpoint_flavor,
                           band_lims_gpt, kmajor, col_mix, fmajor, jeta, tropo, jtemp, jpress, tau);
  // ---------------------
  // Minor Species - lower
  // ---------------------
  int idx_tropo = 1;
  gas_optical_depths_minor(max_gpt_diff_lower, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorlower, nminorklower,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_lower, minor_limits_gpt_lower,
                           minor_scales_with_density_lower, scale_by_complement_lower, idx_minor_lower,
                           idx_minor_scaling_lower, kminor_start_lower, play,  tlay, col_gas, fminor,
                           jeta, itropo_lower, jtemp, tau);
  // ---------------------
  // Minor Species - upper
  // ---------------------
  idx_tropo = 2;
  gas_optical_depths_minor(max_gpt_diff_upper, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorupper, nminorkupper,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_upper, minor_limits_gpt_upper,
                           minor_scales_with_density_upper, scale_by_complement_upper, idx_minor_upper,
                           idx_minor_scaling_upper, kminor_start_upper, play, tlay, col_gas, fminor,
                           jeta, itropo_upper, jtemp, tau);
}



// Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
//   using Rayleigh scattering phase function
void combine_and_reorder_nstr(int ncol, int nlay, int ngpt, int nmom, real3d const &tau_abs, real3d const &tau_rayleigh, real3d const &tau, real3d const &ssa, real4d const &p) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real tiny = std::numeric_limits<real>::min();

  // do icol = 1, ncol
  //   do ilay = 1, nlay
  //     do igpt = 1, ngpt
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
    tau(icol,ilay,igpt) = t;
    if (t > 2. * tiny) {
      ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
    } else {
      ssa(icol,ilay,igpt) = 0;
    }
    for (int imom=1; imom<=nmom; imom++) {
      p(icol,ilay,igpt,imom) = 0;
    }
    if (nmom >= 2) {
      p(icol,ilay,igpt,2) = 0.1;
    }
  });
}

#ifdef RRTMGP_ENABLE_KOKKOS
void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp, int2dk const &flavor,
                   real1dk const &press_ref_log, real1dk const &temp_ref, real press_ref_log_delta, real temp_ref_min,
                   real temp_ref_delta, real press_ref_trop_log, realOff3dk const &vmr_ref, real2dk const &play,
                   real2dk const &tlay, realOff3dk const &col_gas, int2dk const &jtemp, real6dk const &fmajor, real5dk const &fminor,
                   real4dk const &col_mix, bool2dk const &tropo, int4dk const &jeta, int2dk const &jpress) {
  using yakl::intrinsics::merge;

  real2dk ftemp ("ftemp" ,ncol,nlay);
  real2dk fpress("fpress",ncol,nlay);

  real tiny = std::numeric_limits<real>::min();

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}) , KOKKOS_LAMBDA (int ilay, int icol) {
    // index and factor for temperature interpolation
    jtemp(icol,ilay) = (int) ((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta);
    jtemp(icol,ilay) = std::min(ntemp - 1, std::max(1, jtemp(icol,ilay))) - 1; // limit the index range
    ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta;

    // index and factor for pressure interpolation
    real locpress = 1. + (log(play(icol,ilay)) - press_ref_log(0)) / press_ref_log_delta;
    jpress(icol,ilay) = std::min(npres-1, std::max(1, (int)(locpress)));
    fpress(icol,ilay) = locpress - (real)(jpress(icol,ilay));

    // determine if in lower or upper part of atmosphere
    tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log;
  });

  Kokkos::parallel_for( MDRangeP<4>({0,0,0,0}, {nlay,ncol,nflav,2}) , KOKKOS_LAMBDA (int ilay, int icol, int iflav , int itemp) {
    // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
    int itropo = merge(0,1,tropo(icol,ilay));
    auto igases1 = flavor(0,iflav);
    auto igases2 = flavor(1,iflav);

    // compute interpolation fractions needed for lower, then upper reference temperature level
    // compute binary species parameter (eta) for flavor and temperature and
    //  associated interpolation index and factors
    real ratio_eta_half = vmr_ref(itropo,igases1,(jtemp(icol,ilay)+itemp)) /
                          vmr_ref(itropo,igases2,(jtemp(icol,ilay)+itemp));
    col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases1) + ratio_eta_half * col_gas(icol,ilay,igases2);
    real eta = merge(col_gas(icol,ilay,igases1) / col_mix(itemp,iflav,icol,ilay), 0.5,
                     col_mix(itemp,iflav,icol,ilay) > 2. * tiny);
    real loceta = eta * (neta-1.0);
    jeta(itemp,iflav,icol,ilay) = std::min((int)(loceta)+1, neta-1) - 1;
    real feta = fmod(loceta, 1.0);
    // compute interpolation fractions needed for minor species
    real ftemp_term = ((2.0 - (itemp+1)) + (2.0 * (itemp+1) - 3.0 ) * ftemp(icol,ilay));
    fminor(0,itemp,iflav,icol,ilay) = (1. - feta) * ftemp_term;
    fminor(1,itemp,iflav,icol,ilay) =       feta  * ftemp_term;
    // compute interpolation fractions needed for major species
    fmajor(0,0,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(0,itemp,iflav,icol,ilay);
    fmajor(1,0,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay);
    fmajor(0,1,itemp,iflav,icol,ilay) =       fpress(icol,ilay)  * fminor(0,itemp,iflav,icol,ilay);
    fmajor(1,1,itemp,iflav,icol,ilay) =       fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay);
  });
}

void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real3dk const &tau_abs, real3dk const &tau_rayleigh,
                              real3dk const &tau, real3dk const &ssa, real3dk const &g) {
  real tiny = std::numeric_limits<real>::min();

  int constexpr TILE_SIZE=8;
  int colTiles = ncol / TILE_SIZE + 1;
  int gptTiles = ngpt / TILE_SIZE + 1;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int tcol=1; tcol<=colTiles; tcol++) {
  //     for (int tgpt=1; tgpt<=gptTiles; tgpt++) {
  //       for (int itcol=1; itcol<=TILE_SIZE; itcol++) {
  //         for (int itgpt=1; itgpt<=TILE_SIZE; itgpt++) {
  Kokkos::parallel_for( MDRangeP<5>({0,0,0,0,0}, {nlay,colTiles,gptTiles,TILE_SIZE,TILE_SIZE}) , KOKKOS_LAMBDA (int ilay, int tcol, int tgpt, int itcol, int itgpt) {
    int icol = tcol*TILE_SIZE + itcol;
    int igpt = tgpt*TILE_SIZE + itgpt;

    if ( icol < ncol && igpt < ngpt ) {
      real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
      tau(icol,ilay,igpt) = t;
      g  (icol,ilay,igpt) = 0.;
      if(t > 2. * tiny) {
        ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
      } else {
        ssa(icol,ilay,igpt) = 0.;
      }
    }
  });
}

void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
                           real2dk const &tlay, real2dk const &tlev, real1dk const &tsfc, int sfc_lay, real6dk const &fmajor,
                           int4dk const &jeta, bool2dk const &tropo, int2dk const &jtemp, int2dk const &jpress,
                           int1dk const &gpoint_bands, int2dk const &band_lims_gpt, real4dk const &pfracin, real temp_ref_min,
                           real totplnk_delta, real2dk const &totplnk, int2dk const &gpoint_flavor, real2dk const &sfc_src,
                           real3dk const &lay_src, real3dk const &lev_src_inc, real3dk const &lev_src_dec) {
  using yakl::intrinsics::merge;

  real3dk pfrac          ("pfrac"          ,ngpt,nlay,ncol);
  real3dk planck_function("planck_function",nbnd,nlay+1,ncol);
  real1dk one            ("one"            ,2);
  Kokkos::deep_copy(one, 1);

  // Calculation of fraction of band's Planck irradiance associated with each g-point
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0},{nlay,ncol,ngpt}) , KOKKOS_LAMBDA (int ilay, int icol, int igpt) {
    // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
    int itropo = merge(0,1,tropo(icol,ilay));  //WS moved itropo inside loop for GPU
    int iflav = gpoint_flavor(itropo, igpt); //eta interpolation depends on band's flavor
    // interpolation in temperature, pressure, and eta
    int jpress_loc = jpress(icol,ilay)+itropo;
    int jtemp_loc  = jtemp(icol,ilay);

    // inlining interpolate3D
    pfrac(igpt,ilay,icol) = (
       fmajor(0,0,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc  ) +
       fmajor(1,0,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc  ) +
       fmajor(0,1,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc  , jtemp_loc  ) +
       fmajor(1,1,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc  ) ) +
     ( fmajor(0,0,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc+1) +
       fmajor(1,0,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc+1) +
       fmajor(0,1,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc  , jtemp_loc+1) +
       fmajor(1,1,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc+1) );
  });

  //
  // Planck function by band for the surface
  // Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  //
  // for (int icol=1; icol<=ncol; icol++) {
  Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
      interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function, Kokkos::ALL, 0 ,icol),nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // for (int igpt=1; igpt<=ngpt; igpt++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
    sfc_src(igpt,icol) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 0, icol);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {ncol,nlay}) , KOKKOS_LAMBDA (int icol, int ilay) {
    // Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
      interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function, Kokkos::ALL,ilay,icol),nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // Explicitly unroll a time-consuming loop here to increase instruction-level parallelism on a GPU
  // Helps to achieve higher bandwidth
  //
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {ncol,nlay,ngpt}) , KOKKOS_LAMBDA (int icol, int ilay, int igpt) {
    lay_src(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,icol);
  });

  // compute level source irradiances for each g-point, one each for upward and downward paths
  // for (int icol=1; icol<=ncol; icol++) {
  Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
      interpolate1D(tlev(icol,0), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function, Kokkos::ALL,0 ,icol),nPlanckTemp,nbnd);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=2; ilay<=nlay+1; ilay++) {
  Kokkos::parallel_for( MDRangeP<2>({0, 1}, {ncol,nlay+1}) , KOKKOS_LAMBDA (int icol, int ilay) {
      interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function,Kokkos::ALL,ilay,icol),nPlanckTemp,nbnd);
  });

  //
  // Map to g-points
  //
  // Same unrolling as mentioned before
  //
  // for (int icol=1; icol<=ncol; icol+=2) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {ncol,nlay,ngpt}) , KOKKOS_LAMBDA (int icol, int ilay, int igpt) {
    lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  );
    lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  );
    if (icol < ncol-1) {
      lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1);
      lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1);
    }
  });
}



// compute Rayleigh scattering optical depths
void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp,
                          int2dk const &gpoint_flavor, int2dk const &band_lims_gpt, real4dk const &krayl, int idx_h2o,
                          real2dk const &col_dry, realOff3dk const &col_gas, real5dk const &fminor, int4dk const &jeta,
                          bool2dk const &tropo, int2dk const &jtemp, real3dk const &tau_rayleigh) {
  using yakl::intrinsics::merge;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {nlay,ncol,ngpt}) , KOKKOS_LAMBDA (int ilay, int icol, int igpt) {
    int itropo = merge(0,1,tropo(icol,ilay)); // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
    int iflav = gpoint_flavor(itropo, igpt);
    // Inlining interpolate2D
    real k = fminor(0,0,iflav,icol,ilay) * krayl(igpt, jeta(0,iflav,icol,ilay)  , jtemp(icol,ilay)  ,itropo) +
             fminor(1,0,iflav,icol,ilay) * krayl(igpt, jeta(0,iflav,icol,ilay)+1, jtemp(icol,ilay)  ,itropo) +
             fminor(0,1,iflav,icol,ilay) * krayl(igpt, jeta(1,iflav,icol,ilay)  , jtemp(icol,ilay)+1,itropo) +
             fminor(1,1,iflav,icol,ilay) * krayl(igpt, jeta(1,iflav,icol,ilay)+1, jtemp(icol,ilay)+1,itropo);
    tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay));
  });
}

void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, int2dk const &gpt_flv,
                              real3dk const &kminor, int2dk const &minor_limits_gpt, bool1dk const &minor_scales_with_density,
                              bool1dk const &scale_by_complement, int1dk const &idx_minor, int1dk const &idx_minor_scaling,
                              int1dk const &kminor_start, real2dk const &play, real2dk const &tlay, realOff3dk const &col_gas,
                              real5dk const &fminor, int4dk const &jeta, int2dk const &layer_limits, int2dk const &jtemp, real3dk const &tau) {
  real constexpr PaTohPa = 0.01;

  int extent = scale_by_complement.extent(0);

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt0=0; igpt0<=max_gpt_diff; igpt0++) {
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}) , KOKKOS_LAMBDA (int ilay, int icol) {
    // This check skips individual columns with no pressures in range
    if ( layer_limits(icol,0) <= -1 || ilay < layer_limits(icol,0) || ilay > layer_limits(icol,1) ) {
    } else {
      real myplay  = play (icol,ilay);
      real mytlay  = tlay (icol,ilay);
      int  myjtemp = jtemp(icol,ilay);
      real mycol_gas_h2o = col_gas(icol,ilay,idx_h2o);
      real mycol_gas_0   = col_gas(icol,ilay,-1);

      for (int imnr=0; imnr<extent; imnr++) {
        real scaling = col_gas(icol,ilay,idx_minor(imnr));
        if (minor_scales_with_density(imnr)) {
          // NOTE: P needed in hPa to properly handle density scaling.
          scaling = scaling * (PaTohPa * myplay/mytlay);

          if (idx_minor_scaling(imnr) > -1) {  // there is a second gas that affects this gas's absorption
            real mycol_gas_imnr = col_gas(icol,ilay,idx_minor_scaling(imnr));
            real vmr_fact = 1. / mycol_gas_0;
            real dry_fact = 1. / (1. + mycol_gas_h2o * vmr_fact);
            // scale by density of special gas
            if (scale_by_complement(imnr)) { // scale by densities of all gases but the special one
              scaling = scaling * (1. - mycol_gas_imnr * vmr_fact * dry_fact);
            } else {
              scaling = scaling *       mycol_gas_imnr * vmr_fact * dry_fact;
            }
          }
        } // minor_scalse_with_density(imnr)

        // What is the starting point in the stored array of minor absorption coefficients?
        int minor_start = kminor_start(imnr);
        // Which gpoint range does this minor gas affect?
        int gptS = minor_limits_gpt(0,imnr);
        int gptE = minor_limits_gpt(1,imnr);
        for (int igpt=gptS; igpt<=gptE; igpt++) {
          // Interpolation of absorption coefficient and calculation of optical depth
          int iflav = gpt_flv(idx_tropo,igpt); // eta interpolation depends on flavor
          int minor_loc = minor_start + (igpt - gptS); // add offset to starting point
          // Inlined interpolate2D
          real kminor_loc =
            fminor(0,0,iflav,icol,ilay) * kminor(minor_loc, jeta(0,iflav,icol,ilay)  , myjtemp  ) +
            fminor(1,0,iflav,icol,ilay) * kminor(minor_loc, jeta(0,iflav,icol,ilay)+1, myjtemp  ) +
            fminor(0,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)  , myjtemp+1) +
            fminor(1,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)+1, myjtemp+1);
          real tau_minor = kminor_loc * scaling;
          tau(igpt,ilay,icol) += tau_minor;
        }
      }
    }
  });
}

// compute minor species optical depths
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2dk const &gpoint_flavor, int2dk const &band_lims_gpt, real4dk const &kmajor,
                              real4dk const &col_mix, real6dk const &fmajor, int4dk const &jeta, bool2dk const &tropo,
                              int2dk const &jtemp, int2dk const &jpress, real3dk const &tau) {
  using yakl::intrinsics::merge;

  // optical depth calculation for major species
  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     // optical depth calculation for major species
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {nlay,ncol,ngpt}) , KOKKOS_LAMBDA (int ilay, int icol, int igpt) {
    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(0,1,tropo(icol,ilay));  // WS: moved inside innermost loop

    // binary species parameter (eta) and col_mix depend on band flavor
    int iflav = gpoint_flavor(itropo, igpt);
    // interpolation in temperature, pressure, and eta
    int jpress_loc = jpress(icol,ilay)+itropo;
    int jtemp_loc  = jtemp (icol,ilay);

    // inlined interpolate3D
    real tau_major = col_mix(0,iflav,icol,ilay) * (
      fmajor(0,0,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc  ) +
      fmajor(1,0,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc  ) +
      fmajor(0,1,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc  , jtemp_loc  ) +
      fmajor(1,1,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc  ) ) +
      col_mix(1,iflav,icol,ilay) * (
        fmajor(0,0,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc+1) +
        fmajor(1,0,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc+1) +
        fmajor(0,1,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc  , jtemp_loc+1) +
        fmajor(1,1,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc+1) );

    tau(igpt,ilay,icol) += tau_major;
  });
}

// Compute minor and major species opitcal depth from pre-computed interpolation coefficients
//   (jeta,jtemp,jpress)
void compute_tau_absorption(int max_gpt_diff_lower, int max_gpt_diff_upper, int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp, int nminorlower,
                            int nminorklower, int nminorupper, int nminorkupper, int idx_h2o, int2dk const &gpoint_flavor,
                            int2dk const &band_lims_gpt, real4dk const &kmajor, real3dk const &kminor_lower, real3dk const &kminor_upper,
                            int2dk const &minor_limits_gpt_lower, int2dk const &minor_limits_gpt_upper, bool1dk const &minor_scales_with_density_lower,
                            bool1dk const &minor_scales_with_density_upper, bool1dk const &scale_by_complement_lower,
                            bool1dk const &scale_by_complement_upper, int1dk const &idx_minor_lower, int1dk const &idx_minor_upper,
                            int1dk const &idx_minor_scaling_lower, int1dk const &idx_minor_scaling_upper, int1dk const &kminor_start_lower,
                            int1dk const &kminor_start_upper, bool2dk const &tropo, real4dk const &col_mix, real6dk const &fmajor,
                            real5dk const &fminor, real2dk const &play, real2dk const &tlay, realOff3dk const &col_gas, int4dk const &jeta,
                            int2dk const &jtemp, int2dk const &jpress, real3dk const &tau, bool top_at_1) {
  int2dk itropo_lower("itropo_lower",ncol,2);
  int2dk itropo_upper("itropo_upper",ncol,2);

  int huge  = std::numeric_limits<int>::max();
  int small = std::numeric_limits<int>::min();

  if (top_at_1) {

    Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
      itropo_lower(icol,1) = nlay - 1;
      {
        int minloc = huge;
        real mn = (real) huge;
        for (int i=0; i<nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = std::min(mn,play(icol,i));
          }
        }
        itropo_lower(icol,0) = minloc;
      }

      itropo_upper(icol,0) = 0;
      {
        int maxloc = small;
        real mx = (real) small;
        for (int i=0; i<nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = std::max(mx,play(icol,i));
          }
        }
        itropo_upper(icol,1) = maxloc;
      }
    });

  } else {  // top_at_1

    Kokkos::parallel_for( ncol , KOKKOS_LAMBDA ( int icol ) {
      itropo_lower(icol,0) = 0;
      {
        int minloc = huge;
        real mn = (real) huge;
        for (int i=0; i<nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = std::min(mn,play(icol,i));
          }
        }
        itropo_lower(icol,1) = minloc;
      }

      itropo_upper(icol,1) = nlay-1;
      {
        int maxloc = small;
        real mx = (real) small;
        for (int i=0; i<nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = std::max(mx,play(icol,i));
          }
        }
        itropo_upper(icol,0) = maxloc;
      }

    });

  }  // top_at_1

  // ---------------------
  // Major Species
  // ---------------------
  gas_optical_depths_major(ncol, nlay, nbnd, ngpt, nflav, neta, npres, ntemp, gpoint_flavor,
                           band_lims_gpt, kmajor, col_mix, fmajor, jeta, tropo, jtemp, jpress, tau);
  // ---------------------
  // Minor Species - lower
  // ---------------------
  int idx_tropo = 0;
  gas_optical_depths_minor(max_gpt_diff_lower, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorlower, nminorklower,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_lower, minor_limits_gpt_lower,
                           minor_scales_with_density_lower, scale_by_complement_lower, idx_minor_lower,
                           idx_minor_scaling_lower, kminor_start_lower, play,  tlay, col_gas, fminor,
                           jeta, itropo_lower, jtemp, tau);
  // ---------------------
  // Minor Species - upper
  // ---------------------
  idx_tropo = 1;
  gas_optical_depths_minor(max_gpt_diff_upper, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorupper, nminorkupper,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_upper, minor_limits_gpt_upper,
                           minor_scales_with_density_upper, scale_by_complement_upper, idx_minor_upper,
                           idx_minor_scaling_upper, kminor_start_upper, play, tlay, col_gas, fminor,
                           jeta, itropo_upper, jtemp, tau);
}
#endif
