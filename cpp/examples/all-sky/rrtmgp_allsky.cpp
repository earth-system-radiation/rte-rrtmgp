
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "mo_gas_concentrations.h"
#include "mo_garand_atmos_io.h"
#include "mo_gas_optics_rrtmgp.h"
#include "rrtmgp_const.h"
#include "mo_load_coefficients.h"
#include "mo_load_cloud_coefficients.h"
#include "mo_fluxes.h"
#include "mo_fluxes_byband.h"
#include "mo_rte_lw.h"
#include "mo_rte_sw.h"

bool constexpr verbose      = false;
bool constexpr write_fluxes = false;
bool constexpr print_norms  = true;

int main(int argc , char **argv) {

  yakl::init();

  {
    using yakl::intrinsics::size;
    using yakl::intrinsics::sum;
    using yakl::intrinsics::merge;
    using yakl::intrinsics::mod;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    bool constexpr use_luts = true;

    if (argc < 5) { stoprun("Error: Fewer than 4 command line arguments provided"); }
    std::string input_file        =      argv[1];
    std::string k_dist_file       =      argv[2];
    std::string cloud_optics_file =      argv[3];
    int ncol                      = atoi(argv[4]);
    int nloops = 1;
    if (argc >= 6) { nloops       = atoi(argv[5]); }
    if (ncol   <= 0) { stoprun("Error: Number of columns must be > 0"); }
    if (nloops <= 0) { stoprun("Error: Number of loops must be > 0"); }
    if (argc > 6) { std::cout << "WARNING: Using only 5 parameters. Ignoring the rest\n"; }
    if (input_file == "-h" || input_file == "--help") {
      std::cout << "./rrtmgp_allsky  input_file  absorption_coefficients_file  cloud_optics_file  ncol  [nloops]\n\n";
      exit(0);
    }

    if (verbose) {
      std::cout << "Parameters: \n";
      std::cout << "    Input file:        " << input_file        << "\n";
      std::cout << "    k_dist file:       " << k_dist_file       << "\n";
      std::cout << "    Cloud Optics file: " << cloud_optics_file << "\n";
      std::cout << "    ncol:              " << ncol              << "\n";
      std::cout << "    nloops:            " << nloops            << "\n\n";
    }

    // Read temperature, pressure, gas concentrations. Arrays are allocated as they are read
    real2d p_lay;
    real2d t_lay;
    real2d p_lev;
    real2d t_lev;
    GasConcs gas_concs;
    real2d col_dry;

    // Read data from the input file
    if (verbose) std::cout << "Reading input file\n\n";
    read_atmos(input_file, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

    int nlay = size(p_lay,2);

    // load data into classes
    if (verbose) std::cout << "Reading k_dist file\n\n";
    GasOpticsRRTMGP k_dist;
    load_and_init(k_dist, k_dist_file, gas_concs);

    bool is_sw = k_dist.source_is_external();

    if (verbose) std::cout << "Reading cloud optics file\n\n";
    CloudOptics cloud_optics;
    if (use_luts) {
      load_cld_lutcoeff (cloud_optics, cloud_optics_file);
    } else {
      load_cld_padecoeff(cloud_optics, cloud_optics_file);
    }
    cloud_optics.set_ice_roughness(2);

    // Problem sizes
    int nbnd = k_dist.get_nband();
    int ngpt = k_dist.get_ngpt();
    auto p_lay_host = p_lay.createHostCopy();
    bool top_at_1 = p_lay_host(1, 1) < p_lay_host(1, nlay);


    // LW calculations neglect scattering; SW calculations use the 2-stream approximation
    if (is_sw) {  // Shortwave

      if (verbose) std::cout << "This is a shortwave simulation\n\n";
      OpticalProps2str atmos;
      OpticalProps2str clouds;

      // Clouds optical props are defined by band
      clouds.init(k_dist.get_band_lims_wavenumber());

      // Allocate arrays for the optical properties themselves.
      atmos .alloc_2str(ncol, nlay, k_dist);
      clouds.alloc_2str(ncol, nlay);

      //  Boundary conditions depending on whether the k-distribution being supplied
      real2d toa_flux   ("toa_flux"   ,ncol,ngpt);
      real2d sfc_alb_dir("sfc_alb_dir",nbnd,ncol);
      real2d sfc_alb_dif("sfc_alb_dif",nbnd,ncol);
      real1d mu0        ("mu0"        ,ncol);
      // Ocean-ish values for no particular reason
      memset( sfc_alb_dir , 0.06_wp );
      memset( sfc_alb_dif , 0.06_wp );
      memset( mu0         , 0.86_wp );

      // Fluxes
      real2d flux_up ("flux_up" ,ncol,nlay+1);
      real2d flux_dn ("flux_dn" ,ncol,nlay+1);
      real2d flux_dir("flux_dir",ncol,nlay+1);
      real2d flux_net("flux_net",ncol,nlay+1);
      real3d bnd_flux_up ("bnd_flux_up" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_dn ("bnd_flux_dn" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_dir("bnd_flux_dir",ncol,nlay+1,nbnd);
      real3d bnd_flux_net("bnd_flux_net",ncol,nlay+1,nbnd);

      // Clouds
      real2d lwp("lwp",ncol,nlay);
      real2d iwp("iwp",ncol,nlay);
      real2d rel("rel",ncol,nlay);
      real2d rei("rei",ncol,nlay);
      bool2d cloud_mask("cloud_mask",ncol,nlay);

      // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa) and not very close to the ground (< 900 hPa), and
      // put them in 2/3 of the columns since that's roughly the total cloudiness of earth
      real rel_val = 0.5 * (cloud_optics.get_min_radius_liq() + cloud_optics.get_max_radius_liq());
      real rei_val = 0.5 * (cloud_optics.get_min_radius_ice() + cloud_optics.get_max_radius_ice());
      // do ilay=1,nlay
      //   do icol=1,ncol
      parallel_for( KERNEL_NAME() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
        cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100._wp * 100._wp && p_lay(icol,ilay) < 900._wp * 100._wp && mod(icol, 3) != 0;
        // Ice and liquid will overlap in a few layers
        lwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) && t_lay(icol,ilay) > 263._wp);
        iwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) && t_lay(icol,ilay) < 273._wp);
        rel(icol,ilay) = merge(rel_val, 0._wp, lwp(icol,ilay) > 0._wp);
        rei(icol,ilay) = merge(rei_val, 0._wp, iwp(icol,ilay) > 0._wp);
      });

      if (verbose) std::cout << "Running the main loop\n\n";
      for (int iloop = 1 ; iloop <= nloops ; iloop++) {
        yakl::timer_start("shortwave");

        cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel, rei, clouds);

        // Solvers
        FluxesByband fluxes;
        fluxes.flux_up     = flux_up ;
        fluxes.flux_dn     = flux_dn ;
        fluxes.flux_dn_dir = flux_dir;
        fluxes.flux_net    = flux_net;
        fluxes.bnd_flux_up     = bnd_flux_up ;
        fluxes.bnd_flux_dn     = bnd_flux_dn ;
        fluxes.bnd_flux_dn_dir = bnd_flux_dir;
        fluxes.bnd_flux_net = bnd_flux_net;

        k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay, gas_concs, atmos, toa_flux);
        clouds.delta_scale();
        clouds.increment(atmos);
        rte_sw(atmos, top_at_1, mu0, toa_flux, sfc_alb_dir, sfc_alb_dif, fluxes);

        yakl::timer_stop("shortwave");

        if (print_norms) fluxes.print_norms();
      }

      if (verbose) std::cout << "Writing fluxes\n\n";
      if (write_fluxes) write_sw_fluxes(input_file, flux_up, flux_dn, flux_dir, ncol);

      // Hacky "unit" test against pre-computed reference fluxes
      if (ncol == 1 && nloops == 1) {
        if (abs(sum(flux_up )-19104.862129836212)/(19104.862129836212) > 1.e-10) exit(-1);
        if (abs(sum(flux_dn )-38046.157649700355)/(38046.157649700355) > 1.e-10) exit(-1);
        if (abs(sum(flux_dir)-24998.593939345046)/(24998.593939345046) > 1.e-10) exit(-1);
        // And test to make sure our broadband and byband fluxes are consistent
        if (abs(sum(flux_up )-sum(bnd_flux_up ) )/sum(flux_up )        > 1.e-10) exit(-1);
        if (abs(sum(flux_dn )-sum(bnd_flux_dn ) )/sum(flux_dn )        > 1.e-10) exit(-1);
        if (abs(sum(flux_dir)-sum(bnd_flux_dir) )/sum(flux_dir)        > 1.e-10) exit(-1);
        if (abs(sum(flux_net)-sum(bnd_flux_net) )/sum(flux_net)        > 1.e-10) exit(-1);
      }

    } else {  // Longwave

      if (verbose) std::cout << "This is a longwave simulation\n\n";

      // Weights and angle secants for first order (k=1) Gaussian quadrature.
      //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
      //   after Abramowitz & Stegun 1972, page 921
      int constexpr max_gauss_pts = 4;
      realHost2d gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
      gauss_Ds_host(1,1) = 1.66_wp      ; gauss_Ds_host(2,1) =         0._wp; gauss_Ds_host(3,1) =         0._wp; gauss_Ds_host(4,1) =         0._wp;
      gauss_Ds_host(1,2) = 1.18350343_wp; gauss_Ds_host(2,2) = 2.81649655_wp; gauss_Ds_host(3,2) =         0._wp; gauss_Ds_host(4,2) =         0._wp;
      gauss_Ds_host(1,3) = 1.09719858_wp; gauss_Ds_host(2,3) = 1.69338507_wp; gauss_Ds_host(3,3) = 4.70941630_wp; gauss_Ds_host(4,3) =         0._wp;
      gauss_Ds_host(1,4) = 1.06056257_wp; gauss_Ds_host(2,4) = 1.38282560_wp; gauss_Ds_host(3,4) = 2.40148179_wp; gauss_Ds_host(4,4) = 7.15513024_wp;

      realHost2d gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
      gauss_wts_host(1,1) = 0.5_wp         ; gauss_wts_host(2,1) = 0._wp          ; gauss_wts_host(3,1) = 0._wp          ; gauss_wts_host(4,1) = 0._wp          ;
      gauss_wts_host(1,2) = 0.3180413817_wp; gauss_wts_host(2,2) = 0.1819586183_wp; gauss_wts_host(3,2) = 0._wp          ; gauss_wts_host(4,2) = 0._wp          ;
      gauss_wts_host(1,3) = 0.2009319137_wp; gauss_wts_host(2,3) = 0.2292411064_wp; gauss_wts_host(3,3) = 0.0698269799_wp; gauss_wts_host(4,3) = 0._wp          ;
      gauss_wts_host(1,4) = 0.1355069134_wp; gauss_wts_host(2,4) = 0.2034645680_wp; gauss_wts_host(3,4) = 0.1298475476_wp; gauss_wts_host(4,4) = 0.0311809710_wp;

      real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
      real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
      gauss_Ds_host .deep_copy_to(gauss_Ds );
      gauss_wts_host.deep_copy_to(gauss_wts);

      OpticalProps1scl atmos;
      OpticalProps1scl clouds;

      // Clouds optical props are defined by band
      clouds.init(k_dist.get_band_lims_wavenumber());

      // Allocate arrays for the optical properties themselves.
      atmos .alloc_1scl(ncol, nlay, k_dist);
      clouds.alloc_1scl(ncol, nlay);

      //  Boundary conditions depending on whether the k-distribution being supplied
      //   is LW or SW
      SourceFuncLW lw_sources;
      lw_sources.alloc(ncol, nlay, k_dist);

      real1d t_sfc   ("t_sfc"        ,ncol);
      real2d emis_sfc("emis_sfc",nbnd,ncol);
      // Surface temperature
      auto t_lev_host = t_lev.createHostCopy();
      memset( t_sfc    , t_lev_host(1, merge(nlay+1, 1, top_at_1)) );
      memset( emis_sfc , 0.98_wp                                   );

      // Fluxes
      real2d flux_up ( "flux_up" ,ncol,nlay+1);
      real2d flux_dn ( "flux_dn" ,ncol,nlay+1);
      real2d flux_net("flux_net" ,ncol,nlay+1);
      real3d bnd_flux_up ("bnd_flux_up" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_dn ("bnd_flux_dn" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_net("bnd_flux_net" ,ncol,nlay+1,nbnd);

      // Clouds
      real2d lwp("lwp",ncol,nlay);
      real2d iwp("iwp",ncol,nlay);
      real2d rel("rel",ncol,nlay);
      real2d rei("rei",ncol,nlay);
      bool2d cloud_mask("cloud_mask",ncol,nlay);

      // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
      //   and not very close to the ground (< 900 hPa), and
      //   put them in 2/3 of the columns since that's roughly the
      //   total cloudiness of earth
      real rel_val = 0.5 * (cloud_optics.get_min_radius_liq() + cloud_optics.get_max_radius_liq());
      real rei_val = 0.5 * (cloud_optics.get_min_radius_ice() + cloud_optics.get_max_radius_ice());
      // do ilay=1,nlay
      //   do icol=1,ncol
      parallel_for( KERNEL_NAME() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
        cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100._wp * 100._wp && p_lay(icol,ilay) < 900._wp * 100._wp && mod(icol, 3) != 0;
        // Ice and liquid will overlap in a few layers
        lwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) && t_lay(icol,ilay) > 263._wp);
        iwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) && t_lay(icol,ilay) < 273._wp);
        rel(icol,ilay) = merge(rel_val, 0._wp, lwp(icol,ilay) > 0._wp);
        rei(icol,ilay) = merge(rei_val, 0._wp, iwp(icol,ilay) > 0._wp);
      });

      // Multiple iterations for big problem sizes, and to help identify data movement
      //   For CPUs we can introduce OpenMP threading over loop iterations
      if (verbose) std::cout << "Running the main loop\n\n";
      for (int iloop = 1 ; iloop <= nloops ; iloop++) {
        yakl::timer_start("longwave");

        cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel, rei, clouds);

        // Solvers
        FluxesByband fluxes;
        fluxes.flux_up = flux_up;
        fluxes.flux_dn = flux_dn;
        fluxes.flux_net= flux_net;
        fluxes.bnd_flux_up = bnd_flux_up;
        fluxes.bnd_flux_dn = bnd_flux_dn;
        fluxes.bnd_flux_net= bnd_flux_net;

        // Calling with an empty col_dry parameter
        k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay, t_sfc, gas_concs, atmos, lw_sources, real2d(), t_lev);
        clouds.increment(atmos);
        rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, atmos, top_at_1, lw_sources, emis_sfc, fluxes);

        yakl::timer_stop("longwave");

        if (print_norms) fluxes.print_norms();
      }

      if (verbose) std::cout << "Writing fluxes\n\n";
      if (write_fluxes) write_lw_fluxes(input_file, flux_up, flux_dn, ncol);

      // Hacky "unit" test against pre-computed reference fluxes
      if (ncol == 1 && nloops == 1) {
        if (abs(sum(flux_up )-10264.518998579415)/(10264.518998579415) > 1.e-10) exit(-1);
        if (abs(sum(flux_dn )-6853.2350138542843)/(6853.2350138542843) > 1.e-10) exit(-1);
        // And test to make sure our broadband and byband fluxes are consistent
        if (abs(sum(flux_up )-sum(bnd_flux_up ) )/sum(flux_up )        > 1.e-10) exit(-1);
        if (abs(sum(flux_dn )-sum(bnd_flux_dn ) )/sum(flux_dn )        > 1.e-10) exit(-1);
        if (abs(sum(flux_net)-sum(bnd_flux_net) )/sum(flux_net)        > 1.e-10) exit(-1);
      }

    }  // if (is_sw)

  }

  yakl::finalize();
}


