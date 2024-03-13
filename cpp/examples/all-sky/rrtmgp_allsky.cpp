
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

#ifdef RRTMGP_ENABLE_YAKL
  yakl::init();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif

  {
#ifdef RRTMGP_ENABLE_YAKL
    using yakl::intrinsics::size;
    using yakl::intrinsics::sum;
    using yakl::intrinsics::merge;
    using yakl::intrinsics::mod;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;
#endif
#if defined(RRTMGP_ENABLE_KOKKOS) && !defined(RRTMGP_ENABLE_YAKL)
    using conv::merge;
#endif

    bool constexpr use_luts = true;

    if (argc < 5) { stoprun("Error: Fewer than 4 command line arguments provided"); }
    std::string input_file        =      argv[1];
    std::string k_dist_file       =      argv[2];
    std::string cloud_optics_file =      argv[3];
    int ncol                      = std::atoi(argv[4]);
    int nloops = 1;
    if (argc >= 6) { nloops       = std::atoi(argv[5]); }
    if (ncol   <= 0) { stoprun("Error: Number of columns must be > 0"); }
    if (nloops <= 0) { stoprun("Error: Number of loops must be > 0"); }
    if (argc > 6) { std::cout << "WARNING: Using only 5 parameters. Ignoring the rest\n"; }
    if (input_file == "-h" || input_file == "--help") {
      std::cout << "./rrtmgp_allsky  input_file  absorption_coefficients_file  cloud_optics_file  ncol  [nloops]\n\n";
      std::exit(0);
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
#ifdef RRTMGP_ENABLE_YAKL
    real2d p_lay;
    real2d t_lay;
    real2d p_lev;
    real2d t_lev;
    GasConcs gas_concs;
    real2d col_dry;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    real2dk p_lay_k;
    real2dk t_lay_k;
    real2dk p_lev_k;
    real2dk t_lev_k;
    GasConcsK gas_concs_k;
    real2dk col_dry_k;
#endif

    // Read data from the input file
    if (verbose) std::cout << "Reading input file\n\n";
#ifdef RRTMGP_ENABLE_YAKL
    read_atmos(input_file, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    read_atmos(input_file, p_lay_k, t_lay_k, p_lev_k, t_lev_k, gas_concs_k, col_dry_k, ncol);
    COMPARE_ALL_WRAP(std::vector<real2d>({p_lay, t_lay, p_lev, t_lev, col_dry}),
                     std::vector<real2dk>({p_lay_k, t_lay_k, p_lev_k, t_lev_k, col_dry_k}));
    VALIDATE_KOKKOS(gas_concs, gas_concs_k);
#endif

    int nlay = COMPUTE_SWITCH(size(p_lay,2), p_lay_k.extent(1));

    // load data into classes
    if (verbose) std::cout << "Reading k_dist file\n\n";
#ifdef RRTMGP_ENABLE_YAKL
    GasOpticsRRTMGP k_dist;
    load_and_init(k_dist, k_dist_file, gas_concs);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    GasOpticsRRTMGPK k_dist_k;
    load_and_init(k_dist_k, k_dist_file, gas_concs_k);
    VALIDATE_KOKKOS(k_dist, k_dist_k);
#endif

    bool is_sw = COMPUTE_SWITCH(k_dist.source_is_external(), k_dist_k.source_is_external());

    if (verbose) std::cout << "Reading cloud optics file\n\n";
#ifdef RRTMGP_ENABLE_YAKL
    CloudOptics cloud_optics;
    if (use_luts) {
      load_cld_lutcoeff (cloud_optics, cloud_optics_file);
    } else {
      load_cld_padecoeff(cloud_optics, cloud_optics_file);
    }
    cloud_optics.set_ice_roughness(2);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    CloudOpticsK cloud_optics_k;
    if (use_luts) {
      load_cld_lutcoeff (cloud_optics_k, cloud_optics_file);
    } else {
      load_cld_padecoeff(cloud_optics_k, cloud_optics_file);
    }
    cloud_optics_k.set_ice_roughness(2);
    VALIDATE_KOKKOS(cloud_optics, cloud_optics_k);
#endif

    // Problem sizes
    int nbnd = COMPUTE_SWITCH(k_dist.get_nband(), k_dist_k.get_nband());
    int ngpt = COMPUTE_SWITCH(k_dist.get_ngpt(), k_dist_k.get_ngpt());
#ifdef RRTMGP_ENABLE_YAKL
    auto p_lay_host = p_lay.createHostCopy();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    auto p_lay_host_k = Kokkos::create_mirror_view(p_lay_k);
    Kokkos::deep_copy(p_lay_host_k, p_lay_k);
#endif
    bool top_at_1 = COMPUTE_SWITCH(p_lay_host(1, 1) < p_lay_host(1, nlay), p_lay_host_k(0, 0) < p_lay_host_k(0, nlay-1));

    // LW calculations neglect scattering; SW calculations use the 2-stream approximation
    if (is_sw) {  // Shortwave

      if (verbose) std::cout << "This is a shortwave simulation\n\n";
#ifdef RRTMGP_ENABLE_YAKL
      OpticalProps2str atmos;
      OpticalProps2str clouds;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      OpticalProps2strK atmos_k;
      OpticalProps2strK clouds_k;
#endif

      // Clouds optical props are defined by band
#ifdef RRTMGP_ENABLE_YAKL
      clouds.init(k_dist.get_band_lims_wavenumber());
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      clouds_k.init(k_dist_k.get_band_lims_wavenumber());
      VALIDATE_KOKKOS(clouds, clouds_k);
#endif

      // Allocate arrays for the optical properties themselves.
#ifdef RRTMGP_ENABLE_YAKL
      atmos .alloc_2str(ncol, nlay, k_dist);
      clouds.alloc_2str(ncol, nlay);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      atmos_k .alloc_2str(ncol, nlay, k_dist_k);
      clouds_k.alloc_2str(ncol, nlay);
#endif

      //  Boundary conditions depending on whether the k-distribution being supplied
#ifdef RRTMGP_ENABLE_YAKL
      real2d toa_flux   ("toa_flux"   ,ncol,ngpt);
      real2d sfc_alb_dir("sfc_alb_dir",nbnd,ncol);
      real2d sfc_alb_dif("sfc_alb_dif",nbnd,ncol);
      real1d mu0        ("mu0"        ,ncol);
      // Ocean-ish values for no particular reason
      sfc_alb_dir = 0.06;
      sfc_alb_dif = 0.06;
      mu0         = 0.86;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      real2dk toa_flux_k   ("toa_flux"   ,ncol,ngpt);
      real2dk sfc_alb_dir_k("sfc_alb_dir",nbnd,ncol);
      real2dk sfc_alb_dif_k("sfc_alb_dif",nbnd,ncol);
      real1dk mu0_k        ("mu0"        ,ncol);
      // Ocean-ish values for no particular reason
      Kokkos::deep_copy(sfc_alb_dir_k, 0.06);
      Kokkos::deep_copy(sfc_alb_dif_k, 0.06);
      Kokkos::deep_copy(mu0_k        , 0.86);
#endif

      // Fluxes
#ifdef RRTMGP_ENABLE_YAKL
      real2d flux_up ("flux_up" ,ncol,nlay+1);
      real2d flux_dn ("flux_dn" ,ncol,nlay+1);
      real2d flux_dir("flux_dir",ncol,nlay+1);
      real2d flux_net("flux_net",ncol,nlay+1);
      real3d bnd_flux_up ("bnd_flux_up" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_dn ("bnd_flux_dn" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_dir("bnd_flux_dir",ncol,nlay+1,nbnd);
      real3d bnd_flux_net("bnd_flux_net",ncol,nlay+1,nbnd);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      real2dk flux_up_k ("flux_up" ,ncol,nlay+1);
      real2dk flux_dn_k ("flux_dn" ,ncol,nlay+1);
      real2dk flux_dir_k("flux_dir",ncol,nlay+1);
      real2dk flux_net_k("flux_net",ncol,nlay+1);
      real3dk bnd_flux_up_k ("bnd_flux_up" ,ncol,nlay+1,nbnd);
      real3dk bnd_flux_dn_k ("bnd_flux_dn" ,ncol,nlay+1,nbnd);
      real3dk bnd_flux_dir_k("bnd_flux_dir",ncol,nlay+1,nbnd);
      real3dk bnd_flux_net_k("bnd_flux_net",ncol,nlay+1,nbnd);
#endif

      // Clouds
#ifdef RRTMGP_ENABLE_YAKL
      real2d lwp("lwp",ncol,nlay);
      real2d iwp("iwp",ncol,nlay);
      real2d rel("rel",ncol,nlay);
      real2d rei("rei",ncol,nlay);
      bool2d cloud_mask("cloud_mask",ncol,nlay);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      real2dk lwp_k("lwp",ncol,nlay);
      real2dk iwp_k("iwp",ncol,nlay);
      real2dk rel_k("rel",ncol,nlay);
      real2dk rei_k("rei",ncol,nlay);
      bool2dk cloud_mask_k("cloud_mask",ncol,nlay);
#endif

      // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa) and not very close to the ground (< 900 hPa), and
      // put them in 2/3 of the columns since that's roughly the total cloudiness of earth
      real rel_val = COMPUTE_SWITCH(0.5 * (cloud_optics.get_min_radius_liq() + cloud_optics.get_max_radius_liq()),
                                    0.5 * (cloud_optics_k.get_min_radius_liq() + cloud_optics_k.get_max_radius_liq()));
      real rei_val = COMPUTE_SWITCH(0.5 * (cloud_optics.get_min_radius_ice() + cloud_optics.get_max_radius_ice()),
                                    0.5 * (cloud_optics_k.get_min_radius_ice() + cloud_optics_k.get_max_radius_ice()));

      // do ilay=1,nlay
      //   do icol=1,ncol
#ifdef RRTMGP_ENABLE_YAKL
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
        cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100. * 100. && p_lay(icol,ilay) < 900. * 100. && mod(icol, 3) != 0;
        // Ice and liquid will overlap in a few layers
        lwp(icol,ilay) = merge(10.,  0., cloud_mask(icol,ilay) && t_lay(icol,ilay) > 263.);
        iwp(icol,ilay) = merge(10.,  0., cloud_mask(icol,ilay) && t_lay(icol,ilay) < 273.);
        rel(icol,ilay) = merge(rel_val, 0., lwp(icol,ilay) > 0.);
        rei(icol,ilay) = merge(rei_val, 0., iwp(icol,ilay) > 0.);
      });
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}) , KOKKOS_LAMBDA (int ilay, int icol) {
        cloud_mask_k(icol,ilay) = p_lay_k(icol,ilay) > 100. * 100. && p_lay_k(icol,ilay) < 900. * 100. && ((icol+1) % 3) != 0;
        // Ice and liquid will overlap in a few layers
        lwp_k(icol,ilay) = merge(10.,  0., cloud_mask_k(icol,ilay) && t_lay_k(icol,ilay) > 263.);
        iwp_k(icol,ilay) = merge(10.,  0., cloud_mask_k(icol,ilay) && t_lay_k(icol,ilay) < 273.);
        rel_k(icol,ilay) = merge(rel_val, 0., lwp_k(icol,ilay) > 0.);
        rei_k(icol,ilay) = merge(rei_val, 0., iwp_k(icol,ilay) > 0.);
      });
      COMPARE_WRAP(cloud_mask, cloud_mask_k);
      COMPARE_ALL_WRAP(std::vector<real2d>({lwp, iwp, rel, rei}),
                       std::vector<real2dk>({lwp_k, iwp_k, rel_k, rei_k}));
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
      const size_t MaxBlockSize   = 528640000;
      const size_t MinBlockSize   = MaxBlockSize / 20;
      const size_t SuperBlockSize = MaxBlockSize;
      const size_t Capacity       = 4e9;
      conv::MemPoolSingleton::init(Capacity, MinBlockSize, MaxBlockSize, SuperBlockSize,
                                   Capacity, MinBlockSize, MaxBlockSize, SuperBlockSize);
      conv::MemPoolSingleton::print_state();
      realOff3dk col_gas  ("col_gas"     ,std::make_pair(0, ncol-1), std::make_pair(0, nlay-1), std::make_pair(-1, k_dist_k.get_ngas()-1));
#endif

      if (verbose) std::cout << "Running the main loop\n\n";
      auto start_t = std::chrono::high_resolution_clock::now();

      for (int iloop = 1 ; iloop <= nloops ; iloop++) {

#ifdef RRTMGP_ENABLE_YAKL
        cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel, rei, clouds);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        cloud_optics_k.cloud_optics(ncol, nlay, lwp_k, iwp_k, rel_k, rei_k, clouds_k);
        VALIDATE_KOKKOS(cloud_optics, cloud_optics_k);
        VALIDATE_KOKKOS(clouds, clouds_k);
#endif

        // Solvers
#ifdef RRTMGP_ENABLE_YAKL
        FluxesByband fluxes;
        fluxes.flux_up     = flux_up ;
        fluxes.flux_dn     = flux_dn ;
        fluxes.flux_dn_dir = flux_dir;
        fluxes.flux_net    = flux_net;
        fluxes.bnd_flux_up     = bnd_flux_up ;
        fluxes.bnd_flux_dn     = bnd_flux_dn ;
        fluxes.bnd_flux_dn_dir = bnd_flux_dir;
        fluxes.bnd_flux_net = bnd_flux_net;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        FluxesBybandK fluxes_k;
        fluxes_k.flux_up     = flux_up_k ;
        fluxes_k.flux_dn     = flux_dn_k ;
        fluxes_k.flux_dn_dir = flux_dir_k;
        fluxes_k.flux_net    = flux_net_k;
        fluxes_k.bnd_flux_up     = bnd_flux_up_k ;
        fluxes_k.bnd_flux_dn     = bnd_flux_dn_k ;
        fluxes_k.bnd_flux_dn_dir = bnd_flux_dir_k;
        fluxes_k.bnd_flux_net = bnd_flux_net_k;
#endif

#ifdef RRTMGP_ENABLE_YAKL
        k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay, gas_concs, atmos, toa_flux);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        k_dist_k.gas_optics(ncol, nlay, top_at_1, p_lay_k, p_lev_k, t_lay_k, gas_concs_k, col_gas, atmos_k, toa_flux_k);
        VALIDATE_KOKKOS(k_dist, k_dist_k);
        VALIDATE_KOKKOS(gas_concs, gas_concs_k);
        VALIDATE_KOKKOS(atmos, atmos_k);
        COMPARE_WRAP(toa_flux, toa_flux_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
        clouds.delta_scale();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        clouds_k.delta_scale();
        VALIDATE_KOKKOS(clouds, clouds_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
        clouds.increment(atmos);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        clouds_k.increment(atmos_k);
        VALIDATE_KOKKOS(clouds, clouds_k);
        VALIDATE_KOKKOS(atmos, atmos_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
        rte_sw(atmos, top_at_1, mu0, toa_flux, sfc_alb_dir, sfc_alb_dif, fluxes);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        rte_sw(atmos_k, top_at_1, mu0_k, toa_flux_k, sfc_alb_dir_k, sfc_alb_dif_k, fluxes_k);
        VALIDATE_KOKKOS(fluxes, fluxes_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
        if (print_norms) fluxes.print_norms();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        if (print_norms) fluxes_k.print_norms();
#endif
      }

#ifdef RRTMGP_ENABLE_KOKKOS
      conv::MemPoolSingleton::finalize();
#endif

      auto stop_t = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_t - start_t);
      std::cout << "Shortwave did " << nloops << " loops of " << ncol << " cols and " << nlay << " layers in " <<  duration.count() / 1000000.0 << " s" << std::endl;

      if (verbose) std::cout << "Writing fluxes\n\n";
#ifdef RRTMGP_ENABLE_YAKL
      if (write_fluxes) write_sw_fluxes(input_file, flux_up, flux_dn, flux_dir, ncol);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      if (write_fluxes) write_sw_fluxes(input_file, flux_up_k, flux_dn_k, flux_dir_k, ncol);
#endif

      // Hacky "unit" test against pre-computed reference fluxes
      if (ncol == 1 && nloops == 1) {
#ifdef RRTMGP_ENABLE_YAKL
        if (std::abs(sum(flux_up )-19104.862129836212)/(19104.862129836212) > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_dn )-38046.157649700355)/(38046.157649700355) > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_dir)-24998.593939345046)/(24998.593939345046) > 1.e-10) std::exit(-1);
        // And test to make sure our broadband and byband fluxes are consistent
        if (std::abs(sum(flux_up )-sum(bnd_flux_up ) )/sum(flux_up )        > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_dn )-sum(bnd_flux_dn ) )/sum(flux_dn )        > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_dir)-sum(bnd_flux_dir) )/sum(flux_dir)        > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_net)-sum(bnd_flux_net) )/sum(flux_net)        > 1.e-10) std::exit(-1);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        if (std::abs(conv::sum(flux_up_k )-19104.862129836212)/(19104.862129836212)          > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_dn_k )-38046.157649700355)/(38046.157649700355)          > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_dir_k)-24998.593939345046)/(24998.593939345046)          > 1.e-10) std::exit(-1);
        // And test to make sure our broadband and byband fluxes are consistent
        if (std::abs(conv::sum(flux_up_k )-conv::sum(bnd_flux_up_k ) )/conv::sum(flux_up_k ) > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_dn_k )-conv::sum(bnd_flux_dn_k ) )/conv::sum(flux_dn_k ) > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_dir_k)-conv::sum(bnd_flux_dir_k) )/conv::sum(flux_dir_k) > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_net_k)-conv::sum(bnd_flux_net_k) )/conv::sum(flux_net_k) > 1.e-10) std::exit(-1);
#endif
      }

    } else {  // Longwave

      if (verbose) std::cout << "This is a longwave simulation\n\n";

      // Weights and angle secants for first order (k=1) Gaussian quadrature.
      //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
      //   after Abramowitz & Stegun 1972, page 921
      int constexpr max_gauss_pts = 4;
#ifdef RRTMGP_ENABLE_YAKL
      realHost2d gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
      gauss_Ds_host(1,1) = 1.66      ; gauss_Ds_host(2,1) =         0.; gauss_Ds_host(3,1) =         0.; gauss_Ds_host(4,1) =         0.;
      gauss_Ds_host(1,2) = 1.18350343; gauss_Ds_host(2,2) = 2.81649655; gauss_Ds_host(3,2) =         0.; gauss_Ds_host(4,2) =         0.;
      gauss_Ds_host(1,3) = 1.09719858; gauss_Ds_host(2,3) = 1.69338507; gauss_Ds_host(3,3) = 4.70941630; gauss_Ds_host(4,3) =         0.;
      gauss_Ds_host(1,4) = 1.06056257; gauss_Ds_host(2,4) = 1.38282560; gauss_Ds_host(3,4) = 2.40148179; gauss_Ds_host(4,4) = 7.15513024;

      realHost2d gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
      gauss_wts_host(1,1) = 0.5         ; gauss_wts_host(2,1) = 0.          ; gauss_wts_host(3,1) = 0.          ; gauss_wts_host(4,1) = 0.          ;
      gauss_wts_host(1,2) = 0.3180413817; gauss_wts_host(2,2) = 0.1819586183; gauss_wts_host(3,2) = 0.          ; gauss_wts_host(4,2) = 0.          ;
      gauss_wts_host(1,3) = 0.2009319137; gauss_wts_host(2,3) = 0.2292411064; gauss_wts_host(3,3) = 0.0698269799; gauss_wts_host(4,3) = 0.          ;
      gauss_wts_host(1,4) = 0.1355069134; gauss_wts_host(2,4) = 0.2034645680; gauss_wts_host(3,4) = 0.1298475476; gauss_wts_host(4,4) = 0.0311809710;

      real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
      real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
      gauss_Ds_host .deep_copy_to(gauss_Ds );
      gauss_wts_host.deep_copy_to(gauss_wts);

      OpticalProps1scl atmos;
      OpticalProps1scl clouds;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      realHost2dk gauss_Ds_host_k ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
      gauss_Ds_host_k(0,0) = 1.66      ; gauss_Ds_host_k(1,0) =         0.; gauss_Ds_host_k(2,0) =         0.; gauss_Ds_host_k(3,0) =         0.;
      gauss_Ds_host_k(0,1) = 1.18350343; gauss_Ds_host_k(1,1) = 2.81649655; gauss_Ds_host_k(2,1) =         0.; gauss_Ds_host_k(3,1) =         0.;
      gauss_Ds_host_k(0,2) = 1.09719858; gauss_Ds_host_k(1,2) = 1.69338507; gauss_Ds_host_k(2,2) = 4.70941630; gauss_Ds_host_k(3,2) =         0.;
      gauss_Ds_host_k(0,3) = 1.06056257; gauss_Ds_host_k(1,3) = 1.38282560; gauss_Ds_host_k(2,3) = 2.40148179; gauss_Ds_host_k(3,3) = 7.15513024;

      realHost2dk gauss_wts_host_k("gauss_wts",max_gauss_pts,max_gauss_pts);
      gauss_wts_host_k(0,0) = 0.5         ; gauss_wts_host_k(1,0) = 0.          ; gauss_wts_host_k(2,0) = 0.          ; gauss_wts_host_k(3,0) = 0.          ;
      gauss_wts_host_k(0,1) = 0.3180413817; gauss_wts_host_k(1,1) = 0.1819586183; gauss_wts_host_k(2,1) = 0.          ; gauss_wts_host_k(3,1) = 0.          ;
      gauss_wts_host_k(0,2) = 0.2009319137; gauss_wts_host_k(1,2) = 0.2292411064; gauss_wts_host_k(2,2) = 0.0698269799; gauss_wts_host_k(3,2) = 0.          ;
      gauss_wts_host_k(0,3) = 0.1355069134; gauss_wts_host_k(1,3) = 0.2034645680; gauss_wts_host_k(2,3) = 0.1298475476; gauss_wts_host_k(3,3) = 0.0311809710;

      real2dk gauss_Ds_k ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
      real2dk gauss_wts_k("gauss_wts",max_gauss_pts,max_gauss_pts);
      Kokkos::deep_copy(gauss_Ds_k, gauss_Ds_host_k);
      Kokkos::deep_copy(gauss_wts_k, gauss_wts_host_k);
      COMPARE_WRAP(gauss_Ds, gauss_Ds_k);
      COMPARE_WRAP(gauss_wts, gauss_wts_k);

      OpticalProps1sclK atmos_k;
      OpticalProps1sclK clouds_k;
#endif

      // Clouds optical props are defined by band
#ifdef RRTMGP_ENABLE_YAKL
      clouds.init(k_dist.get_band_lims_wavenumber());
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      clouds_k.init(k_dist_k.get_band_lims_wavenumber());
      VALIDATE_KOKKOS(clouds, clouds_k);
#endif

      // Allocate arrays for the optical properties themselves.
#ifdef RRTMGP_ENABLE_YAKL
      atmos .alloc_1scl(ncol, nlay, k_dist);
      clouds.alloc_1scl(ncol, nlay);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      atmos_k .alloc_1scl(ncol, nlay, k_dist_k);
      clouds_k.alloc_1scl(ncol, nlay);
#endif

      //  Boundary conditions depending on whether the k-distribution being supplied
      //   is LW or SW
#ifdef RRTMGP_ENABLE_YAKL
      SourceFuncLW lw_sources;
      lw_sources.alloc(ncol, nlay, k_dist);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      SourceFuncLWK lw_sources_k;
      lw_sources_k.alloc(ncol, nlay, k_dist_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
      real1d t_sfc   ("t_sfc"        ,ncol);
      real2d emis_sfc("emis_sfc",nbnd,ncol);
      // Surface temperature
      auto t_lev_host = t_lev.createHostCopy();
      t_sfc    = t_lev_host(1, merge(nlay+1, 1, top_at_1));
      emis_sfc = 0.98;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      real1dk t_sfc_k   ("t_sfc"        ,ncol);
      real2dk emis_sfc_k("emis_sfc",nbnd,ncol);
      // Surface temperature
      auto t_lev_host_k = Kokkos::create_mirror_view(t_lev_k);
      Kokkos::deep_copy(t_lev_host_k, t_lev_k);
      Kokkos::deep_copy(t_sfc_k, t_lev_host_k(0, merge(nlay, 0, top_at_1)));
      Kokkos::deep_copy(emis_sfc_k, 0.98);
      COMPARE_WRAP(t_sfc, t_sfc_k);
#endif

      // Fluxes
#ifdef RRTMGP_ENABLE_YAKL
      real2d flux_up ( "flux_up" ,ncol,nlay+1);
      real2d flux_dn ( "flux_dn" ,ncol,nlay+1);
      real2d flux_net("flux_net" ,ncol,nlay+1);
      real3d bnd_flux_up ("bnd_flux_up" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_dn ("bnd_flux_dn" ,ncol,nlay+1,nbnd);
      real3d bnd_flux_net("bnd_flux_net" ,ncol,nlay+1,nbnd);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      real2dk flux_up_k ( "flux_up" ,ncol,nlay+1);
      real2dk flux_dn_k ( "flux_dn" ,ncol,nlay+1);
      real2dk flux_net_k("flux_net" ,ncol,nlay+1);
      real3dk bnd_flux_up_k ("bnd_flux_up" ,ncol,nlay+1,nbnd);
      real3dk bnd_flux_dn_k ("bnd_flux_dn" ,ncol,nlay+1,nbnd);
      real3dk bnd_flux_net_k("bnd_flux_net" ,ncol,nlay+1,nbnd);
#endif

      // Clouds
#ifdef RRTMGP_ENABLE_YAKL
      real2d lwp("lwp",ncol,nlay);
      real2d iwp("iwp",ncol,nlay);
      real2d rel("rel",ncol,nlay);
      real2d rei("rei",ncol,nlay);
      bool2d cloud_mask("cloud_mask",ncol,nlay);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      real2dk lwp_k("lwp",ncol,nlay);
      real2dk iwp_k("iwp",ncol,nlay);
      real2dk rel_k("rel",ncol,nlay);
      real2dk rei_k("rei",ncol,nlay);
      bool2dk cloud_mask_k("cloud_mask",ncol,nlay);
#endif

      // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
      //   and not very close to the ground (< 900 hPa), and
      //   put them in 2/3 of the columns since that's roughly the
      //   total cloudiness of earth
      real rel_val = COMPUTE_SWITCH(0.5 * (cloud_optics.get_min_radius_liq() + cloud_optics.get_max_radius_liq()),
                                    0.5 * (cloud_optics_k.get_min_radius_liq() + cloud_optics_k.get_max_radius_liq()));
      real rei_val = COMPUTE_SWITCH(0.5 * (cloud_optics.get_min_radius_ice() + cloud_optics.get_max_radius_ice()),
                                    0.5 * (cloud_optics_k.get_min_radius_ice() + cloud_optics_k.get_max_radius_ice()));

      // do ilay=1,nlay
      //   do icol=1,ncol
#ifdef RRTMGP_ENABLE_YAKL
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
        cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100. * 100. && p_lay(icol,ilay) < 900. * 100. && mod(icol, 3) != 0;
        // Ice and liquid will overlap in a few layers
        lwp(icol,ilay) = merge(10.,  0., cloud_mask(icol,ilay) && t_lay(icol,ilay) > 263.);
        iwp(icol,ilay) = merge(10.,  0., cloud_mask(icol,ilay) && t_lay(icol,ilay) < 273.);
        rel(icol,ilay) = merge(rel_val, 0., lwp(icol,ilay) > 0.);
        rei(icol,ilay) = merge(rei_val, 0., iwp(icol,ilay) > 0.);
      });
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}) , KOKKOS_LAMBDA (int ilay, int icol) {
        cloud_mask_k(icol,ilay) = p_lay_k(icol,ilay) > 100. * 100. && p_lay_k(icol,ilay) < 900. * 100. && ((icol+1) % 3) != 0;
        // Ice and liquid will overlap in a few layers
        lwp_k(icol,ilay) = merge(10.,  0., cloud_mask_k(icol,ilay) && t_lay_k(icol,ilay) > 263.);
        iwp_k(icol,ilay) = merge(10.,  0., cloud_mask_k(icol,ilay) && t_lay_k(icol,ilay) < 273.);
        rel_k(icol,ilay) = merge(rel_val, 0., lwp_k(icol,ilay) > 0.);
        rei_k(icol,ilay) = merge(rei_val, 0., iwp_k(icol,ilay) > 0.);
      });
      COMPARE_WRAP(cloud_mask, cloud_mask_k);
      COMPARE_WRAP(lwp, lwp_k);
      COMPARE_WRAP(iwp, iwp_k);
      COMPARE_WRAP(rel, rel_k);
      COMPARE_WRAP(rei, rei_k);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
      // const size_t MinBlockSize   = 32768;
      // const size_t MaxBlockSize   = 1048576;
      // const size_t SuperBlockSize = 1048576;
      // conv::MemPoolSingleton::init(10000000, MinBlockSize, MaxBlockSize, SuperBlockSize,
      //                              10000000, MinBlockSize, MaxBlockSize, SuperBlockSize);
      realOff3dk col_gas  ("col_gas"     ,std::make_pair(0, ncol-1), std::make_pair(0, nlay-1), std::make_pair(-1, k_dist_k.get_ngas()-1));
#endif

      // Multiple iterations for big problem sizes, and to help identify data movement
      //   For CPUs we can introduce OpenMP threading over loop iterations
      if (verbose) std::cout << "Running the main loop\n\n";
      auto start_t = std::chrono::high_resolution_clock::now();

      for (int iloop = 1 ; iloop <= nloops ; iloop++) {
#ifdef RRTMGP_ENABLE_YAKL
        cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel, rei, clouds);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        cloud_optics_k.cloud_optics(ncol, nlay, lwp_k, iwp_k, rel_k, rei_k, clouds_k);
        VALIDATE_KOKKOS(cloud_optics, cloud_optics_k);
        VALIDATE_KOKKOS(clouds, clouds_k);
#endif

        // Solvers
#ifdef RRTMGP_ENABLE_YAKL
        FluxesByband fluxes;
        fluxes.flux_up = flux_up;
        fluxes.flux_dn = flux_dn;
        fluxes.flux_net= flux_net;
        fluxes.bnd_flux_up = bnd_flux_up;
        fluxes.bnd_flux_dn = bnd_flux_dn;
        fluxes.bnd_flux_net= bnd_flux_net;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        FluxesBybandK fluxes_k;
        fluxes_k.flux_up = flux_up_k;
        fluxes_k.flux_dn = flux_dn_k;
        fluxes_k.flux_net= flux_net_k;
        fluxes_k.bnd_flux_up = bnd_flux_up_k;
        fluxes_k.bnd_flux_dn = bnd_flux_dn_k;
        fluxes_k.bnd_flux_net= bnd_flux_net_k;
#endif

        // Calling with an empty col_dry parameter
#ifdef RRTMGP_ENABLE_YAKL
        k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay, t_sfc, gas_concs, atmos, lw_sources, real2d(), t_lev);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        k_dist_k.gas_optics(ncol, nlay, top_at_1, p_lay_k, p_lev_k, t_lay_k, t_sfc_k, gas_concs_k, col_gas, atmos_k, lw_sources_k, real2dk(), t_lev_k);
        VALIDATE_KOKKOS(k_dist, k_dist_k);
        VALIDATE_KOKKOS(gas_concs, gas_concs_k);
        VALIDATE_KOKKOS(atmos, atmos_k);
        VALIDATE_KOKKOS(lw_sources, lw_sources_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
        clouds.increment(atmos);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        clouds_k.increment(atmos_k);
        VALIDATE_KOKKOS(clouds, clouds_k);
        VALIDATE_KOKKOS(atmos, atmos_k);
#endif

#ifdef RRTMGP_ENABLE_YAKL
        rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, atmos, top_at_1, lw_sources, emis_sfc, fluxes);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        rte_lw(max_gauss_pts, gauss_Ds_k, gauss_wts_k, atmos_k, top_at_1, lw_sources_k, emis_sfc_k, fluxes_k);
        VALIDATE_KOKKOS(fluxes, fluxes_k);
#endif
#ifdef RRTMGP_ENABLE_YAKL
        if (print_norms) fluxes.print_norms();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        if (print_norms) fluxes_k.print_norms();
#endif
      }

      auto stop_t = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_t - start_t);
      std::cout << "Longwave did " << nloops << " loops of " << ncol << " cols and " << nlay << " layers in " <<  duration.count() / 1000000.0 << " s" << std::endl;

      if (verbose) std::cout << "Writing fluxes\n\n";
#ifdef RRTMGP_ENABLE_YAKL
      if (write_fluxes) write_lw_fluxes(input_file, flux_up, flux_dn, ncol);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      if (write_fluxes) write_lw_fluxes(input_file, flux_up_k, flux_dn_k, ncol);
#endif

      // Hacky "unit" test against pre-computed reference fluxes
      if (ncol == 1 && nloops == 1) {
#ifdef RRTMGP_ENABLE_YAKL
        if (std::abs(sum(flux_up )-10264.518998579415)/(10264.518998579415) > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_dn )-6853.2350138542843)/(6853.2350138542843) > 1.e-10) std::exit(-1);
        // And test to make sure our broadband and byband fluxes are consistent
        if (std::abs(sum(flux_up )-sum(bnd_flux_up ) )/sum(flux_up )        > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_dn )-sum(bnd_flux_dn ) )/sum(flux_dn )        > 1.e-10) std::exit(-1);
        if (std::abs(sum(flux_net)-sum(bnd_flux_net) )/sum(flux_net)        > 1.e-10) std::exit(-1);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        if (std::abs(conv::sum(flux_up_k )-10264.518998579415)/(10264.518998579415)          > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_dn_k )-6853.2350138542843)/(6853.2350138542843)          > 1.e-10) std::exit(-1);
        // And test to make sure our broadband and byband fluxes are consistent
        if (std::abs(conv::sum(flux_up_k )-conv::sum(bnd_flux_up_k ) )/conv::sum(flux_up_k ) > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_dn_k )-conv::sum(bnd_flux_dn_k ) )/conv::sum(flux_dn_k ) > 1.e-10) std::exit(-1);
        if (std::abs(conv::sum(flux_net_k)-conv::sum(bnd_flux_net_k) )/conv::sum(flux_net_k) > 1.e-10) std::exit(-1);
#endif
      }

    }  // if (is_sw)

  }

#ifdef RRTMGP_ENABLE_YAKL
  yakl::finalize();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
  Kokkos::finalize();
#endif
}
