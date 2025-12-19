! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-  Atmospheric and Environmental Research,
!    Regents of the University of Colorado,
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ----------------------------------------------------------------------------------
!
! Exercise a range of alternative approaches to RFMIP clear sky problem
!   Accuracy is assessed relative to RFMIP submission with validation-plots.py
!   Serves also to exercise various code paths
! Longwave:
!   omiting level temperatures, use optimal angle, use three-angle integration,
!   two-stream solution; reduced-resolution gas optics
! Shortwave:
!   reduced-resolution gas optics
!
program rte_clear_sky_regression
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optics_utils,       only: gas_optics => gas_optics, load_and_init
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty,  &
                                   read_and_block_lw_bc, read_and_block_sw_bc, determine_gas_names
  use mo_simple_netcdf,      only: get_dim_size, read_field
  use mo_testing_utils,      only: stop_on_err
  use mo_heating_rates,      only: compute_heating_rate
  use netcdf
  use mo_simple_netcdf
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay, exp), as read from RFMIP test cases
  real(wp), dimension(:,:,:), allocatable :: p_lay_3d, t_lay_3d, p_lev_3d
  real(wp), dimension(:,:,:), allocatable :: t_lev_3d
  real(wp), dimension(:,:),   allocatable :: sfc_t_3d
  real(wp), dimension(:,:),   allocatable :: bc_3d, tsi_3d  ! boundary conditions (ncol, nexp)
  real(wp), dimension(:,:),   allocatable :: sza ! Solar zenith angle (ncol, nexp)

  !
  ! Local versions of variables - used
  !
  real(wp), dimension(:,:),     allocatable :: p_lay, t_lay, p_lev
  !
  ! Longwave only
  !
  real(wp), dimension(:,:),   allocatable :: t_lev
  real(wp), dimension(:),     allocatable :: sfc_t
  real(wp), dimension(:,:),   allocatable :: sfc_emis ! First dimension is band
  !
  ! Shortwave only
  !
  real(wp), dimension(:),    allocatable :: mu0 ! solar zenith angle, cosine of same
  real(wp), dimension(:,:),  allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  logical,  dimension(:),    allocatable :: sun_up
  !
  ! Source functions
  !
  !   Longwave
  type(ty_source_func_lw)               :: lw_sources
  !   Shortwave
  real(wp), dimension(:,:), allocatable :: toa_flux
  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_concs)         :: gas_concs
  type(ty_gas_concs), dimension(:), allocatable &
                             :: gas_conc_array
  class(ty_optical_props_arry), &
                 allocatable :: atmos
  type(ty_fluxes_broadband)  :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: is_sw, is_lw

  integer  :: ncol, nlay, nbnd, ngpt, nexp
  integer  :: icol, ilay, ibnd, iloop, igas

  integer  :: nUserArgs=0

  character(len=32 ), &
            dimension(:), allocatable :: kdist_gas_names, rfmip_gas_games

  character(len=512) :: input_file = "", output_file = "", gas_optics_file = "", gas_optics_file_2 = ""
  integer            :: ncid, dimid
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line for any file names, block size
  !
  nUserArgs = command_argument_count()
  if (nUserArgs <  3) call stop_on_err("Need to supply input_file output_file gas_optics_file [gas_optics_file_2]")
  if (nUserArgs >= 1) call get_command_argument(1,input_file)
  if (nUserArgs >= 2) call get_command_argument(2,output_file)
  if (nUserArgs >= 3) call get_command_argument(3,gas_optics_file)
  if (nUserArgs >= 4) call get_command_argument(4,gas_optics_file_2)
  if (nUserArgs >  5) print *, "Ignoring command line arguments beyond the first four..."
  if(trim(input_file) == '-h' .or. trim(input_file) == "--help") then
    call stop_on_err("clear_sky_regression input_file absorption_coefficients_file")
  end if
  !
  ! Read temperature, pressure, gas concentrations.
  !   Arrays are allocated as they are read
  !
  call read_size          (input_file, ncol, nlay, nexp)
  call determine_gas_names(input_file, gas_optics_file, 1, kdist_gas_names, rfmip_gas_games)
  call read_and_block_pt  (input_file, ncol, p_lay_3d, p_lev_3d, t_lay_3d, t_lev_3d)
  !
  ! Only do the first RFMIP experiment
  !
  allocate(p_lay(ncol, nlay), p_lev(ncol, nlay+1), &
           t_lay(ncol, nlay), t_lev(ncol, nlay+1))
  p_lay(:,:) = p_lay_3d(:,:,1)
  p_lev(:,:) = p_lev_3d(:,:,1)
  t_lay(:,:) = t_lay_3d(:,:,1)
  t_lev(:,:) = t_lev_3d(:,:,1)
  deallocate(p_lay_3d, p_lev_3d, t_lay_3d, t_lev_3d)
  !
  ! Read the gas concentrations and surface properties
  !   RFMIP I/O returns an array, we're going to use first ncol values = experiement 1 (present-day)
  !
  call read_and_block_gases_ty(input_file, ncol*nexp, kdist_gas_names, rfmip_gas_games, gas_conc_array)
  !
  ! All the profiles are in the first and only element of the array of gas concentration types
  !   Extract the first ncol profiles (this is part of testing)
  !
  call stop_on_err(gas_conc_array(1)%get_subset(1, ncol, gas_concs))
  call gas_conc_array(1)%reset()
  deallocate(gas_conc_array)
  ! ----------------------------------------------------------------------------
  ! load data into classes
  call load_and_init(gas_optics, gas_optics_file, gas_concs)
  is_sw = gas_optics%source_is_external()
  is_lw = .not. is_sw
  print *, "gas optics is for the " // merge("longwave ", "shortwave", is_lw)
  !
  ! Problem sizes
  !
  nbnd = gas_optics%get_nband()
  ngpt = gas_optics%get_ngpt()
  ! ----------------------------------------------------------------------------
  !
  !  Boundary conditions
  !
  if(is_sw) then
    allocate(toa_flux(ncol, ngpt))
    allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), &
             mu0(ncol), sun_up(ncol))
    call read_and_block_sw_bc(input_file, ncol, &
                              bc_3d, tsi_3d, sza)
    !
    ! Surface albedo is spectrally uniform
    !
    sfc_alb_dir(:,:) = spread(bc_3d(:,1), dim=1, ncopies=nbnd)
    sfc_alb_dif(:,:) = sfc_alb_dir(:,:)
    sun_up(:) = sza(:,1) < 87.5_wp ! Limits is from LBLRTM
    mu0(:) = merge(cos(sza(:,1) * acos(-1._wp)/180._wp), &
                   1._wp,                                &
                   sun_up)
  else
    allocate(sfc_t(ncol), sfc_emis(nbnd, ncol))
    call stop_on_err(lw_sources%alloc(ncol, nlay, gas_optics))
    call read_and_block_lw_bc(input_file, ncol, bc_3d, sfc_t_3d)
    !
    ! Surface emissivity is spectrally uniform
    !
    sfc_emis(:,:) = spread(bc_3d(:,1), dim=1, ncopies=nbnd)

    sfc_t   (:)   = sfc_t_3d(:,1)
  end if
  ! ----------------------------------------------------------------------------
  !
  ! Fluxes
  !
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1), flux_dir(ncol,nlay+1))
  ! ----------------------------------------------------------------------------
  !
  ! Create output file and site, level, layer dimensions
  !
  if(nf90_create(trim(output_file), NF90_CLOBBER, ncid) /= NF90_NOERR) &
    call stop_on_err("write_fluxes: can't create file " // trim(output_file))
  if(nf90_def_dim(ncid, "site", ncol, dimid) /= NF90_NOERR) &
    call stop_on_err("fail to create 'site' dimension")
  if(nf90_def_dim(ncid, "level", nlay+1, dimid) /= NF90_NOERR) &
    call stop_on_err("fail to create 'level' dimension")
  if(nf90_def_dim(ncid, "layer", nlay, dimid) /= NF90_NOERR) &
    call stop_on_err("fail to create 'layer' dimension")
  if(nf90_enddef(ncid) /= NF90_NOERR) &
    call stop_on_err("fail to to end redefinition??")

  ! ----------------------------------------------------------------------------
  !
  ! Solvers
  !
  fluxes%flux_up => flux_up(:,:)
  fluxes%flux_dn => flux_dn(:,:)
  if(is_lw) then
    call make_optical_props_1scl(gas_optics)
    call atmos%set_name("gas only atmosphere")
    call lw_clear_sky_default
    call lw_clear_sky_notlev
    call lw_clear_sky_3ang
    call lw_clear_sky_optangle
    call lw_clear_sky_jaco
    call make_optical_props_2str(gas_optics)
    call lw_clear_sky_2str
    !
    ! Replaces default gas optics with alternative
    !
    if(len_trim(gas_optics_file_2) > 0) then
      call load_and_init(gas_optics, gas_optics_file_2, gas_concs)
      print *, "Alternate gas optics is for the " // merge("longwave ", "shortwave", gas_optics%source_is_internal())
      print *, "  Resolution :", gas_optics%get_nband(), gas_optics%get_ngpt()
      ngpt = gas_optics%get_ngpt()
      call atmos%finalize()
      call make_optical_props_1scl(gas_optics)
      call stop_on_err(lw_sources%alloc(ncol, nlay, gas_optics))
      call lw_clear_sky_alt
    end if
  else
    call make_optical_props_2str(gas_optics)
    call sw_clear_sky_default
    if(len_trim(gas_optics_file_2) > 0) then
      call load_and_init(gas_optics, gas_optics_file_2, gas_concs)
      print *, "Alternate gas optics is for the " // merge("longwave ", "shortwave", gas_optics%source_is_internal())
      print *, "  Resolution :", gas_optics%get_nband(), gas_optics%get_ngpt()
      call atmos%finalize()
      call make_optical_props_2str(gas_optics)
      deallocate(toa_flux)
      allocate(toa_flux(ncol, gas_optics%get_ngpt()))
      call sw_clear_sky_alt
    end if
  end if
  ncid = nf90_close(ncid)
contains
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes,
  !
  subroutine lw_clear_sky_default
    real(wp), dimension(ncol, nlay+1), target :: flux_net
    real(wp), dimension(ncol, nlay)           :: heating_rate

    fluxes%flux_net => flux_net
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes))
    call write_broadband_field(flux_up,  "lw_flux_up",  "LW flux up")
    call write_broadband_field(flux_dn,  "lw_flux_dn",  "LW flux dn")
    call write_broadband_field(flux_net, "lw_flux_net", "LW flux net")

    call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
    call write_broadband_field(heating_rate,  &
                                                     "lw_flux_hr_default",  "LW heating rate", vert_dim_name = "layer")
    !
    ! Test for mo_fluxes_broadband for computing only net flux
    !
    nullify(fluxes%flux_up)
    nullify(fluxes%flux_dn)
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes))
    call write_broadband_field(flux_net, "lw_flux_net_2", "LW flux net, direct")
    fluxes%flux_up => flux_up
    fluxes%flux_dn => flux_dn
    nullify(fluxes%flux_net)

  end subroutine lw_clear_sky_default
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, level temperatures provided
  !
  subroutine lw_clear_sky_notlev
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes))
    call write_broadband_field(flux_up, "lw_flux_up_notlev", "LW flux up, no level temperatures")
    call write_broadband_field(flux_dn, "lw_flux_dn_notlev", "LW flux dn, no level temperatures")
  end subroutine lw_clear_sky_notlev
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !
  subroutine lw_clear_sky_3ang
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,  &
                            fluxes, n_gauss_angles=3))
    call write_broadband_field(flux_up, "lw_flux_up_3ang", "LW flux up, three quadrature angles")
    call write_broadband_field(flux_dn, "lw_flux_dn_3ang", "LW flux dn, three quadrature angles")
  end subroutine lw_clear_sky_3ang
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !
  subroutine lw_clear_sky_optangle
    real(wp), dimension(ncol, ngpt) :: lw_Ds
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(gas_optics%compute_optimal_angles(atmos, lw_Ds))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes, lw_Ds=lw_Ds))
    call write_broadband_field(flux_up, "lw_flux_up_optang", "LW flux up, single optimal angles")
    call write_broadband_field(flux_dn, "lw_flux_dn_optang", "LW flux dn, single optimal angles")
  end subroutine lw_clear_sky_optangle
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, computing surface Jacobian
  !
  subroutine lw_clear_sky_jaco
    real(wp), dimension(ncol,nlay+1) :: jFluxUp

    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes,     &
                            flux_up_Jac = jFluxUp))
    call write_broadband_field(flux_up, "lw_flux_up_jaco", "LW flux up, computing Jaobians")
    call write_broadband_field(flux_dn, "lw_flux_dn_jaco", "LW flux dn, computing Jaobians")
    call write_broadband_field(jFluxUp, "lw_jaco_up"     , "Jacobian of LW flux up to surface temperature")

    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t + 1._wp, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes))
    call write_broadband_field(flux_up, "lw_flux_up_stp1", "LW flux up, surface T+1K")
    call write_broadband_field(flux_dn, "lw_flux_dn_stp1", "LW flux dn, surface T+1K")
  end subroutine lw_clear_sky_jaco
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !   Two variations on a scattering calculation but using purely emitting/absorbing layers
  !   The first uses rescaling and should agree with _1scl answers
  !   The second uses the two-stream solver
  !
  subroutine lw_clear_sky_2str
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes))
    call write_broadband_field(flux_up, "lw_flux_up_1rescl", "LW flux up, clear-sky _1rescl")
    call write_broadband_field(flux_dn, "lw_flux_dn_1rescl", "LW flux dn, clear-sky _1rescl")

    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes, use_2stream=.true.))
    call write_broadband_field(flux_up, "lw_flux_up_2str", "LW flux up, clear-sky _2str")
    call write_broadband_field(flux_dn, "lw_flux_dn_2str", "LW flux dn, clear-sky _2str")

  end subroutine lw_clear_sky_2str
  ! ----------------------------------------------------------------------------
  subroutine lw_clear_sky_alt
    real(wp), dimension(ncol, nlay+1), target :: flux_net
    real(wp), dimension(ncol, nlay)           :: heating_rate
    real(wp), dimension(ncol, ngpt)           :: lw_Ds

    fluxes%flux_net => flux_net
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes))
    call write_broadband_field(flux_up,  "lw_flux_up_alt",  "LW flux up, fewer g-points")
    call write_broadband_field(flux_dn,  "lw_flux_dn_alt",  "LW flux dn, fewer g-points")
    call write_broadband_field(flux_net, "lw_flux_net_alt", "LW flux ne, fewer g-pointst")
    call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
    call write_broadband_field(heating_rate,  &
                                                     "lw_flux_hr_alt",  "LW heating rate, fewer g-points", &
                                                     vert_dim_name = "layer")

    call stop_on_err(gas_optics%compute_optimal_angles(atmos, lw_Ds))
    call stop_on_err(rte_lw(atmos,      &
                            lw_sources, &
                            sfc_emis,   &
                            fluxes, lw_Ds=lw_Ds))
    call write_broadband_field(flux_up,  "lw_flux_up_alt_oa",  "LW flux up, fewer g-points, opt. angle")
    call write_broadband_field(flux_dn,  "lw_flux_dn_alt_oa",  "LW flux dn, fewer g-points, opt. angle")
    call write_broadband_field(flux_net, "lw_flux_net_alt_oa", "LW flux ne, fewer g-points, opt. angle")
    call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
    call write_broadband_field(heating_rate,  &
                                                     "lw_flux_hr_alt_oa",  "LW heating rate, fewer g-points, opt. angle", &
                                                     vert_dim_name = "layer")
    call gas_optics%finalize()
  end subroutine lw_clear_sky_alt
  ! ----------------------------------------------------------------------------
  !
  ! Shortwave tests
  !
  ! ----------------------------------------------------------------------------
  !
  ! Shortwave - default
  !
  subroutine sw_clear_sky_default
    real(wp), dimension(ncol, ngpt) :: rfmip_tsi_scale
    real(wp)                        :: rrtmgp_tsi
    type(ty_optical_props_2str)     :: atmos2

    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    !
    ! Scaling factor for column dependence of TSI in RFMIP
    !
    rrtmgp_tsi = sum(toa_flux(1,:))
    rfmip_tsi_scale(:,:) = spread(tsi_3d(:,1)/rrtmgp_tsi, dim=2, ncopies=ngpt)
    toa_flux(:,:) = toa_flux(:,:) * rfmip_tsi_scale(:,:)

    call stop_on_err(rte_sw(atmos, mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    !
    ! Fluxes were computed for all columns; here mask where sun is down
    !
    where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
      flux_up  = 0._wp
      flux_dn  = 0._wp
    end where
    call write_broadband_field(flux_up,  "sw_flux_up",  "SW flux up")
    call write_broadband_field(flux_dn,  "sw_flux_dn",  "SW flux dn")
  end subroutine sw_clear_sky_default
  ! ----------------------------------------------------------------------------
  subroutine sw_clear_sky_alt
    real(wp), dimension(ncol, gas_optics%get_ngpt()) &
                                    :: rfmip_tsi_scale
    real(wp)                        :: rrtmgp_tsi
    call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    !
    ! Scaling factor for column dependence of TSI in RFMIP
    !
    rrtmgp_tsi = sum(toa_flux(1,:))
    rfmip_tsi_scale(:,:) = spread(tsi_3d(:,1)/rrtmgp_tsi, dim=2, ncopies=gas_optics%get_ngpt())
    toa_flux(:,:) = toa_flux(:,:) * rfmip_tsi_scale(:,:)

    call stop_on_err(rte_sw(atmos, mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    !
    ! Fluxes were computed for all columns; here mask where sun is down
    !
    where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
      flux_up  = 0._wp
      flux_dn  = 0._wp
    end where
    call write_broadband_field(flux_up,  "sw_flux_up_alt",  "SW flux up, fewer g-points")
    call write_broadband_field(flux_dn,  "sw_flux_dn_alt",  "SW flux dn, fewer g-points")
  end subroutine sw_clear_sky_alt
    ! ----------------------------------------------------------------------------
  subroutine make_optical_props_1scl(gas_optics)
    class (ty_optical_props), intent(in) :: gas_optics

    if(allocated(atmos)) then
       call atmos%finalize()
       deallocate(atmos)
     end if
    allocate(ty_optical_props_1scl::atmos)
    ! Clouds optical props are defined by band
    !
    ! Allocate arrays for the optical properties themselves.
    !
    select type(atmos)
      class is (ty_optical_props_1scl)
        call stop_on_err(atmos%alloc_1scl(ncol, nlay, gas_optics))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_1scl
  ! ----------------------------------------------------------------------------
  subroutine make_optical_props_2str(gas_optics)
    class (ty_optical_props), intent(in) :: gas_optics
    if(allocated(atmos)) then
       call atmos%finalize()
       deallocate(atmos)
     end if
    allocate(ty_optical_props_2str::atmos)
    ! Clouds optical props are defined by band
    !
    ! Allocate arrays for the optical properties themselves.
    !
    select type(atmos)
      class is (ty_optical_props_2str)
        call stop_on_err(atmos%alloc_2str(ncol, nlay, gas_optics))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_2str
  ! ----------------------------------------------------------------------------
  subroutine write_broadband_field(field, field_name, field_description, col_dim_name, vert_dim_name)
    !
    ! Write a field defined by column and some vertical dimension (lev or lay))
    !
    real(wp), dimension(:,:), intent(in) :: field
    character(len=*),         intent(in) :: field_name, field_description
    character(len=*), optional, &
                              intent(in) ::col_dim_name, vert_dim_name
    ! -------------------
    integer           :: varid, ncol, nlev
    !
    ! Names of column (first) and vertical (second) dimension.
    !   Because they are used in an array constuctor the need to have the same number of characters
    !
    character(len=32) :: cdim, vdim
    ! -------------------
    cdim = "site "
    vdim = "level"
    if(present(col_dim_name))  cdim = col_dim_name
    if(present(vert_dim_name)) vdim = vert_dim_name

    ncol  = size(field, dim=1)
    nlev  = size(field, dim=2)

    call create_var(ncid, trim(field_name),  [cdim, vdim], [ncol, nlev])
    call stop_on_err(write_field(ncid, trim(field_name),  field))
    !
    ! Adding descriptive text as an attribute means knowing the varid
    !
    if(nf90_inq_varid(ncid, trim(field_name), varid) /= NF90_NOERR) &
      call stop_on_err("Can't find variable " // trim(field_name))
    if(nf90_redef(ncid) /= NF90_NOERR) &
      call stop_on_err("write_broadband_field: can't put file into redefine mode")
    if(nf90_put_att(ncid, varid, "description", trim(field_description)) /= NF90_NOERR) &
      call stop_on_err("Can't write 'description' attribute to variable " // trim(field_name))
    if(nf90_enddef(ncid) /= NF90_NOERR) &
      call stop_on_err("write_broadband_field: fail to to end redefinition??")

  end subroutine write_broadband_field
  ! ----------------------------------------------------------------------------
end program
