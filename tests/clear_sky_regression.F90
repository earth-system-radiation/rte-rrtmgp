subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "clear-sky regression tests stopping"
    error stop 1
  end if
end subroutine stop_on_err
! ----------------------------------------------------------------------------------
program rte_clear_sky_regression
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, &
                                  ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty,  &
                                   read_and_block_lw_bc, read_and_block_sw_bc, determine_gas_names
  use mo_simple_netcdf,      only: get_dim_size, read_field
  use mo_heating_rates,      only: compute_heating_rate
  use mo_testing_io,         only: write_broadband_field
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
  type(ty_gas_optics_rrtmgp) :: k_dist, k_dist_2
  type(ty_gas_concs)         :: gas_concs
  type(ty_gas_concs), dimension(:), allocatable &
                             :: gas_conc_array
  class(ty_optical_props_arry), &
                 allocatable :: atmos
  type(ty_fluxes_broadband)  :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: top_at_1, is_sw, is_lw

  integer  :: ncol, nlay, nbnd, ngpt, nexp
  integer  :: icol, ilay, ibnd, iloop, igas

  integer  :: nUserArgs=0

  character(len=32 ), &
            dimension(:), allocatable :: kdist_gas_names, rfmip_gas_games

  character(len=256) :: input_file = "", k_dist_file = "", k_dist_file_2 = ""
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line for any file names, block size
  !
  nUserArgs = command_argument_count()
  if (nUserArgs <  2) call stop_on_err("Need to supply input_file k_distribution_file [k_dist_file_2]")
  if (nUserArgs >= 1) call get_command_argument(1,input_file)
  if (nUserArgs >= 2) call get_command_argument(2,k_dist_file)
  if (nUserArgs >= 3) call get_command_argument(3,k_dist_file_2)
  if (nUserArgs >  4) print *, "Ignoring command line arguments beyond the first four..."
  if(trim(input_file) == '-h' .or. trim(input_file) == "--help") then
    call stop_on_err("clear_sky_regression input_file absorption_coefficients_file")
  end if
  !
  ! Read temperature, pressure, gas concentrations.
  !   Arrays are allocated as they are read
  !
  call read_size          (input_file, ncol, nlay, nexp)
  call determine_gas_names(input_file, k_dist_file, 1, kdist_gas_names, rfmip_gas_games)
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
  call load_and_init(k_dist, k_dist_file, gas_concs)
  is_sw = k_dist%source_is_external()
  is_lw = .not. is_sw
  print *, "k-distribution is for the " // merge("longwave ", "shortwave", is_lw)
  print *, "  pressure    limits (Pa):", k_dist%get_press_min(), k_dist%get_press_max()
  print *, "  temperature limits (K):", k_dist%get_temp_min(),  k_dist%get_temp_max()
  !
  ! Problem sizes
  !
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
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
    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist))
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
  ! Solvers
  !
  fluxes%flux_up => flux_up(:,:)
  fluxes%flux_dn => flux_dn(:,:)
  if(is_lw) then
    call make_optical_props_1scl(k_dist)
    call atmos%finalize()
    call make_optical_props_1scl(k_dist)
    call atmos%set_name("gas only atmosphere")
    call lw_clear_sky_default
    call lw_clear_sky_notlev
    call lw_clear_sky_3ang
    call lw_clear_sky_optangle
    call lw_clear_sky_jaco
    call lw_clear_sky_subset
    call lw_clear_sky_vr
    call lw_clear_sky_incr
    call make_optical_props_2str(k_dist)
    call lw_clear_sky_2str
    if(len_trim(k_dist_file_2) > 0) then
      call load_and_init(k_dist_2, k_dist_file_2, gas_concs)
      print *, "Alternate k-distribution is for the " // merge("longwave ", "shortwave", .not. k_dist_2%source_is_external())
      print *, "  Resolution :", k_dist_2%get_nband(), k_dist_2%get_ngpt()
      ngpt = k_dist_2%get_ngpt()
      call atmos%finalize()
      call make_optical_props_1scl(k_dist_2)
      call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist_2))
      call lw_clear_sky_alt
    end if
  else
    call make_optical_props_2str(k_dist)
    call atmos%finalize()
    call make_optical_props_2str(k_dist)
    call sw_clear_sky_default
    call sw_clear_sky_tsi
    call sw_clear_sky_vr
    if(len_trim(k_dist_file_2) > 0) then
      call load_and_init(k_dist_2, k_dist_file_2, gas_concs)
      print *, "Alternate k-distribution is for the " // merge("longwave ", "shortwave", .not. k_dist_2%source_is_external())
      print *, "  Resolution :", k_dist_2%get_nband(), k_dist_2%get_ngpt()
      call atmos%finalize()
      call make_optical_props_2str(k_dist_2)
      deallocate(toa_flux)
      allocate(toa_flux(ncol, k_dist_2%get_ngpt()))
      call sw_clear_sky_alt
    end if
  end if
contains
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes,
  !
  subroutine lw_clear_sky_default
    real(wp), dimension(ncol, nlay+1), target :: flux_net
    real(wp), dimension(ncol, nlay)           :: heating_rate

    fluxes%flux_net => flux_net
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up,  "lw_flux_up",  "LW flux up")
    call write_broadband_field(input_file, flux_dn,  "lw_flux_dn",  "LW flux dn")
    call write_broadband_field(input_file, flux_net, "lw_flux_net", "LW flux net")

    call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
    call write_broadband_field(input_file, heating_rate,  &
                                                     "lw_flux_hr_default",  "LW heating rate", vert_dim_name = "layer")
    !
    ! Test for mo_fluxes_broadband for computing only net flux
    !
    nullify(fluxes%flux_up)
    nullify(fluxes%flux_dn)
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_net, "lw_flux_net_2", "LW flux net, direct")
    fluxes%flux_up => flux_up
    fluxes%flux_dn => flux_dn
    nullify(fluxes%flux_net)

  end subroutine lw_clear_sky_default
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, level temperatures provided
  !
  subroutine lw_clear_sky_notlev
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_notlev", "LW flux up, no level temperatures")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_notlev", "LW flux dn, no level temperatures")
  end subroutine lw_clear_sky_notlev
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !
  subroutine lw_clear_sky_3ang
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes, n_gauss_angles=3))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_3ang", "LW flux up, three quadrature angles")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_3ang", "LW flux dn, three quadrature angles")
  end subroutine lw_clear_sky_3ang
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !
  subroutine lw_clear_sky_optangle
    real(wp), dimension(ncol, ngpt) :: lw_Ds
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(k_dist%compute_optimal_angles(atmos, lw_Ds))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes, lw_Ds=lw_Ds))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_optang", "LW flux up, single optimal angles")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_optang", "LW flux dn, single optimal angles")
  end subroutine lw_clear_sky_optangle
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, computing surface Jacobian
  !
  subroutine lw_clear_sky_jaco
    real(wp), dimension(ncol,nlay+1) :: jFluxUp

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes,          &
                            flux_up_Jac = jFluxUp))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_jaco", "LW flux up, computing Jaobians")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_jaco", "LW flux dn, computing Jaobians")
    call write_broadband_field(input_file, jFluxUp, "lw_jaco_up"     , "Jacobian of LW flux up to surface temperature")

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t + 1._wp, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_stp1", "LW flux up, surface T+1K")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_stp1", "LW flux dn, surface T+1K")
  end subroutine lw_clear_sky_jaco
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, half the columns at a time
  !   We're counting on ncol being even
  !
  subroutine lw_clear_sky_subset
    type(ty_optical_props_1scl) :: atmos_subset
    type(ty_source_func_lw)     :: sources_subset
    type(ty_fluxes_broadband)   :: fluxes ! Use local variable
    real(wp), dimension(ncol/2, nlay+1), target &
                                :: up, dn
    integer :: i, colS, colE

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources, &
                                       tlev=t_lev))
    ! Testing init, then alloc for 1scl
    call stop_on_err(atmos_subset%init(atmos))
    fluxes%flux_up => up
    fluxes%flux_dn => dn
    do i = 1, 2
      colS = ((i-1) * ncol/2) + 1
      colE = i * ncol/2
      call stop_on_err(atmos%get_subset     (colS, ncol/2, atmos_subset))
      call stop_on_err(lw_sources%get_subset(colS, ncol/2, sources_subset))
      call stop_on_err(rte_lw(atmos_subset, top_at_1, &
                              sources_subset,         &
                              sfc_emis(:, colS:colE), &
                              fluxes))
      flux_up(colS:colE,:) = up
      flux_dn(colS:colE,:) = dn
    end do

    call write_broadband_field(input_file, flux_up, "lw_flux_up_subset", "LW flux up, done in two parts")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_subset", "LW flux dn, done in two parts")
  end subroutine lw_clear_sky_subset
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !   Reverse orientation in the vertical, compute, un-reverse
  !
  subroutine lw_clear_sky_vr
    type(ty_gas_concs) :: gas_concs_vr
    integer            :: i
    real(wp), dimension(ncol,nlay) :: vmr
    character(32), &
              dimension(gas_concs%get_num_gases()) &
                                  :: gc_gas_names
    !
    ! Reverse the orientation of the problem
    !
    p_lay  (:,:) = p_lay  (:, nlay   :1:-1)
    t_lay  (:,:) = t_lay  (:, nlay   :1:-1)
    p_lev  (:,:) = p_lev  (:,(nlay+1):1:-1)
    t_lev  (:,:) = t_lev  (:,(nlay+1):1:-1)
    top_at_1 = .not. top_at_1
    !
    ! No direct access to gas concentrations so use the classes
    !   This also tests otherwise uncovered routines for ty_gas_concs
    !
    gc_gas_names(:) = gas_concs%get_gas_names()
    call stop_on_err(gas_concs_vr%init(gc_gas_names(:)))
    do i = 1, gas_concs%get_num_gases()
      call stop_on_err(gas_concs%get_vmr(gc_gas_names(i), vmr))
      vmr(:,:)  = vmr(:,nlay:1:-1)
      call stop_on_err(gas_concs_vr%set_vmr(gc_gas_names(i), vmr))
    end do

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs_vr, &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    flux_up(:,:) = flux_up(:,(nlay+1):1:-1)
    flux_dn(:,:) = flux_dn(:,(nlay+1):1:-1)
    call write_broadband_field(input_file, flux_up, "lw_flux_up_vr", "LW flux up, vertically-reversed")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_vr", "LW flux dn, vertically-reversed")

    p_lay  (:,:) = p_lay  (:, nlay   :1:-1)
    t_lay  (:,:) = t_lay  (:, nlay   :1:-1)
    p_lev  (:,:) = p_lev  (:,(nlay+1):1:-1)
    t_lev  (:,:) = t_lev  (:,(nlay+1):1:-1)
    top_at_1 = .not. top_at_1
  end subroutine lw_clear_sky_vr
  ! ----------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, all info, three angles
  !   Two variations on a scattering calculation but using purely emitting/absorbing layers
  !   The first uses rescaling and should agree with _1scl answers
  !   The second uses the two-stream solver
  !
  subroutine lw_clear_sky_2str
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev=t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_1rescl", "LW flux up, clear-sky _1rescl")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_1rescl", "LW flux dn, clear-sky _1rescl")

    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes, use_2stream=.true.))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_2str", "LW flux up, clear-sky _2str")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_2str", "LW flux dn, clear-sky _2str")

  end subroutine lw_clear_sky_2str
  ! ----------------------------------------------------------------------------
  !
  ! Tests for incrementing optical properties
  !
  subroutine lw_clear_sky_incr
    type(ty_optical_props_1scl) :: atmos2
    type(ty_optical_props_2str) :: atmos_2str
    type(ty_optical_props_nstr) :: atmos_nstr
    integer, parameter          :: nmom = 4

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))

    ! Two sets of optical properties, each with half the original optical
    !   depth, should produce the same fluxes
    call stop_on_err(atmos2%alloc_1scl(ncol, nlay, atmos))
    atmos%tau (:,:,:) = atmos%tau(:,:,:) * 0.5_wp
    atmos2%tau(:,:,:) = atmos%tau(:,:,:)
    call stop_on_err(atmos2%increment(atmos))

    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_inc_1scl_with_1scl", "LW flux up, incr. 1scl with 1scl")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_inc_1scl_with_1scl", "LW flux dn, incr. 1scl with 1scl")

    !
    ! Create a transparent set of 2-stream optical properties
    !
    call stop_on_err(atmos_2str%alloc_2str(ncol, nlay, k_dist))
    atmos_2str%tau(:,:,:) = 0._wp
    atmos_2str%ssa(:,:,:) = 0._wp
    atmos_2str%g  (:,:,:) = 0._wp
    !
    ! Add a transparent 2-stream atmosphere to no-scattering
    !
    call stop_on_err(atmos_2str%increment(atmos))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_inc_1scl_with_2str", "LW flux up, incr. 1scl with 2str")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_inc_1scl_with_2str", "LW flux dn, incr. 1scl with 2str")

    !
    ! Create a transparent set of n-stream optical properties
    !
    call stop_on_err(atmos_nstr%alloc_nstr(nmom, ncol, nlay, k_dist))
    atmos_nstr%tau(:,:,:)   = 0._wp
    atmos_nstr%ssa(:,:,:)   = 0._wp
    atmos_nstr%p  (:,:,:,:) = 0._wp
    !
    ! Add a transparent n-stream atmosphere to no-scattering
    !
    call stop_on_err(atmos_nstr%increment(atmos))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up, "lw_flux_up_inc_1scl_with_nstr", "LW flux up, incr. 1scl with nstr")
    call write_broadband_field(input_file, flux_dn, "lw_flux_dn_inc_1scl_with_nstr", "LW flux dn, incr. 1scl with nstr")

  end subroutine lw_clear_sky_incr
  ! ----------------------------------------------------------------------------
  subroutine lw_clear_sky_alt
    real(wp), dimension(ncol, nlay+1), target :: flux_net
    real(wp), dimension(ncol, nlay)           :: heating_rate
    real(wp), dimension(ncol, ngpt)           :: lw_Ds

    fluxes%flux_net => flux_net
    call stop_on_err(k_dist_2%gas_optics(p_lay, p_lev, &
                                         t_lay, sfc_t, &
                                         gas_concs,    &
                                         atmos,        &
                                         lw_sources,   &
                                         tlev = t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    call write_broadband_field(input_file, flux_up,  "lw_flux_up_alt",  "LW flux up, fewer g-points")
    call write_broadband_field(input_file, flux_dn,  "lw_flux_dn_alt",  "LW flux dn, fewer g-points")
    call write_broadband_field(input_file, flux_net, "lw_flux_net_alt", "LW flux ne, fewer g-pointst")
    call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
    call write_broadband_field(input_file, heating_rate,  &
                                                     "lw_flux_hr_alt",  "LW heating rate, fewer g-points", &
                                                     vert_dim_name = "layer")

    call stop_on_err(k_dist_2%compute_optimal_angles(atmos, lw_Ds))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes, lw_Ds=lw_Ds))
    call write_broadband_field(input_file, flux_up,  "lw_flux_up_alt_oa",  "LW flux up, fewer g-points, opt. angle")
    call write_broadband_field(input_file, flux_dn,  "lw_flux_dn_alt_oa",  "LW flux dn, fewer g-points, opt. angle")
    call write_broadband_field(input_file, flux_net, "lw_flux_net_alt_oa", "LW flux ne, fewer g-points, opt. angle")
    call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))
    call write_broadband_field(input_file, heating_rate,  &
                                                     "lw_flux_hr_alt_oa",  "LW heating rate, fewer g-points, opt. angle", &
                                                     vert_dim_name = "layer")
    call k_dist_2%finalize()
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

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
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

    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    !
    ! Fluxes were computed for all columns; here mask where sun is down
    !
    where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
      flux_up  = 0._wp
      flux_dn  = 0._wp
    end where
    call write_broadband_field(input_file, flux_up,  "sw_flux_up",  "SW flux up")
    call write_broadband_field(input_file, flux_dn,  "sw_flux_dn",  "SW flux dn")

    !
    ! Two sets of optical properties, each with half the original optical
    !   depth, should produce the same fluxes
    !
    call stop_on_err(atmos2%alloc_2str(ncol, nlay, atmos))
    atmos%tau(:,:,:)  = atmos%tau(:,:,:) * 0.5_wp
    atmos2%tau(:,:,:) = atmos%tau(:,:,:)
    select type(atmos)
      class is (ty_optical_props_2str)
        atmos2%ssa(:,:,:) = atmos%ssa(:,:,:)
        atmos2%g  (:,:,:) = atmos%g  (:,:,:)
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
    call stop_on_err(atmos2%increment(atmos))
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
      flux_up  = 0._wp
      flux_dn  = 0._wp
    end where
    call write_broadband_field(input_file, flux_up,  "sw_flux_up_incr",  "SW flux up, incr. 2str with 2str")
    call write_broadband_field(input_file, flux_dn,  "sw_flux_dn_incr",  "SW flux dn, incr. 2str with 2str")
  end subroutine sw_clear_sky_default
  ! ----------------------------------------------------------------------------
  !
  ! Shortwave - Set the total solar irradiance
  !
  subroutine sw_clear_sky_tsi
    real(wp), dimension(ncol, ngpt) :: rfmip_tsi_scale
    real(wp)                        :: rrtmgp_tsi
    real(wp), parameter             :: tsi_scale = 0.5_wp

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    !
    ! Scaling factor for column dependence of TSI in RFMIP
    !
    rrtmgp_tsi = sum(toa_flux(1,:))
    rfmip_tsi_scale(:,:) = spread(tsi_3d(:,1)/rrtmgp_tsi, dim=2, ncopies=ngpt)

    ! Set TSI to half the default
    call stop_on_err(k_dist%set_tsi(tsi_scale*rrtmgp_tsi))
    ! Redo gas optics
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    toa_flux(:,:) = toa_flux(:,:) * rfmip_tsi_scale(:,:)
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    !
    ! Fluxes were computed for all columns; here mask where sun is down
    !
    where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
      flux_up  = 0._wp
      flux_dn  = 0._wp
    end where
    call write_broadband_field(input_file, flux_up /tsi_scale, "sw_flux_up_tsi",  "SW flux up, reset TSI")
    call write_broadband_field(input_file, flux_dn /tsi_scale, "sw_flux_dn_tsi",  "SW flux dn, reset TSI")

    call stop_on_err(k_dist%set_tsi(2.0_wp * sum(toa_flux(1,:))))
  end subroutine sw_clear_sky_tsi
  ! ----------------------------------------------------------------------------
  !
  ! Shortwave - vertically reversed
  !
  subroutine sw_clear_sky_vr
    real(wp), dimension(ncol, ngpt) :: rfmip_tsi_scale
    real(wp)                        :: rrtmgp_tsi
    type(ty_gas_concs) :: gas_concs_vr
    integer            :: i
    real(wp), dimension(ncol,nlay) :: vmr
    character(32), &
              dimension(gas_concs%get_num_gases()) &
                                  :: gc_gas_names

  !
  ! Reverse the orientation of the problem
  !
  p_lay  (:,:) = p_lay  (:, nlay   :1:-1)
  t_lay  (:,:) = t_lay  (:, nlay   :1:-1)
  p_lev  (:,:) = p_lev  (:,(nlay+1):1:-1)
  top_at_1 = .not. top_at_1
  !
  ! No direct access to gas concentrations so use the classes
  !   This also tests otherwise uncovered routines for ty_gas_concs
  !
  gc_gas_names(:) = gas_concs%get_gas_names()
  call stop_on_err(gas_concs_vr%init(gc_gas_names(:)))
  do i = 1, gas_concs%get_num_gases()
    call stop_on_err(gas_concs%get_vmr(gc_gas_names(i), vmr))
    vmr(:,:)  = vmr(:,nlay:1:-1)
    call stop_on_err(gas_concs_vr%set_vmr(gc_gas_names(i), vmr))
  end do

  call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                     t_lay,        &
                                     gas_concs_vr, &
                                     atmos,        &
                                     toa_flux))
  !
  ! Scaling factor for column dependence of TSI in RFMIP
  !
  rrtmgp_tsi = sum(toa_flux(1,:))
  rfmip_tsi_scale(:,:) = spread(tsi_3d(:,1)/rrtmgp_tsi, dim=2, ncopies=ngpt)
  toa_flux(:,:) = toa_flux(:,:) * rfmip_tsi_scale(:,:)

  call stop_on_err(rte_sw(atmos, top_at_1, &
                          mu0,   toa_flux, &
                          sfc_alb_dir, sfc_alb_dif, &
                          fluxes))
  !
  ! Fluxes were computed for all columns; here mask where sun is down
  !
  where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
    flux_up  = 0._wp
    flux_dn  = 0._wp
  end where
  flux_up (:,:) = flux_up (:,(nlay+1):1:-1)
  flux_dn (:,:) = flux_dn (:,(nlay+1):1:-1)
  call write_broadband_field(input_file, flux_up,  "sw_flux_up_vr",  "SW flux up, vertically-reversed")
  call write_broadband_field(input_file, flux_dn,  "sw_flux_dn_vr",  "SW flux dn, vertically-reversed")

  p_lay  (:,:) = p_lay  (:, nlay   :1:-1)
  t_lay  (:,:) = t_lay  (:, nlay   :1:-1)
  p_lev  (:,:) = p_lev  (:,(nlay+1):1:-1)
  top_at_1 = .not. top_at_1
end subroutine sw_clear_sky_vr
  ! ----------------------------------------------------------------------------
  subroutine sw_clear_sky_alt
    real(wp), dimension(ncol, k_dist_2%get_ngpt()) &
                                    :: rfmip_tsi_scale
    real(wp)                        :: rrtmgp_tsi
    call stop_on_err(k_dist_2%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    !
    ! Scaling factor for column dependence of TSI in RFMIP
    !
    rrtmgp_tsi = sum(toa_flux(1,:))
    rfmip_tsi_scale(:,:) = spread(tsi_3d(:,1)/rrtmgp_tsi, dim=2, ncopies=k_dist_2%get_ngpt())
    toa_flux(:,:) = toa_flux(:,:) * rfmip_tsi_scale(:,:)

    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    !
    ! Fluxes were computed for all columns; here mask where sun is down
    !
    where(spread(.not. sun_up, dim=2, ncopies=nlay+1))
      flux_up  = 0._wp
      flux_dn  = 0._wp
    end where
    call write_broadband_field(input_file, flux_up,  "sw_flux_up_alt",  "SW flux up, fewer g-points")
    call write_broadband_field(input_file, flux_dn,  "sw_flux_dn_alt",  "SW flux dn, fewer g-points")
  end subroutine sw_clear_sky_alt
    ! ----------------------------------------------------------------------------
  subroutine make_optical_props_1scl(k_dist)
    class (ty_optical_props), intent(in) :: k_dist

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
        call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_1scl
  ! ----------------------------------------------------------------------------
  subroutine make_optical_props_2str(k_dist)
    class (ty_optical_props), intent(in) :: k_dist
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
        call stop_on_err(atmos%alloc_2str(ncol, nlay, k_dist))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_2str
  ! ----------------------------------------------------------------------------
end program rte_clear_sky_regression
