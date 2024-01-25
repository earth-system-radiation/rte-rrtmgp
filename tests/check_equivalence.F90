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
! ----------------------------------------------------------------------------
  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "equivalence tests stopping"
      error stop 1
    end if
  end subroutine stop_on_err
! ----------------------------------------------------------------------------------
program rte_check_equivalence
  use iso_fortran_env, only : error_unit
  !
  ! Exercise various paths through RTE+RRTMGP code that should result in the same answer
  !   Some sections test e.g. initialization and finalization 
  !   Others compare two sets of fluxes which should be the same, e.g. with respect to vertical ordering 
  !
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_util_array,     only: zero_array
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
  ! Local versions of variables
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
                            allocatable :: ref_flux_up, ref_flux_dn, ref_flux_dir, & 
                                           tst_flux_up, tst_flux_dn, tst_flux_dir, &
                                           heating_rate, jFluxUp
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_optics_rrtmgp) :: k_dist
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
  logical  :: failed

  character(len=32 ), &
            dimension(:), allocatable :: kdist_gas_names, rfmip_gas_games

  character(len=256) :: input_file = "", k_dist_file = ""
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line for any file names, block size
  !
  failed = .false. 
  nUserArgs = command_argument_count()
  if (nUserArgs <  2) call stop_on_err("Need to supply input_file k_distribution_file ")
  if (nUserArgs >  3) print *, "Ignoring command line arguments beyond the first three..."
  call get_command_argument(1,input_file)
  call get_command_argument(2,k_dist_file)
  if(trim(input_file) == '-h' .or. trim(input_file) == "--help") then
    call stop_on_err("rte_check_equivalence input_file absorption_coefficients_file")
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
    !
    ! No comparision to outside results so use absolute value of SZA
    !
    mu0(:) = cos(abs(sza(:,1)) * acos(-1._wp)/180._wp)
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
  ! Fluxes, heat rates, Jacobians
  !
  allocate(ref_flux_up(ncol,nlay+1), ref_flux_dn(ncol,nlay+1), & 
           tst_flux_up(ncol,nlay+1), tst_flux_dn(ncol,nlay+1), & 
           heating_rate(ncol, nlay)) 
  if(is_lw) then 
    allocate(jFluxUp(ncol,nlay+1))
  else
    allocate(ref_flux_dir(ncol,nlay+1), tst_flux_dir(ncol,nlay+1))
  end if 

  ! ----------------------------------------------------------------------------
  !
  ! Solvers
  !
  if(is_lw) then
    !
    ! initialization, finalization of optical properties 
    !
    call make_optical_props_1scl(k_dist)
    call atmos%finalize()
    call make_optical_props_1scl(k_dist)
    call atmos%set_name("gas only atmosphere")
    print *, "  Intialized atmosphere twice"
    !
    ! Default calculation 
    !
    fluxes%flux_up => ref_flux_up(:,:)
    fluxes%flux_dn => ref_flux_dn(:,:)
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
    print *, "  Default calculation"
    !
    ! Heating rate calculation
    !
    call stop_on_err(compute_heating_rate(ref_flux_up, ref_flux_dn, p_lev, heating_rate))
    print *, "  Computed heating rates"
    ! -------------------------------------------------------
    !
    ! Net fluxes 
    !
    nullify(fluxes%flux_up)
    nullify(fluxes%flux_dn)
    allocate(fluxes%flux_net(ncol,nlay+1))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    if(.not. allclose(fluxes%flux_net, ref_flux_dn-ref_flux_up) )  &  
      call stop_on_err("Net fluxes don't match when computed alone")
    fluxes%flux_up => tst_flux_up(:,:)
    fluxes%flux_dn => tst_flux_dn(:,:)
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    if(.not. allclose(fluxes%flux_net, ref_flux_dn-ref_flux_up) )  &  
      call report_err("Net fluxes don't match when computed in tandem")
    print *, "  Net fluxes"  
    nullify(fluxes%flux_net)
    ! -------------------------------------------------------
    !
    ! Orientation invariance 
    !
    call lw_clear_sky_vr
    if(.not. allclose(tst_flux_up, ref_flux_up, tol=4._wp) .or. &
       .not. allclose(tst_flux_dn, ref_flux_dn, tol=4._wp) )    &
      call report_err(" Vertical invariance failure")
    print *, "  Vertical orientation invariance"
    ! -------------------------------------------------------
    !
    ! Subsets of atmospheric columns 
    !
    call lw_clear_sky_subset
    if(.not. allclose(tst_flux_up, ref_flux_up) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn) )    & 
      call report_err("  Doing problem in subsets fails")
    print *, "  Subsetting invariance"
    ! -------------------------------------------------------
    !
    ! Incrementing  
    !
    atmos%tau(:,:,:) = 0.5_wp * atmos%tau(:,:,:) 
    call stop_on_err(atmos%increment(atmos))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn) )    & 
      call report_err("  halving/doubling fails")

    call increment_with_1scl
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn) )    & 
      call report_err("  Incrementing with 1scl fails")

    call increment_with_2str
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn) )    & 
      call report_err("  Incrementing with 2str fails")

    call increment_with_nstr
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn) )    & 
      call report_err("  Incrementing with nstr fails")
    print *, "  Incrementing"
    ! -------------------------------------------------------
    !
    ! Computing Jacobian shouldn't change net fluxes 
    !
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes,          &
                            flux_up_Jac = jFluxUp))
    if(.not. allclose(tst_flux_up, ref_flux_up) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn) )    & 
      call report_err("  Computing Jacobian changes fluxes")
    !
    ! Increase surface temperature by 1K and recompute fluxes 
    !
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay, sfc_t + 1._wp, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    !
    ! Comparision of fluxes with increased surface T aren't expected to match 
    !   fluxes + their Jacobian w.r.t. surface T exactly
    !
    if (.not. allclose(tst_flux_up, ref_flux_up + jFluxUp, tol=30._wp)) then
      call report_err("  Jacobian approx. differs from flux with perturbed surface T")
      print *, maxval(abs(tst_flux_up - (ref_flux_up + jFluxUp))/spacing(tst_flux_up))
    end if
    print *, "  Jacobian"
 else
    !
    ! Shortwave  
    !
    !
    ! initialization, finalization of optical properties 
    !
    call make_optical_props_2str(k_dist)
    call atmos%finalize()
    call make_optical_props_2str(k_dist)
    print *, "  Intialized atmosphere twice"

    !
    ! Default calculation 
    !
    fluxes%flux_up     => ref_flux_up (:,:)
    fluxes%flux_dn     => ref_flux_dn (:,:)
    fluxes%flux_dn_dir => ref_flux_dir(:,:)
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    print *, "  Default calculation"
    fluxes%flux_up     => tst_flux_up(:,:)
    fluxes%flux_dn     => tst_flux_dn(:,:)
    fluxes%flux_dn_dir => tst_flux_dir(:,:)

    ! -------------------------------------------------------
    !
    ! Orientation invariance 
    !
    call sw_clear_sky_vr
    if(.not. allclose(tst_flux_up, ref_flux_up, tol = 4._wp) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn, tol = 4._wp) .or. & 
       .not. allclose(tst_flux_dir,ref_flux_dir,tol = 4._wp))    &  
      call report_err(" Vertical invariance failure")
    print *, "  Vertical orientation invariance"
    ! -------------------------------------------------------
    !
    ! Scaling incoming flux
    !   Threshold of 4x spacing() works on CPUs but 8x is needed for GPUs
    !
    call sw_clear_sky_tsi
    if(.not. allclose(tst_flux_up, ref_flux_up, tol = 10._wp) .or. &
       .not. allclose(tst_flux_dn, ref_flux_dn, tol =  8._wp) .or. &
       .not. allclose(tst_flux_dir,ref_flux_dir,tol =  8._wp))     &
      call report_err("  Changing TSI fails")
    print *, "  TSI invariance"
    ! -------------------------------------------------------
    !
    ! Incrementing 
    !   Threshold of 4x spacing() works in double precision 
    !
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    atmos%tau(:,:,:) = 0.5_wp * atmos%tau(:,:,:) 
    call stop_on_err(atmos%increment(atmos))
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up, tol = 8._wp) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn, tol = 6._wp) .or. & 
       .not. allclose(tst_flux_dir,ref_flux_dir,tol = 8._wp))    &  
      call report_err("  halving/doubling fails")

    call increment_with_1scl
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up, tol = 8._wp) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn, tol = 6._wp) .or. & 
       .not. allclose(tst_flux_dir,ref_flux_dir,tol = 6._wp))    &  
      call report_err("  Incrementing with 1scl fails")

     call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
   call increment_with_2str
   if(.not. allclose(tst_flux_up, ref_flux_up, tol = 8._wp) .or. & 
      .not. allclose(tst_flux_dn, ref_flux_dn, tol = 6._wp) .or. & 
      .not. allclose(tst_flux_dir,ref_flux_dir,tol = 6._wp))    &  
      call report_err("  Incrementing with 2str fails")

    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    call increment_with_nstr
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    if(.not. allclose(tst_flux_up, ref_flux_up, tol = 8._wp) .or. & 
       .not. allclose(tst_flux_dn, ref_flux_dn, tol = 6._wp) .or. & 
       .not. allclose(tst_flux_dir,ref_flux_dir,tol = 6._wp))    &  
      call report_err("  Incrementing with nstr fails")
    print *, "  Incrementing"
  end if 
  if(failed) error stop 1
contains
  ! ----------------------------------------------------------------------------
  ! Longwave 
  ! ----------------------------------------------------------------------------
  ! Clear-sky longwave fluxes
  !   Reverse orientation in the vertical, compute, un-reverse
  !
  subroutine lw_clear_sky_vr
    type(ty_gas_concs) :: gas_concs_vr
    integer            :: i
    real(wp), dimension(ncol,nlay) :: vmr
    character(32), &
              dimension(gas_concs%get_num_gases()) &
                                  :: gc_gas_names
    ! ifort was failing this end-to-end test using the accelerator kernels
    !   The failure looks to be in the computation of optical properties. Setting
    !   do_whole_shebang = .false. vertically reverses the existing optical properties.
    !
#ifdef __INTEL_COMPILER
    logical, parameter :: do_whole_shebang = .false.
#else
    logical, parameter :: do_whole_shebang = .true.
#endif

    if (do_whole_shebang) then
        print *, "    Doing the end-to-end problem"
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
    else
      print *, "    Skipping gas optics calculation"
      atmos%tau            (:,:,:) =  atmos%tau            (:, nlay   :1:-1,:)
      lw_sources%lay_source(:,:,:) =  lw_sources%lay_source(:, nlay   :1:-1,:)
      lw_sources%lev_source(:,:,:) =  lw_sources%lev_source(:,(nlay+1):1:-1,:)
      top_at_1 = .not. top_at_1
    end if
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            sfc_emis,        &
                            fluxes))
    tst_flux_up(:,:) = tst_flux_up(:,(nlay+1):1:-1)
    tst_flux_dn(:,:) = tst_flux_dn(:,(nlay+1):1:-1)
    if (do_whole_shebang) then
      p_lay      (:,:) = p_lay      (:, nlay   :1:-1)
      t_lay      (:,:) = t_lay      (:, nlay   :1:-1)
      p_lev      (:,:) = p_lev      (:,(nlay+1):1:-1)
      t_lev      (:,:) = t_lev      (:,(nlay+1):1:-1)
    end do 
    top_at_1 = .not. top_at_1
  end subroutine lw_clear_sky_vr
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
      tst_flux_up(colS:colE,:) = up
      tst_flux_dn(colS:colE,:) = dn
    end do

  end subroutine lw_clear_sky_subset
  ! ----------------------------------------------------------------------------
  !  Shortwave 
  ! ----------------------------------------------------------------------------
  ! Shortwave - vertically reversed
  !
  subroutine sw_clear_sky_vr
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
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    !
    ! Fluxes were computed for all columns; here mask where sun is down
    !
    tst_flux_up (:,:) = tst_flux_up (:,(nlay+1):1:-1)
    tst_flux_dn (:,:) = tst_flux_dn (:,(nlay+1):1:-1)
    tst_flux_dir(:,:) = tst_flux_dir(:,(nlay+1):1:-1)

    p_lay  (:,:) = p_lay  (:, nlay   :1:-1)
    t_lay  (:,:) = t_lay  (:, nlay   :1:-1)
    p_lev  (:,:) = p_lev  (:,(nlay+1):1:-1)
    top_at_1 = .not. top_at_1
  end subroutine sw_clear_sky_vr
  ! ----------------------------------------------------------------------------
  !
  ! Shortwave - Set the total solar irradiance
  !
  subroutine sw_clear_sky_tsi
    real(wp), parameter :: tsi_scale = 0.5_wp
    real(wp)            :: default_tsi

    default_tsi = sum(toa_flux(1, :))
    ! Set TSI to half the default
    call stop_on_err(k_dist%set_tsi(tsi_scale*default_tsi))
    ! Redo gas optics
    call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    tst_flux_up (:,:) = tst_flux_up (:,:) / tsi_scale
    tst_flux_dn (:,:) = tst_flux_dn (:,:) / tsi_scale
    tst_flux_dir(:,:) = tst_flux_dir(:,:) / tsi_scale

    call stop_on_err(k_dist%set_tsi(default_tsi))
    toa_flux    (:,:) = toa_flux(:,:)     / tsi_scale

  end subroutine sw_clear_sky_tsi
  ! ----------------------------------------------------------------------------
  ! General purpose routines
  ! ----------------------------------------------------------------------------
  logical function allclose(array1, array2, tol)
    real(wp), dimension(:,:), intent(in) :: array1, array2
    real(wp), optional,       intent(in) :: tol 
    
    real(wp) :: tolerance 
    if (present(tol)) then 
      tolerance = tol 
    else
      tolerance = 2._wp
    end if 

    allclose = all(abs(array1-array2) <= tolerance * spacing(array1))
  end function allclose 
  ! ----------------------------------------------------------------------------
  subroutine report_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      failed = .true.
    end if
  end subroutine report_err
  ! ----------------------------------------------------------------------------
  subroutine increment_with_1scl 
    type(ty_optical_props_1scl) :: transparent 

    call stop_on_err(transparent%alloc_1scl(ncol, nlay, k_dist))
    call zero_array (ncol, nlay, ngpt, transparent%tau)
    call stop_on_err(transparent%increment(atmos))
    call transparent%finalize() 
  end subroutine increment_with_1scl 
  ! -------
  subroutine increment_with_2str 
    type(ty_optical_props_2str) :: transparent 

    call stop_on_err(transparent%alloc_2str(ncol, nlay, k_dist))
    call zero_array (ncol, nlay, ngpt, transparent%tau)
    call zero_array (ncol, nlay, ngpt, transparent%ssa)
    call zero_array (ncol, nlay, ngpt, transparent%g)
    call stop_on_err(transparent%increment(atmos))
    call transparent%finalize() 
  end subroutine increment_with_2str 
  ! -------
  subroutine increment_with_nstr 
   type(ty_optical_props_nstr) :: transparent 
   integer, parameter :: nmom = 4

    call stop_on_err(transparent%alloc_nstr(nmom, ncol, nlay, k_dist))
    call zero_array (      ncol, nlay, ngpt, transparent%tau)
    call zero_array (      ncol, nlay, ngpt, transparent%ssa)
    call zero_array (nmom, ncol, nlay, ngpt, transparent%p)
    call stop_on_err(transparent%increment(atmos))
    call transparent%finalize() 
  end subroutine increment_with_nstr 
 ! ----------------------------------------------------------------------------
  subroutine make_optical_props_1scl(k_dist)
    class (ty_optical_props), intent(in) :: k_dist

    if(allocated(atmos)) then
       call atmos%finalize()
       deallocate(atmos)
     end if
    allocate(ty_optical_props_1scl::atmos)
    !
    ! Allocate arrays for the optical properties .
    !
    select type(atmos)
      class is (ty_optical_props_1scl)
        call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
      class default
        call stop_on_err("rte_check_equivalence: Don't recognize the kind of optical properties ")
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
    !
    ! Allocate arrays for the optical properties themselves.
    !
    select type(atmos)
      class is (ty_optical_props_2str)
        call stop_on_err(atmos%alloc_2str(ncol, nlay, k_dist))
      class default
        call stop_on_err("rte_check_equivalence: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_2str
  ! ----------------------------------------------------------------------------
end program rte_check_equivalence
