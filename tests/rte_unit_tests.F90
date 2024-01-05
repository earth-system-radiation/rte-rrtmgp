! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2023-  Atmospheric and Environmental Research,
!    Regents of the University of Colorado,
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ----------------------------------------------------------------------------
program rte_unit_tests
  !
  ! Exercise various paths through RTE code including solvers, optical properties, fluxes
  !   Tests are run on idealized problems with analytic solutions (e.g. radiative equilibrium)
  !   Solutions are checked for correctness where possible 
  !   Some tests check invariance, e.g. with respect to vertical ordering 
  !
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_util_array,     only: zero_array
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_testing_utils,      only: allclose, stop_on_err, report_err, check_fluxes, &
                                   vr, & 
                                   increment_with_1scl, increment_with_2str, increment_with_nstr
  implicit none
  ! ----------------------------------------------------------------------------------
  !
  ! Longwave tests use gray radiative equilibrium from
  !   e.g. Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770 
  !   Net flux is constant with height, OLR is known from surface temperature
  ! Tests include 
  !   Solutions match analytic results 
  !   Net fluxes = down-up when computed in various combos (net only, up/down only, all three)
  !     using ty_broadband_flues
  !   Answers are invariant to 
  !     Extracting subsets
  !     Vertical orientation
  !     Adding transparent optical properties 
  ! Longwave specific tests: 
  !   Computing the Jacibian doesn't change fluxes
  !   Fluxes inferred from Jacobian are close to fluxes with perturbed surface T. (TODO)
  !
  ! Other possibilites: 
  !   Vertical discretization? Maybe just check boundary fluxes
  !   Test the application of the boundary condition? 

  real(wp), parameter :: pi = acos(-1._wp)
  integer,  parameter :: ncol = 8, nlay = 16
  integer,  parameter :: nmu0 = 2
  integer             :: icol, ilay, imu0
  !
  ! Longwave tests - gray radiative equilibrium
  !
  real(wp), parameter :: sigma = 5.670374419e-8_wp, & ! Stefan-Boltzmann constant 
                         D     = 1.66_wp              ! Diffusivity angle, from single-angle RRTMGP solver
  real(wp), dimension(  ncol), parameter :: sfc_t     = [(285._wp, icol = 1,ncol/2), & 
                                                         (310._wp, icol = 1,ncol/2)]
  real(wp), dimension(  ncol), parameter :: lw_total_tau = [0.1_wp, 1._wp, 10._wp, 50._wp, &
                                                            0.1_wp, 1._wp, 10._wp, 50._wp] ! Would be nice to parameterize 
  real(wp), dimension(1,ncol), parameter :: sfc_emis   = 1._wp
  real(wp), dimension(ncol,1), parameter :: lw_Ds      = D     ! Diffusivity angle - use default value for all columns

  !
  ! Shorteave tests - thin atmosphere 
  !
  real(wp), dimension(2),      parameter :: g   = [0.85_wp,  0.65_wp], &
                                            tau = [1.e-4_wp, 1.e-2_wp], &
                                            ssa = 1._wp - & 
                                                  [1.e-4_wp, 1.e-2_wp]
  real(wp), dimension(nmu0),   parameter :: mu0 = [1._wp, 0.5_wp]
  real(wp), dimension(1,ncol), parameter :: sfc_albedo   = 0._wp
  real(wp), dimension(ncol,1), parameter :: toa_flux     = 1._wp

  type(ty_optical_props_1scl) :: lw_atmos 
  type(ty_source_func_lw)     :: lw_sources
  type(ty_optical_props_2str) :: sw_atmos 
  type(ty_fluxes_broadband)   :: fluxes
  logical                     :: top_at_1
  real(wp), dimension(ncol,nlay+1), target :: &
                                 ref_flux_up, ref_flux_dn, ref_flux_dir, ref_flux_net, & 
                                 tst_flux_up, tst_flux_dn, tst_flux_dir, tst_flux_net, & 
                                 jFluxUp

  logical :: passed 

  ! ------------------------------------------------------------------------------------------------------
  top_at_1 = .true. 
  ! ------------------------------------------------------------------------------------------------------
  ! 
  ! Longwave tests
  !
  ! ------------------------------------------------------------------------------------------------------
  !
  ! Gray radiative equillibrium 
  !
  print *, "Gray radiative equilibrium"
  call gray_rad_equil(sfc_t(1:ncol), lw_total_tau(1:ncol), nlay, top_at_1, lw_atmos, lw_sources)

  fluxes%flux_up  => ref_flux_up (:,:)
  fluxes%flux_dn  => ref_flux_dn (:,:)
  fluxes%flux_net => ref_flux_net(:,:)
  call stop_on_err(rte_lw(lw_atmos, top_at_1, &
                          lw_sources, &
                          sfc_emis,   &
                          fluxes))
  !
  ! Is the solution correct (does it satisfy the profile for radiative equilibrium?)
  !   Error reporting happens inside check_gray_rad_equil()
  !
  passed = check_gray_rad_equil(sfc_t, lw_total_tau, top_at_1, ref_flux_up, ref_flux_net)
  ! ------------------------------------------------------------------------------------
  !
  ! Net fluxes on- vs off-line
  !  Are the net fluxes correct?
  !
  print *, "  Longwave net flux variants"
  call check_fluxes(ref_flux_net, ref_flux_dn-ref_flux_up, passed, "net fluxes don't match down-up")
  !
  ! Compute only net fluxes 
  !
  nullify(fluxes%flux_up)
  nullify(fluxes%flux_dn)
  call stop_on_err(rte_lw(lw_atmos, top_at_1,  &
                          lw_sources, sfc_emis,&
                          fluxes))
  call check_fluxes(ref_flux_net, ref_flux_dn-ref_flux_up, &
                    passed, "Net fluxes computed alone doesn'tt match down-up computed separately")
  !
  ! Compute only up and down fluxes 
  !
  fluxes%flux_up  => tst_flux_up (:,:)
  fluxes%flux_dn  => tst_flux_dn (:,:)
  call stop_on_err(rte_lw(lw_atmos,   top_at_1, &
                          lw_sources, sfc_emis, &
                          fluxes))
  call check_fluxes(ref_flux_net, tst_flux_dn-tst_flux_up, & 
                    passed, "net fluxes don't match down-up computed together")
  ! -------------------------------------------------------
  !
  ! Subsets of atmospheric columns 
  !
  print *, "  Subsetting invariance"
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, top_at_1, lw_atmos, lw_sources)
  call lw_clear_sky_subset(lw_atmos, lw_sources, sfc_emis)
  call check_fluxes(tst_flux_up, ref_flux_up, &
                    tst_flux_dn, ref_flux_dn, &  
                    passed, "Doing problem in subsets fails")

  ! -------------------------------------------------------
  !
  ! Vertically-reverse
  !
  print *, "  Vertical orientation invariance"
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, top_at_1, lw_atmos, lw_sources)
  call vr(lw_atmos, lw_sources)
  call stop_on_err(rte_lw(lw_atmos,   .not. top_at_1, &
                          lw_sources, sfc_emis, &
                          fluxes))
  call check_fluxes(tst_flux_up(:,nlay+1:1:-1), ref_flux_up, &  
                    tst_flux_dn(:,nlay+1:1:-1), ref_flux_dn, & 
                    passed, "Doing problem upside down fails")

  ! -------------------------------------------------------
  !
  ! Computing Jacobian shouldn't change net fluxes 
  !
  print *, "  Jacobian"
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, top_at_1, lw_atmos, lw_sources)
  call stop_on_err(rte_lw(lw_atmos, top_at_1, &
                          lw_sources,      &
                          sfc_emis,        &
                          fluxes,          &
                          flux_up_Jac = jFluxUp))
  call check_fluxes(tst_flux_up, ref_flux_up, tst_flux_dn, ref_flux_dn, &  
                    passed, "Computing Jacobian changes fluxes fails")
  !
  ! Increase surface temperature in source function by 1K and recompute fluxes
  !
  lw_sources%sfc_source    (:,1) =         sigma/pi * (sfc_t + 1._wp)**4
  lw_sources%sfc_source_Jac(:,1) = 4._wp * sigma/pi * (sfc_t + 1._wp)**3
  call stop_on_err(rte_lw(lw_atmos, top_at_1, &
                          lw_sources,      &
                          sfc_emis,        &
                          fluxes))
  !
  ! Comparision of fluxes with increased surface T aren't expected to match 
  !   fluxes + their Jacobian w.r.t. surface T exactly
  !
  print '("    Jacobian accurate to within ", f6.2, "%")', & 
    maxval((tst_flux_up - ref_flux_up + jFluxUp)/tst_flux_up * 100._wp)
  ! ------------------------------------------------------------------------------------
  !
  ! Using Tang approach for purely absorbing problem should be the same 
  !
  print *, "  Two-stream optical properties"
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, top_at_1, lw_atmos, lw_sources)
  call stop_on_err(sw_atmos%alloc_2str(ncol, nlay, lw_atmos))
  sw_atmos%tau = lw_atmos%tau
  sw_atmos%ssa = 0._wp
  sw_atmos%g   = 0._wp

  call stop_on_err(rte_lw(sw_atmos, top_at_1, &
                          lw_sources,      &
                          sfc_emis,        &
                          fluxes,          &
                          flux_up_Jac = jFluxUp))
  call check_fluxes(tst_flux_up, ref_flux_up, tst_flux_dn, ref_flux_dn, &  
                    passed, "Using two-stream properties fails")
  call sw_atmos%finalize()
  ! ------------------------------------------------------------------------------------
  !
  ! Specifying diffusivity angle
  !
  print *, "  Specified transport angle"
  call stop_on_err(rte_lw(lw_atmos, top_at_1, &
                          lw_sources,      &
                          sfc_emis,        &
                          fluxes,          &
                          lw_Ds = lw_Ds))
  call check_fluxes(tst_flux_up, ref_flux_up, tst_flux_dn, ref_flux_dn, &  
                    passed, "Specifying diffusivity angle D fails")

  ! ------------------------------------------------------------------------------------
  !
  ! Shortwave tests - thin atmospheres
  !
  ! ------------------------------------------------------------------------------------
  print *, "Thin, scattering atmospheres"
  call stop_on_err(sw_atmos%alloc_2str(ncol, nlay, lw_atmos))
  call thin_scattering(tau, ssa, g, nlay, sw_atmos)

  fluxes%flux_up  => ref_flux_up (:,:)
  fluxes%flux_dn  => ref_flux_dn (:,:)
  fluxes%flux_dn_dir => ref_flux_dir(:,:)
  fluxes%flux_net => ref_flux_net(:,:)
  call stop_on_err(rte_sw(sw_atmos, top_at_1, &
                          spread(mu0(1), 1, ncol), & 
                          toa_flux,          &
                          sfc_albedo, sfc_albedo, &
                          fluxes))
  !
  ! Is the solution correct
  !   Error reporting happens inside check_gray_rad_equil()
  !
  passed = check_thin_scattering(sw_atmos, spread(mu0(1), 1, ncol), top_at_1, & 
                                  ref_flux_up, ref_flux_dn, ref_flux_dir)
  ! ------------------------------------------------------------------------------------
  !
  ! Net fluxes on- vs off-line
  !  Are the net fluxes correct?
  !
  print *, "  Shortwave net flux variants"
  call check_fluxes(ref_flux_net, ref_flux_dn-ref_flux_up, passed, "net fluxes don't match down-up")
  !
  ! Compute only net fluxes 
  !
  nullify(fluxes%flux_up)
  nullify(fluxes%flux_dn)
  call stop_on_err(rte_sw(sw_atmos, top_at_1, &
                          spread(mu0(1), 1, ncol), & 
                          toa_flux,          &
                          sfc_albedo, sfc_albedo, &
                          fluxes))
  call check_fluxes(ref_flux_net, ref_flux_dn-ref_flux_up, &
                    passed, "Net fluxes computed alone doesn't match down-up computed separately")
  !
  ! Compute only up and down fluxes 
  !
  fluxes%flux_up  => tst_flux_up (:,:)
  fluxes%flux_dn  => tst_flux_dn (:,:)
  call stop_on_err(rte_sw(sw_atmos, top_at_1, &
                          spread(mu0(1), 1, ncol), & 
                          toa_flux,          &
                          sfc_albedo, sfc_albedo, &
                          fluxes))
  call check_fluxes(ref_flux_net, tst_flux_dn-tst_flux_up, & 
                    passed, "net fluxes don't match down-up computed together")
  ! -------------------------------------------------------
  !
  ! Subsets of atmospheric columns 
  !
  print *, "  Subsetting invariance"
  call sw_clear_sky_subset(sw_atmos, sfc_albedo)
  call check_fluxes(tst_flux_up, ref_flux_up, &
                    tst_flux_dn, ref_flux_dn, &  
                    passed, "Doing problem in subsets fails")
  ! -------------------------------------------------------
  !
  ! Vertically-reverse
  !
  print *, "  Vertical orientation invariance"
  call thin_scattering(tau, ssa, g, nlay, sw_atmos)
  call vr(sw_atmos)
  call stop_on_err(rte_sw(sw_atmos, .not. top_at_1, &
                          spread(mu0(1), 1, ncol), & 
                          toa_flux,          &
                          sfc_albedo, sfc_albedo, &
                          fluxes))
  call check_fluxes(tst_flux_up(:,nlay+1:1:-1), ref_flux_up, &  
                    tst_flux_dn(:,nlay+1:1:-1), ref_flux_dn, & 
                    passed, "Doing problem upside down fails")
  call vr(sw_atmos)
  ! ------------------------------------------------------------------------------------
  ! Done
  !
  print *, "Unit tests done"
  print *
  if(.not. passed) error stop 1
  ! ------------------------------------------------------------------------------------
contains 
  ! ------------------------------------------------------------------------------------
  !
  ! Define an atmosphere in gray radiative equillibrium 
  !   See, for example, section 2 of Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770
  !
  subroutine gray_rad_equil(sfc_t, total_tau, nlay, top_at_1, atmos, sources)
    real(wp), dimension(:), intent(in) :: sfc_t, total_tau
    integer,                intent(in) :: nlay 
    logical,                intent(in) :: top_at_1
    type(ty_optical_props_1scl), intent(inout) :: atmos 
    type(ty_source_func_lw),     intent(inout) :: sources 

    integer                          :: ncol
    real(wp), dimension(size(sfc_t)) :: olr

    ncol = size(sfc_t)
    !
    ! Set up a gray spectral distribution - one band, one g-point
    !
    call stop_on_err(atmos%init(band_lims_wvn = reshape([0._wp, 3250._wp], shape = [2, 1]), & 
                                band_lims_gpt = reshape([1,     1],        shape = [2, 1]), & 
                                name = "Gray atmosphere"))
    call stop_on_err(atmos%alloc_1scl(ncol, nlay))

    !
    ! Divide optical depth evenly among layers 
    !
    atmos%tau(1:ncol,1:nlay,1) = spread(total_tau(1:ncol)/real(nlay, wp), dim=2, ncopies=nlay)

    !
    ! Longwave sources - for broadband these are sigma/pi T^4
    !   (isotropic radiation)
    !
    olr(:) = gray_rad_equil_olr(sfc_t, total_tau)

    call stop_on_err(sources%alloc(ncol, nlay, atmos))
    sources%sfc_source    (:,1) =         sigma/pi * sfc_t**4
    sources%sfc_source_Jac(:,1) = 4._wp * sigma/pi * sfc_t**3
    !
    ! Calculation with top_at_1
    !
      ilay = 1
        sources%lev_source(:,ilay,  1) = 0.5_wp/pi * olr(:)
      do ilay = 2, nlay+1
        sources%lev_source(:,ilay,  1) = 0.5_wp/pi * olr(:) * & 
                                           (1._wp + D * sum(atmos%tau(:,:ilay-1,1),dim=2))
        !
        ! The source is linear in optical depth so layer source is average of edges
        !
        sources%lay_source(:,ilay-1,1) = 0.5_wp * (sources%lev_source(:,ilay,  1) + & 
                                                   sources%lev_source(:,ilay-1,1))
      end do
    if (.not. top_at_1) then
      !
      ! Reverse vertical ordering of source functions
      !
      sources%lev_source(:,1:nlay+1,1) = sources%lev_source(:,nlay+1:1:-1,1)
      sources%lay_source(:,1:nlay,  1) = sources%lay_source(:,nlay  :1:-1,1)
    end if 
  end subroutine gray_rad_equil
  ! ------------------------------------------------------------------------------------
  !
  ! Check that solutions are in gray radiative equilibrium 
  !   We could use this to check heating rates but we'd have to make up pressure levels... 
  !
  function check_gray_rad_equil(sfc_T, lw_tau, top_at_1, up_flux, net_flux)
    real(wp), dimension(:),   intent(in) :: sfc_T, lw_tau
    real(wp), dimension(:,:), intent(in) :: up_flux, net_flux
    logical,                  intent(in) :: top_at_1
    logical                              :: check_gray_rad_equil

    logical :: passed
    integer :: toa 
    ! ------------------------------
    check_gray_rad_equil = .true. 
    toa = merge(1, size(up_flux, 2), top_at_1)

    !
    ! Check top-of-atmosphere energy balance 
    !
    if(.not. allclose(up_flux(:,toa), &
                      gray_rad_equil_olr(sfc_t, lw_tau), tol=4._wp)) then
      call report_err("OLR is not consistent with gray radiative equilibrium")
      check_gray_rad_equil = .false.
    end if 
    !
    ! Check that net fluxes are constant with height
    !  Fairly relaxed threshold w.r.t. spacing() because net flux is small relative to 
    !  large up and down fluxes that vary with tau
    !
    if(.not. allclose(net_flux(:,:), & 
                      spread(net_flux(:,1), dim=2, ncopies=size(net_flux,2)), &
                      tol = 70._wp)) then 
      call report_err("Net flux not constant with tau in gray radiative equilibrium")
      check_gray_rad_equil = .false.
    end if 
  end function check_gray_rad_equil
  ! ------------------------------------------------------------------------------------
  !
  ! Incoming energy = OLR in gray radiative equilibirum
  !   Equation 6b of Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770 with with f0 = OLR 
  !
  function gray_rad_equil_olr(T, tau)
    real(wp), dimension(:), intent(in) :: T, tau
    real(wp), dimension(size(T))       :: gray_rad_equil_olr

    gray_rad_equil_olr(:) = (2._wp * sigma * T(:)**4)/(2 + D * tau(:)) 
  end function gray_rad_equil_olr
  ! ------------------------------------------------------------------------------------
  !
  ! Define an atmosphere in gray radiative equillibrium 
  !   See, for example, section 2 of Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770
  !
  subroutine thin_scattering(tau, ssa, g, nlay, atmos)
    real(wp), dimension(:), intent(in) :: tau, ssa, g
    integer,                intent(in) :: nlay 
    type(ty_optical_props_2str), intent(inout) :: atmos 

    integer :: ntau, nssa, ng, ncol
    integer :: i, j, k
    real(wp), dimension(size(tau)*size(ssa)*size(g)) & 
            :: temp

    ntau = size(tau); nssa = size(ssa); ng = size(g)
    ncol = ntau*nssa*ng
    if(ncol /= atmos%get_ncol()) call stop_on_err("Number of SW columns incompatible")
    !
    ! Set up a gray spectral distribution - one band, one g-point
    !
    call stop_on_err(atmos%init(band_lims_wvn = reshape([3250._wp, 1.e5_wp], shape = [2, 1]), & 
                                band_lims_gpt = reshape([1,     1],          shape = [2, 1]), & 
                                name = "Gray SW atmosphere"))
    call stop_on_err(atmos%alloc_2str(ncol, nlay))

    temp = [(((tau(i), k = 1, 1      ), i = 1, ntau), j = 1, nssa*ng)]
    !
    ! Divide optical depth evenly among layers 
    !
    atmos%tau(1:ncol,1:nlay,1) = spread(temp(1:ncol)/real(nlay, wp), dim=2, ncopies=nlay)
    !
    ! ... and make the medium uniform 
    !
    temp = [(((ssa(i), k = 1, ntau   ), i = 1, nssa), j = 1, ng)]
    atmos%ssa(1:ncol,1:nlay,1) = spread(temp(1:ncol),                dim=2, ncopies=nlay)
    temp = [(((g  (i), k = 1, ntau*ng), i = 1, ng  ), j = 1, 1 )]
    atmos%g  (1:ncol,1:nlay,1) = spread(temp(1:ncol),                dim=2, ncopies=nlay)

    print *, "Original values"
    print '("tau: ", 8(e9.3,2x))', sum(atmos%tau(:,:,1),dim=2)
    print '("ssa: ", 8(e9.3,2x))',     atmos%ssa(:,1,1)
    print '("g  : ", 8(e9.3,2x))',     atmos%g  (:,1,1)
    print *
    !
    ! Delta-scale 
    !
    call stop_on_err(atmos%delta_scale())

  end subroutine thin_scattering
  ! ------------------------------------------------------------------------------------
  function check_thin_scattering(atmos, mu0, top_at_1, ref_flux_up, ref_flux_dn, ref_flux_dir)
    type(ty_optical_props_2str), intent(in) :: atmos 
    real(wp), dimension(:),      intent(in) :: mu0
    logical,                     intent(in) :: top_at_1
    real(wp), dimension(:,:),    intent(in) :: ref_flux_up, ref_flux_dn, ref_flux_dir
    logical                                 :: check_thin_scattering

    real(wp), dimension(atmos%get_ncol()) :: gamma3, R, T ! Reflectance, transmittance 
    
    check_thin_scattering = .true.
    !
    ! Solutions for small tau
    !   Meador and Weaver 1980, 10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
    !   Equation 19 using the same gamma3 as in RTE (and gamma3+gamma4=1)
    ! ssa and g are assumed to be vertically uniform
    !
    gamma3(:) = (2._wp - 3._wp * mu0(:) * atmos%g(:,1,1)) * .25_wp
    R(:) = (atmos%ssa(:,1,1)*sum(atmos%tau(:,:,1),dim=2))/mu0(:) * gamma3(:)

    print '("tau: ", 8(e9.3,2x))', sum(atmos%tau(:,:,1),dim=2)
    print '("ssa: ", 8(e9.3,2x))',     atmos%ssa(:,1,1)
    print '("g  : ", 8(e9.3,2x))',     atmos%g  (:,1,1)
    print *
    print '("R  : ", 8(e9.3,2x))', R
    print '("RTE: ", 8(e9.3,2x))', ref_flux_up(:,1)
    print '("Dif: ", 8(e9.3,2x))', R(:) - ref_flux_up(:,1)
    print '("Rel: ", 8(f9.2,2x))', (R(:) - ref_flux_up(:,1))/ref_flux_up(:,1) * 100._wp

  end function check_thin_scattering
  ! ------------------------------------------------------------------------------------
  !
  ! Invariance tests
  !
  ! ------------------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, half the columns at a time
  !   We're counting on ncol being even
  !
  subroutine lw_clear_sky_subset(lw_atmos, lw_sources, sfc_emis)
    type(ty_optical_props_1scl), intent(inout) :: lw_atmos
    type(ty_source_func_lw),     intent(inout) :: lw_sources
    real(wp), dimension(:,:),    intent(in   ) :: sfc_emis

    type(ty_optical_props_1scl) :: atmos_subset
    type(ty_source_func_lw)     :: sources_subset
    type(ty_fluxes_broadband)   :: fluxes ! Use local variable
    real(wp), dimension(lw_atmos%get_ncol()/2,       & 
                        lw_atmos%get_nlay()+1), target &
                                :: up, dn
    integer :: i, colS, colE
    integer :: ncol
    ! ------------------------------
    ncol = lw_atmos%get_ncol()
    call stop_on_err(atmos_subset%init(lw_atmos))
    fluxes%flux_up => up
    fluxes%flux_dn => dn

    do i = 1, 2
      colS = ((i-1) * ncol/2) + 1
      colE = i * ncol/2
      call stop_on_err(lw_atmos%get_subset  (colS, ncol/2, atmos_subset))
      call stop_on_err(lw_sources%get_subset(colS, ncol/2, sources_subset))
      call stop_on_err(rte_lw(atmos_subset, top_at_1,  &
                              sources_subset,          &
                              sfc_emis(:,colS:colE), &
                              fluxes))
      tst_flux_up(colS:colE,:) = up
      tst_flux_dn(colS:colE,:) = dn
    end do
  end subroutine lw_clear_sky_subset
  ! ------------------------------------------------------------------------------------
  !
  ! Clear-sky longwave fluxes, half the columns at a time
  !   We're counting on ncol being even
  !
  subroutine sw_clear_sky_subset(sw_atmos, sfc_albedo)
    type(ty_optical_props_2str), intent(inout) :: sw_atmos
    real(wp), dimension(:,:),    intent(in   ) :: sfc_albedo

    type(ty_optical_props_2str) :: atmos_subset
    type(ty_fluxes_broadband)   :: fluxes ! Use local variable
    real(wp), dimension(sw_atmos%get_ncol()/2,       & 
                        sw_atmos%get_nlay()+1), target &
                                :: up, dn
    integer :: i, colS, colE
    integer :: ncol
    ! ------------------------------
    ncol = sw_atmos%get_ncol()
    call stop_on_err(atmos_subset%init(sw_atmos))
    fluxes%flux_up => up
    fluxes%flux_dn => dn

    do i = 1, 2
      colS = ((i-1) * ncol/2) + 1
      colE = i * ncol/2
      call stop_on_err(sw_atmos%get_subset(colS, ncol/2, atmos_subset))
      call stop_on_err(rte_sw(atmos_subset, top_at_1, &
                              spread(mu0(1), 1, ncol/2), & 
                              toa_flux(colS:colE,:),     &
                              sfc_albedo(:,colS:colE), sfc_albedo(:,colS:colE), &
                              fluxes))
      tst_flux_up(colS:colE,:) = up
      tst_flux_dn(colS:colE,:) = dn
    end do
  end subroutine sw_clear_sky_subset
  ! ------------------------------------------------------------------------------------

end program rte_unit_tests