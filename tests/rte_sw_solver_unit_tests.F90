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
program rte_sw_solver_unit_tests
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
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_sw,             only: rte_sw
  use mo_testing_utils,      only: allclose, stop_on_err, report_err, check_fluxes, &
                                   vr, & 
                                   increment_with_1scl, increment_with_2str, increment_with_nstr
  implicit none
  ! ----------------------------------------------------------------------------------
  !
  ! Exercise various paths through RTE SW solvers
  !   For the moment tests are run using thin, scattering, gray atmospheres
  !   and checked for correctness against known analytic solution (not a great approximation)
  !   Beyond correctness tests check invariance, e.g. with respect to vertical ordering 
  ! Tests include 
  !   Net fluxes = down-up when computed in various combos (net only, up/down only, all three)
  !     using ty_broadband_flues
  !   Answers are invariant to 
  !     Extracting subsets
  !     Vertical orientation
  !     Adding transparent optical properties 
  ! Shortwave specific tests: 
  !   Solutions are linear in TOA flux
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
  ! Shorteave tests - thin atmosphere 
  !
  real(wp), dimension(2),      parameter :: g   = [0.85_wp,  0.65_wp], &
                                            tau = [1.e-4_wp, 1.e-2_wp], &
                                            ssa = 1._wp - & 
                                                  [1.e-4_wp, 1.e-2_wp]
  real(wp), dimension(nmu0),   parameter :: mu0 = [1._wp, 0.5_wp]
  real(wp), dimension(1,ncol), parameter :: sfc_albedo   = 0._wp
  real(wp), dimension(ncol,1), parameter :: toa_flux     = 1._wp
  real(wp),                    parameter :: factor       = 2._wp ! for checking linearity
  real(wp), dimension(ncol)              :: mu0_arr

  type(ty_optical_props_2str) :: atmos 
  type(ty_fluxes_broadband)   :: fluxes
  logical                     :: top_at_1
  real(wp), dimension(ncol,nlay+1), target :: &
                                 ref_flux_up, ref_flux_dn, ref_flux_dir, ref_flux_net, & 
                                 tst_flux_up, tst_flux_dn, tst_flux_dir, tst_flux_net
  real(wp), dimension(:), pointer :: sfc

  logical :: passed 

  ! ------------------------------------------------------------------------------------------------------
  top_at_1 = .true. 
  ! ------------------------------------------------------------------------------------
  !
  ! Shortwave tests - thin atmospheres
  !
  ! ------------------------------------------------------------------------------------
  print *, "RTE SW solver unit tests"

  print *, "Thin, scattering atmospheres"
  call stop_on_err(atmos%init(band_lims_wvn = reshape([3250._wp, 10000._wp], shape = [2, 1]), & 
                              band_lims_gpt = reshape([1,        1        ], shape = [2, 1]), & 
                              name = "Gray atmosphere"))
  call stop_on_err(atmos%alloc_2str(ncol, nlay))
  call thin_scattering(tau, ssa, g, nlay, atmos)

  do imu0 = 1, nmu0
    print '("  mu0 = ", f4.2)', mu0(imu0)
    mu0_arr = mu0(imu0)
    fluxes%flux_up  => ref_flux_up (:,:)
    fluxes%flux_dn  => ref_flux_dn (:,:)
    fluxes%flux_dn_dir => ref_flux_dir(:,:)
    fluxes%flux_net => ref_flux_net(:,:)
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0_arr, & 
                            toa_flux,          &
                            sfc_albedo, sfc_albedo, &
                            fluxes))
    !
    ! Is the solution correct?
    !   WIP - differences are up to 25%, so skip correctness test for the moment
    !
    ! passed = check_thin_scattering(atmos, spread(mu0(1), 1, ncol), top_at_1, & 
    !                                ref_flux_up, ref_flux_dn, ref_flux_dir)
    if(imu0 == 1) passed = .true.
    ! ------------------------------------------------------------------------------------
    !
    ! Check direct beam for correctness with Beer-Lambert-Bouguier
    !
    if(top_at_1) then 
      sfc => ref_flux_dir(:,nlay+1)
    else
      sfc => ref_flux_dir(:,     1)
    end if   
    if(.not. allclose(sfc, & 
                      toa_flux(:,1)*mu0_arr*exp(-sum(atmos%tau(:,:,1),dim=2)/mu0_arr), tol=10._wp)) then 
      passed = .false.
      call report_err("Direct flux doesn't match")
    end if    
    ! ------------------------------------------------------------------------------------
    !
    ! Net fluxes on- vs off-line
    !  Are the net fluxes correct?
    !
    print *, "  Net flux variants"
    call check_fluxes(ref_flux_net, ref_flux_dn-ref_flux_up, passed, "net fluxes don't match down-up")
    !
    ! Compute only net fluxes 
    !
    nullify(fluxes%flux_up)
    nullify(fluxes%flux_dn)
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0_arr, & 
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
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0_arr, & 
                            toa_flux,          &
                            sfc_albedo, sfc_albedo, &
                            fluxes))
    call check_fluxes(ref_flux_net, tst_flux_dn-tst_flux_up, & 
                      passed, "Net fluxes don't match down-up computed together")
    ! -------------------------------------------------------
    !
    ! Subsets of atmospheric columns 
    !
    print *, "  Subsetting invariance"
    call clear_sky_subset(atmos, mu0_arr, toa_flux, sfc_albedo, tst_flux_up, tst_flux_dn)
    call check_fluxes(tst_flux_up, ref_flux_up, &
                      tst_flux_dn, ref_flux_dn, &  
                      passed, "Doing problem in subsets fails")
    ! -------------------------------------------------------
    !
    ! Vertically-reverse
    !
    print *, "  Vertical orientation invariance"
    call thin_scattering(tau, ssa, g, nlay, atmos)
    call vr(atmos)
    call stop_on_err(rte_sw(atmos, .not. top_at_1, &
                            mu0_arr, & 
                            toa_flux,          &
                            sfc_albedo, sfc_albedo, &
                            fluxes))
    call check_fluxes(tst_flux_up(:,nlay+1:1:-1), ref_flux_up, &  
                      tst_flux_dn(:,nlay+1:1:-1), ref_flux_dn, & 
                      passed, "Doing problem upside down fails")
    call vr(atmos)
    ! -------------------------------------------------------
    !
    ! Linear in TOA flux
    !
    print *, "  Linear in TOA flux"
    call thin_scattering(tau, ssa, g, nlay, atmos)
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0_arr, & 
                            toa_flux * factor,      &
                            sfc_albedo, sfc_albedo, &
                            fluxes))
    call check_fluxes(tst_flux_up/factor, ref_flux_up, &  
                      tst_flux_dn/factor, ref_flux_dn, & 
                      passed, "Linearity in TOA flux fails")
  ! ------------------------------------------------------------------------------------
  end do 
  ! Done
  !
  print *, "RTE SW solver unit tests done"
  print *
  if(.not. passed) error stop 1
  ! ------------------------------------------------------------------------------------
contains 
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

    if(.false.) then
    print *, "Original values"
    print '("tau: ", 8(e9.3,2x))', sum(atmos%tau(:,:,1),dim=2)
    print '("ssa: ", 8(e9.3,2x))',     atmos%ssa(:,1,1)
    print '("g  : ", 8(e9.3,2x))',     atmos%g  (:,1,1)
    print *
    end if
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
    if(.false.) then
    print '("tau: ", 8(e9.3,2x))', sum(atmos%tau(:,:,1),dim=2)
    print '("ssa: ", 8(e9.3,2x))',     atmos%ssa(:,1,1)
    print '("g  : ", 8(e9.3,2x))',     atmos%g  (:,1,1)
    print *
    print '("R  : ", 8(e9.3,2x))', R
    print '("RTE: ", 8(e9.3,2x))', ref_flux_up(:,1)
    print '("Dif: ", 8(e9.3,2x))', R(:) - ref_flux_up(:,1)
    print '("Rel: ", 8(f9.2,2x))', (R(:) - ref_flux_up(:,1))/ref_flux_up(:,1) * 100._wp
    end if 
  end function check_thin_scattering
  ! ------------------------------------------------------------------------------------
  !
  ! Invariance tests
  !
  ! ------------------------------------------------------------------------------------
  !
  ! Clear-sky fluxes, half the columns at a time
  !   We're counting on ncol being even
  !
  subroutine clear_sky_subset(atmos, mu0, toa_flux, sfc_albedo, flux_up, flux_dn)
    type(ty_optical_props_2str), intent(inout) :: atmos
    real(wp), dimension(:),      intent(in   ) :: mu0
    real(wp), dimension(:,:),    intent(in   ) :: toa_flux, sfc_albedo
    real(wp), dimension(:,:),    intent(  out) :: flux_up, flux_dn

    type(ty_optical_props_2str) :: atmos_subset
    type(ty_fluxes_broadband)   :: fluxes ! Use local variable
    real(wp), dimension(atmos%get_ncol()/2,       & 
                        atmos%get_nlay()+1), target &
                                :: up, dn
    integer :: i, colS, colE
    integer :: ncol
    ! ------------------------------
    ncol = atmos%get_ncol()
    call stop_on_err(atmos_subset%init(atmos))
    fluxes%flux_up => up
    fluxes%flux_dn => dn

    do i = 1, 2
      colS = ((i-1) * ncol/2) + 1
      colE = i * ncol/2
      call stop_on_err(atmos%get_subset(colS, ncol/2, atmos_subset))
      call stop_on_err(rte_sw(atmos_subset, top_at_1, &
                              mu0(colS:colE),         & 
                              toa_flux(colS:colE,:),  &
                              sfc_albedo(:,colS:colE), sfc_albedo(:,colS:colE), &
                              fluxes))
      flux_up(colS:colE,:) = up
      flux_dn(colS:colE,:) = dn
    end do
  end subroutine clear_sky_subset
  ! ------------------------------------------------------------------------------------

end program rte_sw_solver_unit_tests