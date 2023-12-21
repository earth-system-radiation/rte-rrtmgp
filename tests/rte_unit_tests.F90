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
  !   Some sections test e.g. initialization and finalization 
  !   Others compare two sets of fluxes which should be the same, e.g. with respect to vertical ordering 
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
  use mo_testing_utils,      only: allclose, stop_on_err, report_err, in_rad_eq
  implicit none
  ! ----------------------------------------------------------------------------------
  !
  ! Do we vary the number of columns? Layers? 
  !   For gray radiative equilibrium: 
  !   We want 
  !     several optical thickness values and several surface temperatures - let's say 8 in total 
  !     even and uneven discretizations in tau?
  !   Want to check correctness
  !     OLR, surface T, and optical depth have the relationships we expect 
  !     net fluxes are constant with height 
  !     net fluxes on- and off-line are the same (this tests ty_fluxes_broadband)
  !     heating rates? couldn't check correctness 
  !   Then use an example to test invariance to vertical orientation 
  !     subsetting
  !     maybe not to incrementing since we're restrictied to a single solver 
  !     vertical discretization? Maybe just check boundary fluxes
  !     *Jacobian
  !     *Computing Jacobian shouldn't change net fluxes 
  !   Expand testing to optical properties 1_scl? 
  !     Increment with 3 variants of transparent media
  !     Divide by two and increment 
  !   Test the application of the boundary condition? 

  real(wp), parameter :: pi = acos(-1._wp)
  integer,  parameter :: ncol = 8, nlay = 10
  integer             :: icol, ilay
  !
  ! Longwave tests - gray radiative equilibrium
  !
  real(wp), parameter :: sigma = 5.670374419e-8_wp, & ! Stefan-Boltzmann constant 
                         D     = 1.66_wp              ! Diffusivity angle, from single-angle RRTMGP solver
  real(wp), dimension(ncol), parameter :: sfc_t     = [(285._wp, icol = 1,ncol/2), & 
                                                       (310._wp, icol = 1,ncol/2)]
  real(wp), dimension(ncol), parameter :: lw_total_tau = [([0.1_wp, 1._wp, 10._wp, 50._wp], icol = 1, ncol/4)]
  real(wp), dimension(ncol), parameter :: sfc_emis     = 1._wp 

  type(ty_optical_props_1scl) :: lw_atmos 
  type(ty_source_func_lw)     :: lw_sources
  type(ty_fluxes_broadband)   :: fluxes
  logical                     :: top_at_1
  real(wp), dimension(ncol,nlay+1), target :: &
                                 ref_flux_up, ref_flux_dn, ref_flux_net

  logical :: passed 

  ! ------------------------------------------------------------------------------------------------------
  
  top_at_1 = .true. 
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, top_at_1, lw_atmos, lw_sources)

  fluxes%flux_up  => ref_flux_up (:,:)
  fluxes%flux_dn  => ref_flux_dn (:,:)
  fluxes%flux_net => ref_flux_net(:,:)
  call stop_on_err(rte_lw(lw_atmos, top_at_1, &
                          lw_sources,         &
                          sfc_emis = reshape(sfc_emis, [1,ncol]), &
                          fluxes = fluxes))

  passed = check_gray_rad_equil(sfc_t, lw_total_tau, top_at_1, &
                                ref_flux_up, ref_flux_net)

  print *, "Unit tests done"
  if(.not. passed) error stop 1

  ! ------------------------------------------------------------------------------------
contains 
  ! ------------------------------------------------------------------------------------
  !
  ! Incoming energy = OLR in gray radiative equilibirum
  !   Equation 6b of Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770 with with f0 = OLR 
  !
  elemental function gray_rad_equil_olr(T, tau)
    real(wp), intent(in) :: T, tau
    real(wp)             :: gray_rad_equil_olr

    gray_rad_equil_olr = (2._wp * sigma * T**4)/(2 + D * tau) 
  end function gray_rad_equil_olr
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
    real(wp), dimension(size(sfc_t)) :: t_lay, olr

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
    olr(:) = gray_rad_equil_olr(sfc_t, lw_total_tau)

    call stop_on_err(sources%alloc(ncol, nlay, atmos))
    sources%sfc_source(:,1) = sigma/pi * sfc_t**4
    if (top_at_1) then
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
    else
      call stop_on_err("Can't do top /= 1 yet")
    end if 
  end subroutine gray_rad_equil
  ! ------------------------------------------------------------------------------------
  !
  ! Check that solutions are in gray radiative equilibrium 
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
                      gray_rad_equil_olr(sfc_t, lw_total_tau))) then 
      call report_err("OLR is not consistent with gray radiative equilibrium")
      check_gray_rad_equil = .false.
    end if 
    !
    ! Check that net fluxes are constant with height -  fairly relaxed threshold w.r.t. spacing() 
    !
    if(.not. allclose(net_flux(:,:), & 
                      spread(net_flux(:,1), dim=2, ncopies=size(net_flux,2)), &
                      tol = 70._wp)) then 
      call report_err("Net flux not constant with tau in gray radiative equilibrium")
      check_gray_rad_equil = .false.
    end if 
  end function check_gray_rad_equil

  ! ------------------------------------------------------------------------------------
end program rte_unit_tests