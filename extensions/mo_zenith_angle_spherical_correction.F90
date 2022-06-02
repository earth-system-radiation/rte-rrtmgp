! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2020-,  Atmospheric and Environmental Research,
! Regents of the University of Colorado, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Routines for solar zenith angle accounting for spherical geometry
!
! -------------------------------------------------------------------------------------------------

module mo_zenith_angle_spherical_correction
  use mo_rte_kind,       only: wp, wl
  use mo_rte_config,     only: check_extents, check_values
  use mo_rte_util_array, only: extents_are, any_vals_outside, any_vals_less_than
  implicit none
  private
  public :: zenith_angle_with_height

  real(wp), protected :: planet_radius = 6371.23_wp * 1.e3_wp ! Planetary radius [m]
contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! Compute cosine of solar zenith angle (min 0) as a function of height given value at a specified height.
  !
  function zenith_angle_with_height(ref_alt, ref_mu, alt, mu) result(error_msg)
    real(wp), dimension(:),   intent(in   ) :: ref_alt ! reference altitude at which solar zenith angle is provided [m]
    real(wp), dimension(:),   intent(in   ) :: ref_mu  ! cosine of the solar zenith angle at reference altitude
    real(wp), dimension(:,:), intent(in   ) :: alt    ! level altitude grid [m]
    real(wp), dimension(:,:), intent(  out) :: mu
    character(len=128)                      :: error_msg

    integer  :: ncol, nlay, icol, ilay
    real(wp) :: sin_theta2
    ! ------------------------------------
    error_msg = ""
    ncol = size(alt,1)
    nlay = size(alt,2)

    if(check_extents) then
      if(.not. extents_are(ref_alt, ncol)) &
        error_msg = "zenith_angle_with_height: ref_alt, alt have different number of columns"
      if(.not. extents_are(ref_mu, ncol)) &
        error_msg = "zenith_angle_with_height: ref_mu, alt have different number of columns"
      if(.not. extents_are(mu,    ncol, nlay)) &
        error_msg = "zenith_angle_with_height: mu, alt have different number of columns"
    end if
    if(len_trim(error_msg) /= 0) return

    if(check_extents) then
      if(any_vals_less_than(ref_alt, -planet_radius)) &
        error_msg = "zenith_angle_with_height: values of ref_alt must be larger than -the planetary radius"
      if(any_vals_outside(ref_mu, -1._wp, 1._wp)) &
        error_msg = "zenith_angle_with_height: values of ref_mu must be in [-1, 1]"
      if(any_vals_less_than(alt, -planet_radius)) &
        error_msg = "zenith_angle_with_height: values of alt must be larger than -the planetary radius"
    end if
    if(len_trim(error_msg) /= 0) return
    ! ------------------------------------
    !$acc                         parallel loop    collapse(2) &
    !$acc copyin(ref_alt, ref_mu, alt) copyout(mu)
    !$omp target teams distribute parallel do simd collapse(2) &
    !$omp map(to:ref_alt, ref_mu, alt) map(from:mu)
    do ilay=1, nlay
      do icol = 1, ncol
        sin_theta2 = (1-ref_mu(icol)**2) * &
                     ((planet_radius + ref_alt(icol)) / &
                      (planet_radius + alt(icol,ilay)))**2
        if(sin_theta2 < 1._wp) then
          mu(icol, ilay) = sqrt(1._wp - sin_theta2)
        else
          mu(icol, ilay) = 0._wp
        end if
      end do
    end do
  end function zenith_angle_with_height
  ! -------------------------------------------------------------------------------------------------
  !
  ! Set the planetary radius used for computing solar zenith angle
  !
  function set_planet_radius(radius) result (error_msg)
    real(wp), intent(in) :: radius
    character(len=128)   :: error_msg

    error_msg = ""
    if(radius .le. 0._wp) then
      error_msg = "set_planet_radius: radius must be > 0"
      return
    end if

    planet_radius = radius
  end function set_planet_radius
  ! -------------------------------------------------------------------------------------------------
end module mo_zenith_angle_spherical_correction
