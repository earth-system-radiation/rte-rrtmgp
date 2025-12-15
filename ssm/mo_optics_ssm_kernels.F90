! See documentation in other modules...
!
! Contacts: could provide here... names, emails or web site
!
! Copyright 2025-,  ... Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Kernels
!>
!> Details here
! -------------------------------------------------------------------------------------------------
module mo_optics_ssm_kernels
  use mo_rte_kind,      only : wp, wl
  use mo_rte_util_array,only : zero_array
  implicit none
  interface compute_Planck_source
    module procedure compute_Planck_source_1D, compute_Planck_source_2D
  end interface

  private
  public :: compute_tau, compute_Planck_source

  !
  ! Physical constants
  !
  real(wp), parameter :: planck_h     = 6.626075540e-34_wp    ! Planck's constant
  real(wp), parameter :: lightspeed   = 2.99792458e8_wp       ! Speed of light
  real(wp), parameter :: boltzmann_k  = 1.38065812e-23_wp     ! Boltzman thermodynamic constant

contains
  !
  ! Doesn't account for pressure broadening yet...
  !
  subroutine compute_tau(ncol, nlay, nnu, ngas,   &
                         absorption_coeffs, layer_mass, p_scaling, &
                         tau) bind(C, name="ssm_compute_tau_absorption")
    integer,  &
      intent(in ) :: ncol, nlay, nnu, ngas
    real(wp), dimension(ngas, nnu), &
      intent(in ) :: absorption_coeffs
    real(wp), dimension(ngas, ncol, nlay), &
      intent(in ) :: layer_mass
    real(wp), dimension(ncol, nlay), &
      intent(in ) :: p_scaling
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: tau

    ! Local variables
    integer :: icol, ilay, inu, igas


    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          tau(icol, ilay, inu) = &
            p_scaling(icol, ilay) * sum( [(layer_mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas)] )
        end do
      end do
    end do

  end subroutine compute_tau
  ! -------------------------------------------------------------------------------------------------
  !
  ! Maybe make an elemental function then call with 1D and 2D wrappers?
  !
  elemental function B_nu(T, nu)
    real(wp), intent(in) :: T, nu
    real(wp)             :: B_nu
    B_nu = 100._wp*2._wp*planck_h*((nu*100._wp)**3)*(lightspeed**2) / &
                  exp( (planck_h * lightspeed * nu * 100._wp) / (boltzmann_k * T) - 1._wp)
  end function
  ! -------
  subroutine compute_Planck_source_2D(&
      ncol, nlay, nnu, &
      nus, dnus, T, &
      source) bind(C, name="ssm_compute_Planck_source_2D")
    integer,  &
      intent(in ) :: ncol, nlay, nnu
    real(wp), dimension(nnu), &
      intent(in ) :: nus, dnus
    real(wp), dimension(ncol, nlay), &
      intent(in ) :: T
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: source

     ! Local variables
     integer :: icol, ilay, inu

    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          source(icol, ilay, inu) = B_nu(T(icol, ilay), nus(inu)) * dnus(inu)
        end do
      end do
    end do

  end subroutine compute_Planck_source_2D
  ! -------
  subroutine compute_Planck_source_1D(&
      ncol, nnu, &
      nus, dnus, T, &
      source) bind(C, name="ssm_compute_Planck_source_1D")
    integer,  &
      intent(in ) :: ncol, nnu
    real(wp), dimension(nnu), &
      intent(in ) :: nus, dnus
    real(wp), dimension(ncol), &
      intent(in ) :: T
    real(wp), dimension(ncol, nnu), &
      intent(out) :: source

     ! Local variables
     integer :: icol, ilay, inu

    do inu = 1, nnu
      do icol = 1, ncol
        source(icol, inu) = B_nu(T(icol), nus(inu)) * dnus(inu)
      end do
    end do

  end subroutine compute_Planck_source_1D
  ! -------

end module mo_optics_ssm_kernels
