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
  private
  public :: compute_tau, compute_Planck_source
contains
  !
  ! Doesn't account for pressure broadening yet...
  !
  subroutine compute_tau(ncol, nlay, nnu, ngas,   &
                         absorption_coeffs, mass, &
                         tau) bind(C, name="ssm_compute_tau_absorption")
    integer,  &
      intent(in ) :: ncol, nlay, nnu, ngas
    real(wp), dimension(ngases, nnu), &
      intent(in ) :: absorption_coeffs
    real(wp), dimension(ngas, ncol, nlay), &
      intent(in ) :: mass
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: tau

    ! Local variables
    integer :: icol, ilay, inu, igas


    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          tau(icol, ilay) = &
            sum([(mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas))])
        end do
      end do
    end do

  end subroutine compute_tau
  ! -------------------------------------------------------------------------------------------------
  !
  ! Maybe make an elemental function then call with 1D and 2D wrappers?
  !
  subroutine compute_Planck_source(&
      ncol, nlay, nnu, &
      nus, dnus, tlay, &
      source) bind(C, name="ssm_compute_Planck_source")

  end subroutine compute_Planck_source

end module mo_optics_ssm_kernels
