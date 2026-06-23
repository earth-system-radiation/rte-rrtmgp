! This is an implementation of a simple spectral model
!
! Contacts: Andrew Williams and Robert Pincus
! email:  andrewwilliams@ucsd.edu
!!
! Copyright 2025-,  ... and Trustees of Columbia University.  All right reserved.
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
  implicit none
  private
  public :: compute_tau
contains
  !
  ! Computes optical depth from atmospheric state and gas optics
  !
  subroutine compute_tau(ncol, nlay, nnu, ngas,   &
                         absorption_coeffs, play, pref, layer_mass, &
                         tau) bind(C, name="ssm_compute_tau_absorption")
    integer,  &
      intent(in ) :: ncol, nlay, nnu, ngas
    real(wp), dimension(ngas, nnu), &
      intent(in ) :: absorption_coeffs
    real(wp), dimension(ncol, nlay), &
      intent(in   ) :: play     !! layer pressures [Pa]; (ncol,nlay)
    real(wp), dimension(ngas, ncol, nlay), &
      intent(in ) :: layer_mass !! !! mass of atm layer [kg/m2]; (ngas, ncol, nlay)
    real(wp), &
      intent(in ) :: pref       !! reference pressure [Pa], do foreign broadening if this is non-zero
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: tau

    ! Local variables
    integer :: icol, ilay, inu, igas

    ! Apply pressure broadening if pref input is non-zero
    if (pref /= 0._wp) then
      !$acc                         parallel loop    collapse(3)
      !$omp target teams distribute parallel do simd collapse(3)
      do inu = 1, nnu
        do ilay = 1, nlay
          do icol = 1, ncol
            tau(icol, ilay, inu) = &
              sum( [(layer_mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas)] ) &
              * play(icol, ilay) / pref
          end do
        end do
      end do
    else
      !$acc                         parallel loop    collapse(3)
      !$omp target teams distribute parallel do simd collapse(3)
      do inu = 1, nnu
        do ilay = 1, nlay
          do icol = 1, ncol
            tau(icol, ilay, inu) = &
              sum( [(layer_mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas)] )
          end do
        end do
      end do
    end if
  end subroutine compute_tau

end module mo_optics_ssm_kernels
