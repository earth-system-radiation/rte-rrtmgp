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
  use mo_rte_util_array,only : zero_array, set_to_scalar
  use mo_rte_gas_optics_utils, &
                        only: grav, compute_Planck_source, compute_layer_mass

  implicit none

  private
  public :: compute_tau

  !
  ! Physical constants
  !
  real(wp), parameter :: planck_h     = 6.626075540e-34_wp    ! Planck's constant
  real(wp), parameter :: lightspeed   = 2.99792458e8_wp       ! Speed of light
  real(wp), parameter :: boltzmann_k  = 1.38065812e-23_wp     ! Boltzman thermodynamic constant

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
    real(wp), dimension(size(play,1), size(play,2)) :: p_scaling

    !$acc        data create(   p_scaling)
    !$omp target data map(alloc:p_scaling)

    ! Apply pressure broadening if pref input is non-zero
    if (pref /= 0._wp) then
      !$acc                         parallel loop    collapse(2)
      !$omp target teams distribute parallel do simd collapse(2)
      do ilay = 1, nlay
        do icol = 1, ncol
          p_scaling(icol, ilay) =  play(icol, ilay) / pref
        end do
      end do
    else
      call set_to_scalar(ncol, nlay, p_scaling, 1._wp)
    end if

    !$acc                         parallel loop    collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          tau(icol, ilay, inu) = &
            p_scaling(icol, ilay) * sum( [(layer_mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas)] )
        end do
      end do
    end do

    !$acc end data
    !$omp end target data

  end subroutine compute_tau

end module mo_optics_ssm_kernels
