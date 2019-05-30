! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2019,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_util_array
!
! This module provide utilites for sanitizing input arrays:
!    checking values and sizes
! These are in a module so code can be written for both CPUs and GPUs
! Used only by Fortran classes so routines don't need C bindings and can use assumed-shape
! Currently only for 3D arrays; could extend through overloading to other ranks
!
  use mo_rte_kind,      only: wp
  implicit none
  interface zero_array
    module procedure zero_array_3D, zero_array_4D
  end interface

  private
  public :: zero_array, any_vals_less_than, any_vals_outside
contains
!-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than(array, minVal)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: minVal

    any_vals_less_than = any(array < minVal)
  end function any_vals_less_than
  ! ---------------------------------
  logical function any_vals_outside(array, minVal, maxVal)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: minVal, maxVal

    any_vals_outside = any(array < minVal .or. array > maxVal)
  end function any_vals_outside
! ---------------------------------
! ----------------------------------------------------------
subroutine zero_array_3D(ni, nj, nk, array) bind(C, name="zero_array_3D")
  integer, intent(in) :: ni, nj, nk
  real(wp), dimension(ni, nj, nk), intent(out) :: array
  ! -----------------------
  integer :: i,j,k
  ! -----------------------
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        array(i,j,k) = 0.0_wp
      end do
    end do
  end do

end subroutine zero_array_3D
! ----------------------------------------------------------
subroutine zero_array_4D(ni, nj, nk, nl, array) bind(C, name="zero_array_4D")
  integer, intent(in) :: ni, nj, nk, nl
  real(wp), dimension(ni, nj, nk, nl), intent(out) :: array
  ! -----------------------
  integer :: i,j,k,l
  ! -----------------------
  do l = 1, nl
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          array(i,j,k,l) = 0.0_wp
        end do
      end do
    end do
  end do

end subroutine zero_array_4D
! ----------------------------------------------------------
end module mo_util_array
