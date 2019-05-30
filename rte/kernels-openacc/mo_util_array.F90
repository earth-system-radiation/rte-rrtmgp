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

    real(wp) :: minValue
    integer  :: i, j, k

    ! This could be written far more compactly as
    !       any_vals_less_than = any(array < minVal)
    ! but an explicit loop also works on GPUs
    minValue = minVal
    !$acc parallel loop collapse(3) copyin(array) reduction(min:minValue)
    do k = 1, size(array,3)
      do j = 1, size(array,2)
        do i = 1, size(array,1)
          minValue = min(array(i,j,k), minValue)
        end do
      end do
    end do
    any_vals_less_than = (minValue < minVal)
  end function any_vals_less_than
  ! ---------------------------------
  logical function any_vals_outside(array, minVal, maxVal)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: minVal, maxVal

    real(wp) :: minValue, maxValue
    integer  :: i, j, k

    ! This could be written far more compactly as
    !   any_vals_outside = any(array < minVal .or. array > maxVal)
    ! but an explicit loop also works on GPUs
    minValue = minVal
    maxValue = maxVal
    !$acc parallel loop collapse(3) copyin(array) reduction(min:minValue) reduction(max:maxValue)
    do k = 1, size(array,3)
      do j = 1, size(array,2)
        do i = 1, size(array,1)
          minValue = min(array(i,j,k), minValue)
          maxValue = max(array(i,j,k), minValue)
        end do
      end do
    end do
    any_vals_outside = (minValue < minVal .or. maxValue > maxVal)
  end function any_vals_outside
! ---------------------------------
  ! ----------------------------------------------------------
  subroutine zero_array_3D(ni, nj, nk, array) bind(C, name="zero_array_3D")
    integer, intent(in) :: ni, nj, nk
    real(wp), dimension(ni, nj, nk), intent(out) :: array
    ! -----------------------
    integer :: i,j,k
    ! -----------------------
    !$acc parallel loop collapse(3) &
    !$acc&     copyout(array(:ni,:nj,:nk))
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
    !$acc parallel loop collapse(4) &
    !$acc&     copyout(array(:ni,:nj,:nk,:nl))
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
  !----------------------------------------------------------
end module mo_util_array
