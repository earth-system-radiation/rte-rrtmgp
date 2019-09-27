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
module mo_rte_util_array
!
! This module provide utilites for sanitizing input arrays:
!    checking values and sizes
! These are in a module so code can be written for both CPUs and GPUs
! Used only by Fortran classes so routines don't need C bindings and can use assumed-shape
!
  use mo_rte_kind,      only: wp
  implicit none
  interface any_vals_less_than
    module procedure any_vals_less_than_1D, any_vals_less_than_2D, any_vals_less_than_3D
  end interface
  interface any_vals_outside
    module procedure any_vals_outside_1D, any_vals_outside_2D, any_vals_outside_3D
  end interface
  interface zero_array
    module procedure zero_array_1D, zero_array_3D, zero_array_4D
  end interface

  private
  public :: zero_array, any_vals_less_than, any_vals_outside
contains
  !-------------------------------------------------------------------------------------------------
  ! Values less than a floor
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_1D(array, check_value)
    real(wp), dimension(:), intent(in) :: array
    real(wp),               intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    !$acc end kernels

    any_vals_less_than_1D = (minValue < check_value)

  end function any_vals_less_than_1D
!-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_2D(array, check_value)
    real(wp), dimension(:,:), intent(in) :: array
    real(wp),                 intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    !$acc end kernels

    any_vals_less_than_2D = (minValue < check_value)

  end function any_vals_less_than_2D
!-------------------------------------------------------------------------------------------------
  logical function any_vals_less_than_3D(array, check_value)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: check_value

    real(wp) :: minValue

    !$acc kernels copyin(array)
    minValue = minval(array)
    !$acc end kernels

    any_vals_less_than_3D = (minValue < check_value)

  end function any_vals_less_than_3D
  !-------------------------------------------------------------------------------------------------
  ! Values outside a range
  !-------------------------------------------------------------------------------------------------
  logical function any_vals_outside_1D(array, minVal, maxVal)
    real(wp), dimension(:), intent(in) :: array
    real(wp),               intent(in) :: minVal, maxVal

#ifdef _OPENACC
      ! Compact version using intrinsics below
      ! but an explicit loop is the only current solution on GPUs
    real(wp) :: minValue, maxValue
    integer  :: i
    minValue = minVal
    maxValue = maxVal
    !$acc parallel loop copyin(array) reduction(min:minValue) reduction(max:maxValue)
    do i = 1, size(array,1)
      minValue = min(array(i), minValue)
      maxValue = max(array(i), maxValue)
    end do
    any_vals_outside_1D = (minValue < minVal .or. maxValue > maxVal)
#else
    any_vals_outside_1D = any(array < minVal .or. array > maxVal)
#endif
  end function any_vals_outside_1D
! ----------------------------------------------------------
  logical function any_vals_outside_2D(array, minVal, maxVal)
    real(wp), dimension(:,:), intent(in) :: array
    real(wp),                 intent(in) :: minVal, maxVal

#ifdef _OPENACC
      ! Compact version using intrinsics below
      ! but an explicit loop is the only current solution on GPUs
    real(wp) :: minValue, maxValue
    integer  :: i, j
    minValue = minVal
    maxValue = maxVal
    !$acc parallel loop collapse(2) copyin(array) reduction(min:minValue) reduction(max:maxValue)
    do j = 1, size(array,2)
      do i = 1, size(array,1)
        minValue = min(array(i,j), minValue)
        maxValue = max(array(i,j), maxValue)
      end do
    end do
    any_vals_outside_2D = (minValue < minVal .or. maxValue > maxVal)
#else
    any_vals_outside_2D = any(array < minVal .or. array > maxVal)
#endif
  end function any_vals_outside_2D
! ----------------------------------------------------------
  logical function any_vals_outside_3D(array, minVal, maxVal)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp),                   intent(in) :: minVal, maxVal

#ifdef _OPENACC
      ! Compact version using intrinsics below
      ! but an explicit loop is the only current solution on GPUs
    real(wp) :: minValue, maxValue
    integer  :: i, j, k
    minValue = minVal
    maxValue = maxVal
    !$acc parallel loop collapse(3) copyin(array) reduction(min:minValue) reduction(max:maxValue)
    do k = 1, size(array,3)
      do j = 1, size(array,2)
        do i = 1, size(array,1)
          minValue = min(array(i,j,k), minValue)
          maxValue = max(array(i,j,k), maxValue)
        end do
      end do
    end do
    any_vals_outside_3D = (minValue < minVal .or. maxValue > maxVal)
#else
    any_vals_outside_3D = any(array < minVal .or. array > maxVal)
#endif
  end function any_vals_outside_3D
  !-------------------------------------------------------------------------------------------------
  ! Initializing arrays to 0
  !-------------------------------------------------------------------------------------------------
  subroutine zero_array_1D(ni, array) bind(C, name="zero_array_1D")
    integer,                 intent(in ) :: ni
    real(wp), dimension(ni), intent(out) :: array
    ! -----------------------
    integer :: i
    ! -----------------------
    !$acc parallel loop copyout(array)
    do i = 1, ni
      array(i) = 0.0_wp
    end do
  end subroutine zero_array_1D
  ! ----------------------------------------------------------
  subroutine zero_array_3D(ni, nj, nk, array) bind(C, name="zero_array_3D")
    integer,                         intent(in ) :: ni, nj, nk
    real(wp), dimension(ni, nj, nk), intent(out) :: array
    ! -----------------------
    !$acc kernels
    array = 0.0_wp
    !$acc end kernels

  end subroutine zero_array_3D
  ! ----------------------------------------------------------
  subroutine zero_array_4D(ni, nj, nk, nl, array) bind(C, name="zero_array_4D")
    integer,                             intent(in ) :: ni, nj, nk, nl
    real(wp), dimension(ni, nj, nk, nl), intent(out) :: array
    ! -----------------------
    !$acc kernels
    array = 0.0_wp
    !$acc end kernels

  end subroutine zero_array_4D
! ----------------------------------------------------------
end module mo_rte_util_array
