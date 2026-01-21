! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-  Atmospheric and Environmental Research,
!    Regents of the University of Colorado,
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_rte_util_array
  use mo_rte_kind,      only: wp, wl
  implicit none
  private
  !>
  !> Efficiently set arrays to zero
  !>
  interface set_to_scalar
    module procedure set_to_scalar_1D, set_to_scalar_2D, set_to_scalar_3D, set_to_scalar_4D
  end interface
  interface zero_array
    module procedure zero_array_1D, zero_array_2D, zero_array_3D, zero_array_4D
  end interface
  public :: zero_array, set_to_scalar
contains
 !-------------------------------------------------------------------------------------------------
  ! Initializing arrays to 0
  !-------------------------------------------------------------------------------------------------
  subroutine set_to_scalar_1D(ni, array, value) bind(C, name="set_to_scalar_1D")
    integer,                 intent(in ) :: ni
    real(wp), dimension(ni), intent(out) :: array
    real(wp),                intent(in ) :: value
    ! -----------------------
    integer :: i
    ! -----------------------
    !$acc parallel loop copyout(array)
    !$omp target teams distribute parallel do simd map(from:array)
    do i = 1, ni
      array(i) = value
    end do
  end subroutine set_to_scalar_1D
  ! ----------------------------------------------------------
  subroutine set_to_scalar_2D(ni, nj, array, value) bind(C, name="set_to_scalar_2D")
    integer,                     intent(in ) :: ni, nj
    real(wp), dimension(ni, nj), intent(out) :: array
    real(wp),                    intent(in ) :: value
    ! -----------------------
    integer :: i,j
    ! -----------------------
    !$acc parallel loop collapse(2) copyout(array)
    !$omp target teams distribute parallel do simd collapse(2) map(from:array)
    do j = 1, nj
      do i = 1, ni
        array(i,j) = value
      end do
    end do
  end subroutine set_to_scalar_2D
  ! ----------------------------------------------------------
  subroutine set_to_scalar_3D(ni, nj, nk, array, value) bind(C, name="set_to_scalar_3D")
    integer,                         intent(in ) :: ni, nj, nk
    real(wp), dimension(ni, nj, nk), intent(out) :: array
    real(wp),                        intent(in ) :: value
    ! -----------------------
    integer :: i,j,k
    ! -----------------------
    !$acc parallel loop collapse(3) copyout(array)
    !$omp target teams distribute parallel do simd collapse(3) map(from:array)
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          array(i,j,k) = value
        end do
      end do
    end do
  end subroutine set_to_scalar_3D
  ! ----------------------------------------------------------
  subroutine set_to_scalar_4D(ni, nj, nk, nl, array, value) bind(C, name="set_to_scalar_4D")
    integer,                             intent(in ) :: ni, nj, nk, nl
    real(wp), dimension(ni, nj, nk, nl), intent(out) :: array
    real(wp),                            intent(in ) :: value
    ! -----------------------
    integer :: i,j,k,l
    ! -----------------------
    !$acc parallel loop collapse(4) copyout(array)
    !$omp target teams distribute parallel do simd collapse(4) map(from:array)
    do l = 1, nl
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            array(i,j,k,l) = value
          end do
        end do
      end do
    end do
  end subroutine set_to_scalar_4D
  !-------------------------------------------------------------------------------------------------
  ! Initializing arrays to 0
  !-------------------------------------------------------------------------------------------------
  subroutine zero_array_1D(ni, array) bind(C, name="zero_array_1D")
    integer,                 intent(in ) :: ni
    real(wp), dimension(ni), intent(out) :: array
    ! -----------------------

    call set_to_scalar(ni, array, 0._wp)
  end subroutine zero_array_1D
  ! ----------------------------------------------------------
  subroutine zero_array_2D(ni, nj, array) bind(C, name="zero_array_2D")
    integer,                     intent(in ) :: ni, nj
    real(wp), dimension(ni, nj), intent(out) :: array
    ! -----------------------

    call set_to_scalar(ni, nj, array, 0._wp)
  end subroutine zero_array_2D
  ! ----------------------------------------------------------
  subroutine zero_array_3D(ni, nj, nk, array) bind(C, name="zero_array_3D")
    integer,                         intent(in ) :: ni, nj, nk
    real(wp), dimension(ni, nj, nk), intent(out) :: array
    ! -----------------------

    call set_to_scalar(ni, nj, nk, array, 0._wp)
  end subroutine zero_array_3D
  ! ----------------------------------------------------------
  subroutine zero_array_4D(ni, nj, nk, nl, array) bind(C, name="zero_array_4D")
    integer,                             intent(in ) :: ni, nj, nk, nl
    real(wp), dimension(ni, nj, nk, nl), intent(out) :: array
    ! -----------------------

    call set_to_scalar(ni, nj, nk, nl, array, 0._wp)
  end subroutine zero_array_4D
end module mo_rte_util_array
