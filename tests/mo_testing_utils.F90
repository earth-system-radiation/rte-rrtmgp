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
module mo_testing_utils
  use iso_fortran_env,   only: error_unit
  use mo_rte_kind,       only: wp
  use mo_rte_util_array, only: zero_array
  implicit none
contains
  ! ----------------------------------------------------------------------------
  !
  ! Compare two arrays; return false if abs(x-y) > tol*spacing(x) for any element
  !
  logical function allclose(array1, array2, tol)
    real(wp), dimension(:,:), intent(in) :: array1, array2
    real(wp), optional,       intent(in) :: tol 
    
    real(wp) :: tolerance 
    if (present(tol)) then 
      tolerance = tol 
    else
      tolerance = 2._wp
    end if 

    allclose = all(abs(array1-array2) <= tolerance * spacing(array1))
  end function allclose 
  ! ----------------------------------------------------------------------------
  subroutine report_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      failed = .true.
    end if
  end subroutine report_err
  ! ----------------------------------------------------------------------------
  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "unit tests stopping"
      error stop 1
    end if
  end subroutine stop_on_err
  ! ----------------------------------------------------------------------------
end module mo_testing_utils
