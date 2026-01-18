module mo_testing_utils
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: stop_on_err, report_err
  ! ----------------------------------------------------------------------------
  !
  ! Error report - print to screen with or without exit
  !
  ! ----------------------------------------------------------------------------
contains
  subroutine report_err(error_msg)
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
    end if
  end subroutine report_err
  ! ----------------------------------------------------------------------------
  subroutine stop_on_err(error_msg)
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "unit tests stopping"
      error stop 1
    end if
  end subroutine stop_on_err
end module mo_testing_utils
