program crash_reproducer
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_long, c_int, c_bool
  implicit none
  integer, parameter :: wp = c_double

  integer,  parameter :: ncol = 8, nlay=16
  integer             :: icol
  real(wp), parameter :: sigma = 5.670374419e-8_wp, & ! Stefan-Boltzmann constant 
                         D     = 1.66_wp              ! Diffusivity angle, from single-angle RRTMGP solver
  real(wp), dimension(  ncol), parameter :: sfc_t     = [(285._wp, icol = 1,ncol/2), & 
                                                         (310._wp, icol = 1,ncol/2)]
  real(wp), dimension(  ncol), parameter :: lw_total_tau = [0.1_wp, 1._wp, 10._wp, 50._wp, &
                                                            0.1_wp, 1._wp, 10._wp, 50._wp] ! Would be nice to parameterize 
  
  !
  ! Problem demonstration - the vector-valued function gray_rad_equil_olr works when called directly...
  !
  print '("Return vector-valued function: ", 8(f6.2, 2x))', gray_rad_equil_olr(sfc_t, lw_total_tau)
  !
  ! ... but not when called from within a subroutine, returning instead a bounds error 
  !
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, .true.)
contains 
  ! ------------------------------------------------------------------------------------
  !
  ! Function version 
  !
  function gray_rad_equil_olr(T, tau)
    real(wp), dimension(:), intent(in) :: T, tau
    real(wp), dimension(size(T))       :: gray_rad_equil_olr

    gray_rad_equil_olr = (2._wp * sigma * T**4)/(2 + D*tau) 
  end function gray_rad_equil_olr
  ! ------------------------------------------------------------------------------------
  !
  ! Subroutine version 
  !
  subroutine gray_rad_equil(sfc_t, total_tau, nlay, top_at_1)
    real(wp), dimension(:), intent(in) :: sfc_t, total_tau
    integer,                intent(in) :: nlay 
    logical,                intent(in) :: top_at_1
 
    integer :: i 

    print '("Looped from within subroutine: ", 8(f6.2, 2x))', & 
       [(gray_rad_equil_olr(sfc_t(i), lw_total_tau(i)), i = 1, size(sfc_t))]
    print '("Called from within subroutine: ", 8(f6.2, 2x))', gray_rad_equil_olr(sfc_t, lw_total_tau)
  end subroutine gray_rad_equil
  ! ------------------------------------------------------------------------------------
  !
  ! Function version - scalar
  !
  function gray_rad_equil_olr(T, tau)
    real(wp),  intent(in) :: T, tau
    real(wp),             :: gray_rad_equil_olr

    gray_rad_equil_olr = (2._wp * sigma * T**4)/(2 + D*tau) 
  end function gray_rad_equil_olr
  ! ------------------------------------------------------------------------------------
end program crash_reproducer 