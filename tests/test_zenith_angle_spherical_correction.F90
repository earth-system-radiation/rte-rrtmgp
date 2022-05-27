subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "test_solar_zenith_angle stopping"
    error stop 1
  end if
end subroutine stop_on_err
! ------------------------------------------------------------------------------
program test_solar_zenith_angle
  use mo_rte_kind,                          only: wp, wl
  use mo_zenith_angle_spherical_correction, only: zenith_angle_with_height

  integer, parameter :: nlay = 80
  real(wp), dimension(1        ) :: refAlt, refMu
  real(wp), dimension(1, 0:nlay) ::    alt,    mu
  integer            :: i

  refAlt(:) = nlay * 1000._wp   ! Reference altitude is TOA
  refMu (:) = .02
  alt(1,:) =  [(i, i = 0, nlay)] * 1000._wp
  call stop_on_err(zenith_angle_with_height(refAlt, refMu, alt, mu))
  print *, "alt     mu0"
  do i = nlay, 0, -1
    print *, int(alt(1,i)/1000._wp), mu(1,i)
  end do
end program test_solar_zenith_angle
