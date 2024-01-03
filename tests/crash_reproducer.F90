program crash_reproducer
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_long, c_int, c_bool
  implicit none
  integer, parameter :: dp = c_double, sp = c_float, i8 = c_long, i4 = c_int
  integer, parameter :: wp = dp

  real(wp), parameter :: pi = acos(-1._wp)
  integer,  parameter :: ncol = 8, nlay = 16
  integer             :: icol, ilay
  logical             :: top_at_1 = .true. 
  !
  ! Longwave tests - gray radiative equilibrium
  !

  real(wp), parameter :: sigma = 5.670374419e-8_wp, & ! Stefan-Boltzmann constant 
                         D     = 1.66_wp              ! Diffusivity angle, from single-angle RRTMGP solver
  real(wp), dimension(  ncol), parameter :: sfc_t     = [(285._wp, icol = 1,ncol/2), & 
                                                         (310._wp, icol = 1,ncol/2)]
  real(wp), dimension(  ncol), parameter :: lw_total_tau = [0.1_wp, 1._wp, 10._wp, 50._wp, &
                                                            0.1_wp, 1._wp, 10._wp, 50._wp] ! Would be nice to parameterize 
  
  call gray_rad_equil(sfc_t, lw_total_tau, nlay, top_at_1)
  print *, gray_rad_equil_olr(sfc_t, lw_total_tau)
contains 
  ! ------------------------------------------------------------------------------------
  !
  ! Define an atmosphere in gray radiative equillibrium 
  !   See, for example, section 2 of Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770
  !
  subroutine gray_rad_equil(sfc_t, total_tau, nlay, top_at_1)
    real(wp), dimension(:), intent(in) :: sfc_t, total_tau
    integer,                intent(in) :: nlay 
    logical,                intent(in) :: top_at_1

    integer                          :: ncol
    real(wp), dimension(size(sfc_t)) :: t_lay, olr

    ncol = size(sfc_t)

    !
    ! Longwave sources - for broadband these are sigma/pi T^4
    !   (isotropic radiation)
    !
    olr = gray_rad_equil_olr(sfc_t, lw_total_tau)

  end subroutine gray_rad_equil
  ! ------------------------------------------------------------------------------------
  !
  ! Incoming energy = OLR in gray radiative equilibirum
  !   Equation 6b of Weaver and Rmanathan 1995 https://doi.org/10.1029/95JD00770 with with f0 = OLR 
  !
  function gray_rad_equil_olr(T, tau)
    real(wp), dimension(:), intent(in) :: T, tau
    real(wp), dimension(size(T))       :: gray_rad_equil_olr

    gray_rad_equil_olr(:) = (2._wp * sigma * T(:)**4)/(2 + D * tau(:)) 
  end function gray_rad_equil_olr

end program crash_reproducer 