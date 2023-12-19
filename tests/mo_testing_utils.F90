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
  use mo_optical_props,  only: ty_optical_props_arry, ty_optical_props_1scl, & 
                               ty_optical_props_2str, ty_optical_props_nstr
  implicit none
  public 

  interface allclose
    module procedure allclose_1, allclose_2
  end interface allclose
contains
  ! ----------------------------------------------------------------------------
  !
  ! Compare two arrays; return false if abs(x-y) > tol*spacing(x) for any element
  !
  ! ----------------------------------------------------------------------------
  logical function allclose_1(array1, array2, tol)
    real(wp), dimension(:), intent(in) :: array1, array2
    real(wp), optional,     intent(in) :: tol 
    
    real(wp) :: tolerance 
    if (present(tol)) then 
      tolerance = tol 
    else
      tolerance = 2._wp
    end if 

    allclose_1 = all(abs(array1-array2) <= tolerance * spacing(array1))
  end function allclose_1
  ! ----------------------------------------------------------------------------
  logical function allclose_2(array1, array2, tol)
    real(wp), dimension(:,:), intent(in) :: array1, array2
    real(wp), optional,       intent(in) :: tol 
    
    real(wp) :: tolerance 
    if (present(tol)) then 
      tolerance = tol 
    else
      tolerance = 2._wp
    end if 

    allclose_2= all(abs(array1-array2) <= tolerance * spacing(array1))
  end function allclose_2
  ! ----------------------------------------------------------------------------
  !
  ! Error report - print to screen with or without exit
  !
  ! ----------------------------------------------------------------------------
  subroutine report_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
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
  !
  ! Does this set of surface temperatures and OLRs satisfy gray radiative equillibrium? 
  !
  logical function in_rad_eq(sfc_t, total_tau, OLR)
    ! Assumed rank arrays - Fortran 2018 
    real(wp), dimension(:), intent(in) :: sfc_t, total_tau, OLR

    ! Approximate value of Stefan-Boltzmann constant 
    real(wp), parameter :: sigma = 5.670374419e-8_wp

    in_rad_eq = allclose(sqrt(sqrt((OLR * (1+total_tau))/(2._wp * sigma)) ), & 
                         sfc_t)

  end function in_rad_eq
  ! ----------------------------------------------------------------------------
  !
  ! Adding transparent (tau = 0) optical properties 
  !   These routines test allocation, incrementing, and 
  !   finalization for optical properties 
  !   Fluxes should not change 
  ! Should these be extended to test with GPUs? 
  !
  ! ----------------------------------------------------------------------------
  subroutine increment_with_1scl(atmos)
    class(ty_optical_props_arry), intent(inout) :: atmos

    ! Local variable 
    type(ty_optical_props_1scl) :: transparent
    integer :: ncol, nlay, ngpt 
    ncol = atmos%get_ncol()
    nlay = atmos%get_nlay()
    ngpt = atmos%get_ngpt()

    call stop_on_err(transparent%alloc_1scl(ncol, nlay, atmos))
    call zero_array (ncol, nlay, ngpt, transparent%tau)
    call stop_on_err(transparent%increment(atmos))
    call transparent%finalize() 
  end subroutine increment_with_1scl 
  ! -------
  subroutine increment_with_2str(atmos)
    class(ty_optical_props_arry), intent(inout) :: atmos

    ! Local variable 
    type(ty_optical_props_2str) :: transparent 
    integer :: ncol, nlay, ngpt 
    ncol = atmos%get_ncol()
    nlay = atmos%get_nlay()
    ngpt = atmos%get_ngpt()

    call stop_on_err(transparent%alloc_2str(ncol, nlay, atmos))
    call zero_array (ncol, nlay, ngpt, transparent%tau)
    call zero_array (ncol, nlay, ngpt, transparent%ssa)
    call zero_array (ncol, nlay, ngpt, transparent%g)
    call stop_on_err(transparent%increment(atmos))
    call transparent%finalize() 
  end subroutine increment_with_2str 
  ! -------
  subroutine increment_with_nstr(atmos)
    class(ty_optical_props_arry), intent(inout) :: atmos

    ! Local variable 
    type(ty_optical_props_nstr) :: transparent 
    integer, parameter :: nmom = 4
    integer :: ncol, nlay, ngpt 
    ncol = atmos%get_ncol()
    nlay = atmos%get_nlay()
    ngpt = atmos%get_ngpt()

    call stop_on_err(transparent%alloc_nstr(nmom, ncol, nlay, atmos))
    call zero_array (      ncol, nlay, ngpt, transparent%tau)
    call zero_array (      ncol, nlay, ngpt, transparent%ssa)
    call zero_array (nmom, ncol, nlay, ngpt, transparent%p)
    call stop_on_err(transparent%increment(atmos))
    call transparent%finalize() 
  end subroutine increment_with_nstr 
  ! ----------------------------------------------------------------------------

end module mo_testing_utils
