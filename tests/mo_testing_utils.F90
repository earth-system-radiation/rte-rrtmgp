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
  use iso_fortran_env,     only: error_unit
  use mo_rte_kind,         only: wp
  use mo_rte_util_array,   only: zero_array
  use mo_optical_props,    only: ty_optical_props_arry, ty_optical_props_1scl, & 
                                 ty_optical_props_2str, ty_optical_props_nstr
  use mo_source_functions, only: ty_source_func_lw
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
  subroutine check_fluxes(flux_1, flux_2, status, message)  
    real(wp), dimension(:,:), intent(in) :: flux_1, flux_2
    logical                              :: status
    character(len=*),         intent(in) :: message

    if(.not. allclose(flux_1, flux_2))  then 
      status = .false. 
      call report_err("    " // trim(message))
    end if 
  end subroutine check_fluxes
  ! ----------------------------------------------------------------------------
  !
  ! Adding transparent (tau = 0) optical properties 
  !   These routines test allocation, validation, incrementing, and 
  !   finalization for optical properties 
  !   Fluxes should not change 
  ! Should these be extended to test end-to-end with GPUs? 
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
    call stop_on_err(atmos%validate())
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
    call stop_on_err(atmos%validate())
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
    call stop_on_err(atmos%validate())
    call transparent%finalize() 
  end subroutine increment_with_nstr 
  ! ----------------------------------------------------------------------------
  !
  ! Vertically reverse
  ! 
  subroutine vr(atmos, sources)
    class(ty_optical_props_arry), intent(inout) :: atmos
    type(ty_source_func_lw), optional, & 
                                  intent(inout) :: sources

    integer :: nlay
    nlay = atmos%get_nlay()

    atmos%tau = atmos%tau(:,nlay:1:-1,:)

    select type (atmos)
      type is (ty_optical_props_2str)
        atmos%ssa = atmos%ssa(:,nlay:1:-1,:)
        atmos%g   = atmos%g  (:,nlay:1:-1,:)
      type is (ty_optical_props_nstr)
        atmos%ssa = atmos%ssa(:,nlay:1:-1,:)
        atmos%p   = atmos%p(:,:,nlay:1:-1,:)
    end select    

    if(present(sources)) then 
      sources%lev_source = sources%lev_source(:,nlay+1:1:-1,:)
      sources%lay_source = sources%lay_source(:,nlay  :1:-1,:)
    end if                           
  end subroutine vr
  ! ----------------------------------------------------------------------------
end module mo_testing_utils
