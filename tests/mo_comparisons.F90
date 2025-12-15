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
module mo_comparisons
  use mo_rte_kind,         only: wp
  use mo_rte_util_array,   only: zero_array
  use mo_optical_props,    only: ty_optical_props_arry, ty_optical_props_1scl, &
                                 ty_optical_props_2str, ty_optical_props_nstr
  use mo_source_functions, only: ty_source_func_lw
  use mo_testing_utils,    only: stop_on_err, report_err
  implicit none
  private
  public :: allclose, ops_match, check_fluxes
  public :: increment_with_1scl, increment_with_2str, increment_with_nstr, vr

  interface allclose
    module procedure allclose_1, allclose_2, allclose_3, allclose_4
  end interface allclose

  interface ops_match
    module procedure ops_match_1scl, ops_match_2str, ops_match_nstr
  end interface ops_match

  interface check_fluxes
    module procedure check_fluxes_1pair, check_fluxes_2pair
  end interface check_fluxes
contains
  ! ----------------------------------------------------------------------------
  !
  ! Compare two arrays; return false if abs(x-y) > tol*spacing(x) for any element
  !
  ! ----------------------------------------------------------------------------
  logical function allclose_1(tst_array, ref_array, tol)
    real(wp), dimension(:), intent(in) :: tst_array, ref_array
    real(wp), optional,     intent(in) :: tol

    real(wp) :: tolerance
    if (present(tol)) then
      tolerance = tol
    else
      tolerance = 2._wp
    end if

    allclose_1 = all(abs(tst_array-ref_array) <= tolerance * spacing(ref_array))
  end function allclose_1
  ! ----------------------------------------------------------------------------
  logical function allclose_2(tst_array, ref_array, tol)
    real(wp), dimension(:,:), intent(in) :: tst_array, ref_array
    real(wp), optional,       intent(in) :: tol

    real(wp) :: tolerance
    if (present(tol)) then
      tolerance = tol
    else
      tolerance = 2._wp
    end if

    allclose_2= all(abs(tst_array-ref_array) <= tolerance * spacing(ref_array))
  end function allclose_2
  ! ----------------------------------------------------------------------------
  logical function allclose_3(tst_array, ref_array, tol)
    real(wp), dimension(:,:,:), intent(in) :: tst_array, ref_array
    real(wp), optional,         intent(in) :: tol

    real(wp) :: tolerance
    if (present(tol)) then
      tolerance = tol
    else
      tolerance = 2._wp
    end if

    allclose_3= all(abs(tst_array-ref_array) <= tolerance * spacing(ref_array))
  end function allclose_3
  ! ----------------------------------------------------------------------------
  logical function allclose_4(tst_array, ref_array, tol)
    real(wp), dimension(:,:,:,:), intent(in) :: tst_array, ref_array
    real(wp), optional,           intent(in) :: tol

    real(wp) :: tolerance
    if (present(tol)) then
      tolerance = tol
    else
      tolerance = 2._wp
    end if

    allclose_4= all(abs(tst_array-ref_array) <= tolerance * spacing(ref_array))
  end function allclose_4
  ! ----------------------------------------------------------------------------
  !
  ! Compare two sets of optical properties; return false if abs(x-y) > tol*spacing(x) for any element
  !
  ! ----------------------------------------------------------------------------
  logical function ops_match_1scl(tst_values, ref_values, tol)
    class(ty_optical_props_1scl), intent(in) :: tst_values, ref_values
    real(wp), optional,           intent(in) :: tol

    ops_match_1scl = allclose(tst_values%tau, ref_values%tau, tol)
  end function ops_match_1scl
  ! ----------------------------------------------------------------------------
  logical function ops_match_2str(tst_values, ref_values, tol)
    class(ty_optical_props_2str), intent(in) :: tst_values, ref_values
    real(wp), optional,           intent(in) :: tol

    ops_match_2str = allclose(tst_values%tau, ref_values%tau, tol) .and. &
                     allclose(tst_values%ssa, ref_values%ssa, tol) .and. &
                     allclose(tst_values%g  , ref_values%g  , tol)
  end function ops_match_2str
  ! ----------------------------------------------------------------------------
  logical function ops_match_nstr(tst_values, ref_values, tol)
    class(ty_optical_props_nstr), intent(in) :: tst_values, ref_values
    real(wp), optional,           intent(in) :: tol

    ops_match_nstr = allclose(tst_values%tau, ref_values%tau, tol) .and. &
                     allclose(tst_values%ssa, ref_values%ssa, tol) .and. &
                     allclose(tst_values%p  , ref_values%p  , tol)
  end function ops_match_nstr
  ! ----------------------------------------------------------------------------
  subroutine check_fluxes_1pair(flux_1, flux_2, status, message)
    real(wp), dimension(:,:), intent(in) :: flux_1, flux_2
    logical                              :: status
    character(len=*),         intent(in) :: message

    if(.not. allclose(flux_1, flux_2))  then
      status = .false.
      print *, "check_fluxes: max diffs rel. to scaling: ", &
        maxval(abs(flux_1 - flux_2)/spacing(flux_1))
      call report_err("    " // trim(message))
    end if
  end subroutine check_fluxes_1pair
  ! ----------------------------------------------------------------------------
  subroutine check_fluxes_2pair(flux_1, flux_2, flux_3, flux_4, status, message)
    real(wp), dimension(:,:), intent(in) :: flux_1, flux_2, flux_3, flux_4
    logical                              :: status
    character(len=*),         intent(in) :: message

    if(.not. (allclose(flux_1, flux_2) .and. &
              allclose(flux_3, flux_4))) then
      status = .false.
      print *, "check_fluxes: max diffs rel. to scaling: ", &
        maxval(abs(flux_1 - flux_2)/spacing(flux_1)), &
        maxval(abs(flux_3 - flux_4)/spacing(flux_3))
      call report_err("    " // trim(message))
    end if
  end subroutine check_fluxes_2pair
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
  ! Vertically reverse optical properties
  !
  subroutine vr(atmos, sources)
    class(ty_optical_props_arry), intent(inout) :: atmos
    type(ty_source_func_lw), optional, &
                                  intent(inout) :: sources

    integer :: nlay
    ! -----------------------
    nlay = atmos%get_nlay()

    call atmos%set_top_at_1(.not. atmos%top_is_at_1())
    atmos%tau(:,:,:) = atmos%tau(:,nlay:1:-1,:)

    select type (atmos)
      type is (ty_optical_props_2str)
        atmos%ssa(:,:,:) = atmos%ssa(:,nlay:1:-1,:)
        atmos%g  (:,:,:)= atmos%g   (:,nlay:1:-1,:)
      type is (ty_optical_props_nstr)
        atmos%ssa(:,:,:) = atmos%ssa(:,nlay:1:-1,:)
        atmos%p(:,:,:,:) = atmos%p(:,:,nlay:1:-1,:)
    end select

    if(present(sources)) then
      sources%lev_source(:,:,:) = sources%lev_source(:,nlay+1:1:-1,:)
      sources%lay_source(:,:,:) = sources%lay_source(:,nlay  :1:-1,:)
    end if
  end subroutine vr
  ! ----------------------------------------------------------------------------
end module mo_comparisons
