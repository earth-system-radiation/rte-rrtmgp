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
program optical_prop_unit_tests
  !
  ! Unit tests for RTE optical properties
  !   Incrementing with tranparent medium (tau=0) doesn't change optical props
  !
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props_arry, & 
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_util_array,     only: zero_array
  use mo_testing_utils,      only: allclose, ops_match, stop_on_err, report_err, &
                                   increment_with_1scl, increment_with_2str, increment_with_nstr

  type(ty_optical_props_1scl) :: ref_1scl, tst_1scl
  type(ty_optical_props_2str) :: ref_2str, tst_2str
  type(ty_optical_props_nstr) :: ref_nstr, tst_nstr
  integer,  parameter :: ncol = 4, nlay = 8, nmom = 4
  integer             :: icol, ilay, imom
  logical             :: passed
  real(wp), dimension(  ncol), parameter :: total_tau = [0.1_wp, 1._wp, 10._wp, 50._wp] 
  real(wp),                    parameter :: g = 0.85_wp, ssa = 1._wp - 1.e-4_wp
  ! ----------------------------------------------------------------------------
  print *, "Optical properties unit testing"
  !
  ! Set up a gray spectral distribution - one band, one g-point
  !
  call stop_on_err(ref_1scl%init(band_lims_wvn = reshape([0._wp, 3250._wp], shape = [2, 1]), & 
                                 band_lims_gpt = reshape([1,     1],        shape = [2, 1]), & 
                                 name = "Gray atmosphere"))
  call stop_on_err(ref_1scl%alloc_1scl(ncol, nlay))
  print '("  Problem size: (ncol, nlay, nband, ngpt): ", 4(i2, 2x))', & 
    ref_1scl%get_ncol(), ref_1scl%get_nlay(), ref_1scl%get_nband(), ref_1scl%get_ngpt()
  !
  ! Divide optical depth evenly among layers 
  !
  ref_1scl%tau(1:ncol,1:nlay,1) = spread(total_tau(1:ncol)/real(nlay, wp), dim=2, ncopies=nlay)
  !
  ! 2- and n-stream optical properties 
  !
  call stop_on_err(ref_2str%alloc_2str(ncol, nlay, ref_1scl))
  ref_2str%tau = ref_1scl%tau
  ref_2str%ssa = ssa
  ref_2str%g   = g

  call stop_on_err(ref_nstr%alloc_nstr(nmom, ncol, nlay, ref_1scl))
  ref_nstr%tau = ref_1scl%tau
  ref_nstr%ssa = ssa
  ! Henyey-Greenstein phase function 
  do imom = 1, nmom
    ref_nstr%p(imom,:,:,:) = g**imom
  end do 

  passed = .true. 
  ! ----------------------------------------------------------------------------
  !
  ! Incrementing with transparent (tau=0) sets of optical properties 
  !
  ! ----------------------------------------------------------------------------
  print *, "  Incrementing 1scl"
  !
  ! Increment 1scl 
  !
  call make_copy_1scl
  call increment_with_1scl(tst_1scl)
  if(.not. ops_match(tst_1scl, ref_1scl)) then 
    call report_err("1scl+1scl fails")
    passed = .false. 
  end if 

  call make_copy_1scl
  call increment_with_2str(tst_1scl)
  if(.not. ops_match(tst_1scl, ref_1scl)) then 
    call report_err("1scl+2str fails")
    passed = .false. 
  end if 

  call make_copy_1scl
  call increment_with_nstr(tst_1scl)
  if(.not. ops_match(tst_1scl, ref_1scl)) then 
    call report_err("1scl+nstr fails")
    passed = .false. 
  end if 

  call tst_1scl%finalize()
  ! ----------------------------------------------------------------------------
  print *, "  Incrementing 2str"
  !
  ! Increment 2str 
  !
  call make_copy_2str
  call increment_with_1scl(tst_2str)
  if(.not. ops_match(tst_2str, ref_2str)) then 
    call report_err("2str+1scl fails")
    passed = .false. 
  end if 

  call make_copy_2str
  call increment_with_2str(tst_2str)
  if(.not. ops_match(tst_2str, ref_2str)) then 
    call report_err("2str+2str fails")
    passed = .false. 
  end if 

  call make_copy_2str
  call increment_with_nstr(tst_2str)
  if(.not. ops_match(tst_2str, ref_2str)) then 
    call report_err("2str+nstr fails")
    passed = .false. 
  end if 

  call tst_2str%finalize()
  ! ----------------------------------------------------------------------------
  print *, "  Incrementing nstr"
  !
  ! Increment nstr 
  !
  call make_copy_nstr
  call increment_with_1scl(tst_nstr)
  if(.not. ops_match(tst_nstr, ref_nstr)) then 
    call report_err("nstr+1scl fails")
    passed = .false. 
  end if 

  call make_copy_nstr
  call increment_with_2str(tst_nstr)
  if(.not. ops_match(tst_nstr, ref_nstr)) then 
    call report_err("nstr+2str fails")
    passed = .false. 
  end if 

  call make_copy_nstr
  call increment_with_nstr(tst_nstr)
  if(.not. ops_match(tst_nstr, ref_nstr)) then 
    call report_err("nstr+nstr fails")
    passed = .false. 
  end if 

  call tst_2str%finalize()
  ! ----------------------------------------------------------------------------
  print *, "  Halving/doubling optical thickness"
  !
  ! Adding two media of half optical thickness to recover original values 
  !
  call make_copy_1scl
  tst_1scl%tau = 0.5_wp * tst_1scl%tau 
  call stop_on_err(tst_1scl%increment(tst_1scl))
  if(.not. ops_match(tst_1scl, ref_1scl)) then 
    call report_err("1scl half/double fails")
    passed = .false. 
  end if 

  call make_copy_2str
  tst_2str%tau = 0.5_wp * tst_2str%tau 
  call stop_on_err(tst_2str%increment(tst_2str))
  if(.not. ops_match(tst_2str, ref_2str)) then 
    call report_err("2str half/double fails")
    passed = .false. 
  end if 

  call make_copy_nstr
  tst_nstr%tau = 0.5_wp * tst_nstr%tau 
  call stop_on_err(tst_nstr%increment(tst_nstr))
  if(.not. ops_match(tst_nstr, ref_nstr)) then 
    call report_err("nstr half/double fails")
    passed = .false. 
  end if 
  ! ----------------------------------------------------------------------------
  print *, "  Delta scaling"
  !
  ! Delta-scale with forward-fraction f=0 (i.e. Rayleigh scattering)
  !
  call make_copy_2str
  call stop_on_err(tst_2str%delta_scale(spread(spread(spread(0._wp, 1, ncol), 2, nlay), 3, 1)))
  if(.not. ops_match(tst_2str, ref_2str)) then 
    call report_err("2str half/double fails")
    passed = .false. 
  end if 
  ! ----------------------------------------------------------------------------
  if (.not. passed) call stop_on_err("Optical props unit tests fail")
  print *, "Optical properties unit testing finished"
  print *
  ! ----------------------------------------------------------------------------
contains 
  ! ----------------------------------------------------------------------------
  !
  ! Make copies of the existing optical depth 
  !
  subroutine make_copy_1scl 
    call tst_1scl%finalize()
    call stop_on_err(tst_1scl%alloc_1scl(ref_1scl%get_ncol(), ref_1scl%get_nlay(), & 
                    ref_1scl))
    tst_1scl%tau = ref_1scl%tau
  end subroutine make_copy_1scl 
  ! ----------------------------------------------------------------------------
  !
  ! Make copies of the existing optical depth 
  !
  subroutine make_copy_2str
    call tst_2str%finalize()
    call stop_on_err(tst_2str%alloc_2str(ref_2str%get_ncol(), ref_2str%get_nlay(), & 
                     ref_2str))
    tst_2str%tau = ref_2str%tau
    tst_2str%ssa = ref_2str%ssa 
    tst_2str%g   = ref_2str%g 
  end subroutine make_copy_2str 
  ! ----------------------------------------------------------------------------
  !
  ! Make copies of the existing optical depth 
  !
  subroutine make_copy_nstr
    call tst_nstr%finalize()
    call stop_on_err(tst_nstr%alloc_nstr(ref_nstr%get_nmom(), ref_nstr%get_ncol(), ref_nstr%get_nlay(), & 
                     ref_2str))
    tst_nstr%tau = ref_nstr%tau
    tst_nstr%ssa = ref_nstr%ssa 
    tst_nstr%p   = ref_nstr%p 
  end subroutine make_copy_nstr 
  ! ----------------------------------------------------------------------------
end program optical_prop_unit_tests