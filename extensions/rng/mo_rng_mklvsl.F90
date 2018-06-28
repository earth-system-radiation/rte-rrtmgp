! Module: mo_rng_mklvsl

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

!
! Implementation of random number interface using Intel Vector Statistics Library
!

include 'mkl_vsl.f90' ! Needed to generate MKL_VSL modules below
module mo_rng_mklvsl
  use mo_rte_kind,   only : dp
  use mo_rrtmgp_errors, only: write_message
  use mo_rng
  USE MKL_VSL_TYPE
  USE MKL_VSL
  implicit none

  integer, parameter :: rng_type = VSL_BRNG_MT19937             ! Merseene Twister
  ! Alternatives are VSL_BRNG_SFMT19937, maybe VSL_BRNG_MT2203?
  integer, parameter :: rng_method = VSL_RNG_METHOD_UNIFORM_STD ! Uniform distribution

  ! extended random number generator class
  type, extends(ty_rng), public :: ty_rng_mklvsl
    TYPE (VSL_STREAM_STATE) :: stream
  contains
    procedure, private :: get_random_vec => get_rng_mkl_vec
    procedure, private :: get_random_vec_mask &
                                         => get_rng_mkl_mask
    procedure, private :: init_rng       => init_rng_mkl
    procedure, private :: end_rng        => end_rng_mkl
  end type
contains
  ! -------------------------------------------------------------------------------------
  ! Provide num random numbers following a uniform distribution between 0 and 1
  !
  function get_rng_mkl_vec(this, num)
    class(ty_rng_mklvsl) :: this
    integer,  intent(in   ) :: num
    real(DP), dimension(num) :: get_rng_mkl_vec

    integer :: status

    status = vdrnguniform(rng_method, this%stream, num, get_rng_mkl_vec, 0._dp, 1._dp)
    if(status /= VSL_STATUS_OK) call write_message("Error getting random numbers")
  end function get_rng_mkl_vec
  ! -------------------------------------------------------------------------------------
  ! Provide random numbers for the TRUE elements of MASK
  !
  function get_rng_mkl_mask(this, mask)
    class(ty_rng_mklvsl)               :: this
    logical, dimension(:), intent(in) :: mask
    real(DP), dimension(size(mask)) :: get_rng_mkl_mask

    get_rng_mkl_mask(:) = UNPACK(get_rng_mkl_vec(this, COUNT(mask)), MASK = mask, FIELD = 0._dp)
  end function get_rng_mkl_mask
  ! -------------------------------------------------------------------------------------
  ! Initialize the random number state.
  !
  subroutine init_rng_mkl(this, seeds)
    class(ty_rng_mklvsl) :: this
    integer,  dimension(:), intent(in) :: seeds

    integer :: status
    status = vslnewstream(this%stream, rng_type, seeds(1))
    if(status /= VSL_STATUS_OK) call write_message("Error initializing random number stream")
  end subroutine init_rng_mkl
  ! -------------------------------------------------------------------------------------
  ! Release any resources associated with the RNG.
  !
  subroutine end_rng_mkl(this)
    class(ty_rng_mklvsl) :: this

    integer :: status
    status = vsldeletestream(this%stream)
    if(status /= VSL_STATUS_OK) call write_message("Error finalizing random number stream")
  end subroutine end_rng_mkl
  ! -------------------------------------------------------------------------------------
end module mo_rng_mklvsl
