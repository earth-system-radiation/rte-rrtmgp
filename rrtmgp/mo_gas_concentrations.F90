! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Encapsulates a collection of volume mixing ratios (concentrations) of gases.
!   Each concentration is associated with a name, normally the chemical formula.
!
! Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
!   (nlay and ncol are determined from the input arrays; self-consistency is enforced)
!   example:
!   error_msg = gas_concs%set_vmr('h2o', values(:,:))
!   error_msg = gas_concs%set_vmr('o3' , values(:)  )
!   error_msg = gas_concs%set_vmr('co2', value      )
!
! Values can be requested as profiles (valid only if there are no 2D fields present in the object)
!   or as 2D fields. Values for all columns are returned although the entire collection
!   can be subsetted in the column dimension
!
! Subsets can be extracted in the column dimension
!
! Functions return strings. Non-empty strings indicate an error.
!
! -------------------------------------------------------------------------------------------------

module mo_gas_concentrations
  use mo_rte_kind,           only: wp
  use mo_rrtmgp_util_string, only: lower_case
  use mo_rte_util_array,     only: any_vals_outside
  implicit none
  integer, parameter :: GAS_NOT_IN_LIST = -1

  type, private :: conc_field
    real(wp), dimension(:,:), pointer :: conc => NULL()
  end type conc_field

  type, public :: ty_gas_concs
    !
    ! Data
    !
    character(len=32), dimension(:), allocatable :: gas_name
    type(conc_field),  dimension(:), allocatable :: concs
    integer :: ncol = 0, nlay = 0
    contains
      !
      ! Procedures
      !
      procedure, private :: find_gas
      procedure, private :: set_vmr_scalar
      procedure, private :: set_vmr_1d
      procedure, private :: set_vmr_2d
      procedure, private :: get_vmr_1d
      procedure, private :: get_vmr_2d
      procedure, private :: get_subset_range
      final :: del
      !
      ! public interface
      !
      procedure, public :: init
      procedure, public :: reset
      generic,   public :: set_vmr => set_vmr_scalar, &
                                      set_vmr_1d, &
                                      set_vmr_2d
      generic,   public :: get_vmr => get_vmr_1d, &
                                      get_vmr_2d
      generic,   public :: get_subset => get_subset_range
      procedure, public :: get_num_gases
      procedure, public :: get_gas_names
  end type ty_gas_concs
contains
  ! -------------------------------------------------------------------------------------
  function init(this, gas_names) result(error_msg)
    class(ty_gas_concs),            intent(inout) :: this
    character(len=*), dimension(:), intent(in   ) :: gas_names
    character(len=128)                            :: error_msg
    ! ---------
    integer :: i, j, ngas
    ! ---------
    error_msg = ''
    ngas = size(gas_names)
    !
    ! Check for no duplicate gas names, no empty names
    !
    if(any(len_trim(gas_names) == 0)) &
      error_msg = "ty_gas_concs%init(): must provide non-empty gas names"

    do i = 1, ngas-1
      do j = i+1, ngas
        if (lower_case(trim(gas_names(i))) == lower_case(trim(gas_names(j)))) then
          error_msg = "ty_gas_concs%init(): duplicate gas names aren't allowed"
          exit
        end if
      end do
    end do
    if(error_msg /= "") return
    !
    ! Allocate fixed-size arrays
    !
    call this%reset()
    allocate(this%gas_name(ngas), this%concs(ngas))
    !$acc enter data copyin(this)
    !$acc enter data copyin(this%concs)

    this%gas_name(:) = gas_names(:)
  end function
  ! -------------------------------------------------------------------------------------
  !
  ! Set concentrations --- scalar, 1D, 2D
  !
  ! -------------------------------------------------------------------------------------
  function set_vmr_scalar(this, gas, w) result(error_msg)
    ! In OpenACC context scalar w always assumed to be on the CPU
    class(ty_gas_concs), intent(inout) :: this
    character(len=*),    intent(in   ) :: gas
    real(wp),            intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''
    if (w < 0._wp .or. w > 1._wp) then
      error_msg = 'ty_gas_concs%set_vmr(): concentrations should be >= 0, <= 1'
      return
    endif

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%set_vmr(): trying to set ' // trim(gas) // ' but name not provided at initialization'
      return
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    ! This cannot be made a function, because we need all the hierarchy for the correct OpenACC attach
    if (associated(this%concs(igas)%conc)) then
      if ( any(shape(this%concs(igas)%conc) /= [1, 1]) ) then
        !$acc exit data delete(this%concs(igas)%conc)
        deallocate(this%concs(igas)%conc)
        nullify   (this%concs(igas)%conc)
      end if
    end if
    if (.not. associated(this%concs(igas)%conc)) then
      allocate(this%concs(igas)%conc(1,1))
      !$acc enter data create(this%concs(igas)%conc)
    end if

    !$acc kernels
    this%concs(igas)%conc(:,:) = w
    !$acc end kernels
  end function set_vmr_scalar
  ! -------------------------------------------------------------------------------------
  function set_vmr_1d(this, gas, w) result(error_msg)
    ! In OpenACC context w assumed to be either on the CPU or on the GPU
    class(ty_gas_concs), intent(inout) :: this
    character(len=*),    intent(in   ) :: gas
    real(wp), dimension(:), &
                         intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''

    if (any_vals_outside(w, 0._wp, 1._wp)) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1'
    endif
    if(this%nlay > 0) then
      if(size(w) /= this%nlay) error_msg = 'ty_gas_concs%set_vmr: different dimension (nlay)'
    else
      this%nlay = size(w)
    end if
    if(error_msg /= "") return

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%set_vmr(): trying to set ' // trim(gas) // ' but name not provided at initialization'
      return
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    ! This cannot be made a function, because we need all the hierarchy for the correct OpenACC attach
    if (associated(this%concs(igas)%conc)) then
      if ( any(shape(this%concs(igas)%conc) /= [1, this%nlay]) ) then
        !$acc exit data delete(this%concs(igas)%conc)
        deallocate(this%concs(igas)%conc)
        nullify   (this%concs(igas)%conc)
      end if
    end if
    if (.not. associated(this%concs(igas)%conc)) then
      allocate(this%concs(igas)%conc(1,this%nlay))
      !$acc enter data create(this%concs(igas)%conc)
    end if

    !$acc kernels copyin(w)
    this%concs(igas)%conc(1,:) = w
    !$acc end kernels

    !$acc exit data delete(w)
  end function set_vmr_1d
  ! -------------------------------------------------------------------------------------
  function set_vmr_2d(this, gas, w) result(error_msg)
    ! In OpenACC context w assumed to be either on the CPU or on the GPU
    class(ty_gas_concs), intent(inout) :: this
    character(len=*),    intent(in   ) :: gas
    real(wp), dimension(:,:),  &
                         intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''

    if (any_vals_outside(w, 0._wp, 1._wp)) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1'
    endif

    if(this%ncol > 0 .and. size(w, 1) /= this%ncol) then
      error_msg = 'ty_gas_concs%set_vmr: different dimension (ncol)'
    else
      this%ncol = size(w, 1)
    end if

    if(this%nlay > 0 .and. size(w, 2) /= this%nlay) then
      error_msg = 'ty_gas_concs%set_vmr: different dimension (nlay)'
    else
      this%nlay = size(w, 2)
    end if
    if(error_msg /= "") return

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%set_vmr(): trying to set ' // trim(gas) // 'but name not provided at initialization'
      return
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    ! This cannot be made a function, because we need all the hierarchy for the correct OpenACC attach
    if (associated(this%concs(igas)%conc)) then
      if ( any(shape(this%concs(igas)%conc) /= [this%ncol,this%nlay]) ) then
        !$acc exit data delete(this%concs(igas)%conc)
        deallocate(this%concs(igas)%conc)
        nullify   (this%concs(igas)%conc)
      end if
    end if
    if (.not. associated(this%concs(igas)%conc)) then
      allocate(this%concs(igas)%conc(this%ncol,this%nlay))
      !$acc enter data create(this%concs(igas)%conc)
    end if

    !$acc kernels copyin(w)
    this%concs(igas)%conc(:,:) = w(:,:)
    !$acc end kernels
  end function set_vmr_2d
  ! -------------------------------------------------------------------------------------
  !
  ! Return volume mixing ratio as 1D or 2D array
  !
  ! -------------------------------------------------------------------------------------
  !
  ! 1D array ( lay depdendence only)
  !
  function get_vmr_1d(this, gas, array) result(error_msg)
    class(ty_gas_concs) :: this
    character(len=*),         intent(in ) :: gas
    real(wp), dimension(:),   intent(out) :: array
    character(len=128) :: error_msg
    ! ---------------------
    integer :: igas
    ! ---------------------
    error_msg = ''

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' not found'
    else if(.not. associated(this%concs(igas)%conc)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // " concentration hasn't been set"
    else if(size(this%concs(igas)%conc, 1) > 1) then ! Are we requesting a single profile when many are present?
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' requesting single profile but many are available'
    end if

    if(this%nlay > 0 .and. this%nlay /= size(array)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (nlay)'
    end if
    if(error_msg /= "") return

    !$acc data copyout (array) present(this)
    if(size(this%concs(igas)%conc, 2) > 1) then
      !$acc kernels default(none)
      array(:) = this%concs(igas)%conc(1,:)
      !$acc end kernels
    else
      !$acc kernels default(none)
      array(:) = this%concs(igas)%conc(1,1)
      !$acc end kernels
    end if
    !$acc end data

  end function get_vmr_1d
  ! -------------------------------------------------------------------------------------
  !
  ! 2D array (col, lay)
  !
  function get_vmr_2d(this, gas, array) result(error_msg)
    class(ty_gas_concs) :: this
    character(len=*),         intent(in ) :: gas
    real(wp), dimension(:,:), intent(out) :: array
    character(len=128)                    :: error_msg
    ! ---------------------
    integer :: icol, ilay, igas
    ! ---------------------
    error_msg = ''

    igas = this%find_gas(gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' not found'
    else if(.not. associated(this%concs(igas)%conc)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // " concentration hasn't been set"
    end if
    !
    ! Is the requested array the correct size?
    !
    if(this%ncol > 0 .and. this%ncol /= size(array,1)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (ncol)'
    end if
    if(this%nlay > 0 .and. this%nlay /= size(array,2)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (nlay)'
    end if
    if(error_msg /= "") return

    !$acc data copyout (array) present(this, this%concs)
    if(size(this%concs(igas)%conc, 1) > 1) then      ! Concentration stored as 2D
      !$acc parallel loop collapse(2) default(none)
      do ilay = 1, size(array,2)
        do icol = 1, size(array,1)
          !print *, (size(this%concs))
          array(icol,ilay) = this%concs(igas)%conc(icol,ilay)
        end do
      end do
    else if(size(this%concs(igas)%conc, 2) > 1) then ! Concentration stored as 1D
      !$acc parallel loop collapse(2) default(none)
      do ilay = 1, size(array,2)
        do icol = 1, size(array,1)
         array(icol, ilay) = this%concs(igas)%conc(1,ilay)
        end do
      end do
    else                                             ! Concentration stored as scalar
      !$acc parallel loop collapse(2) default(none)
      do ilay = 1, size(array,2)
        do icol = 1, size(array,1)
          array(icol,ilay) = this%concs(igas)%conc(1,1)
        end do
      end do
    end if
    !$acc end data

  end function get_vmr_2d
  ! -------------------------------------------------------------------------------------
  !
  ! Extract a subset of n columns starting with column 'start'
  !
  ! -------------------------------------------------------------------------------------
  function get_subset_range(this, start, n, subset) result(error_msg)
    class(ty_gas_concs),      intent(in   ) :: this
    integer,                  intent(in   ) :: start, n
    class(ty_gas_concs),      intent(inout) :: subset
    character(len=128)                      :: error_msg
    ! ---------------------
    integer :: i
    ! ---------------------
    error_msg = ''
    if(n <= 0) &
       error_msg = "gas_concs%get_vmr: Asking for 0 or fewer columns "
    if(start < 1 ) &
       error_msg = "gas_concs%get_vmr: Asking for columns outside range"
    if(this%ncol > 0 .and. start > this%ncol .or. start+n-1 > this%ncol ) &
       error_msg = "gas_concs%get_vmr: Asking for columns outside range"
    if(error_msg /= "") return

    call subset%reset()
    allocate(subset%gas_name(size(this%gas_name)), &
             subset%concs   (size(this%concs))) ! These two arrays should be the same length
    !$acc enter data create(subset, subset%concs)
    subset%nlay = this%nlay
    subset%ncol = merge(n, 0, this%ncol > 0)
    subset%gas_name(:)  = this%gas_name(:)

    do i = 1, size(this%gas_name)
      !
      ! Preserve scalar/1D/2D representation in subset,
      !   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
      !
      allocate(subset%concs(i)%conc(min(max(subset%ncol,1), size(this%concs(i)%conc, 1)), &
                                    min(    subset%nlay,    size(this%concs(i)%conc, 2))))
      !$acc enter data create(subset%concs(i)%conc)
      if(size(this%concs(i)%conc, 1) > 1) then      ! Concentration stored as 2D
        !$acc kernels
        subset%concs(i)%conc(:,:) = this%concs(i)%conc(start:(start+n-1),:)
        !$acc end kernels
      else
        !$acc kernels
        subset%concs(i)%conc(:,:) = this%concs(i)%conc(:,:)
        !$acc end kernels
      end if
    end do

  end function get_subset_range
  ! -------------------------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  ! -------------------------------------------------------------------------------------
  subroutine reset(this)
    class(ty_gas_concs), intent(inout) :: this
    ! -----------------
    integer :: i
    ! -----------------
    this%nlay = 0
    this%ncol = 0
    if(allocated(this%gas_name)) deallocate(this%gas_name)
    if (allocated(this%concs)) then
      do i = 1, size(this%concs)
        if(associated(this%concs(i)%conc)) then
          !$acc exit data delete(this%concs(i)%conc)
          deallocate(this%concs(i)%conc)
          nullify(this%concs(i)%conc)
        end if
      end do
      !$acc exit data delete(this%concs)
      deallocate(this%concs)
    end if
  end subroutine reset
  ! -------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  ! -------------------------------------------------------------------------------------
  pure function get_num_gases(this)
    class(ty_gas_concs), intent(in) :: this
    integer :: get_num_gases

    get_num_gases = size(this%gas_name)
    return
  end function get_num_gases
  ! -------------------------------------------------------------------------------------
  pure function get_gas_names(this)
    class(ty_gas_concs), intent(in) :: this
    character(len=32), dimension(this%get_num_gases()) :: get_gas_names

    get_gas_names(:) = this%gas_name(:)
    return
  end function get_gas_names
  ! -------------------------------------------------------------------------------------
  !
  ! Private procedures
  !
  ! -------------------------------------------------------------------------------------
  !
  ! find gas in list; GAS_NOT_IN_LIST if not found
  !
  function find_gas(this, gas)
    character(len=*),    intent(in) :: gas
    class(ty_gas_concs), intent(in) :: this
    integer                         :: find_gas
    ! -----------------
    integer :: igas
    ! -----------------
    find_gas = GAS_NOT_IN_LIST
    if(.not. allocated(this%gas_name)) return
    do igas = 1, size(this%gas_name)
      if (lower_case(trim(this%gas_name(igas))) == lower_case(trim(gas))) then
        find_gas = igas
      end if
    end do
  end function
  ! -------------------------------------------------------------------------------------
  subroutine del(this)
    type(ty_gas_concs), intent(inout) :: this
    call this%reset()
    !$acc exit data delete(this)
  end subroutine del
  ! -------------------------------------------------------------------------------------
end module mo_gas_concentrations
