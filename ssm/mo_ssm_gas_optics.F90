! See documentation in other modules... 
!
! Contacts: could provide here... names, emails or web site 
!
! Copyright 2025-,  ... Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Class implementing Simple Spectral model gas optics
!>
!> Details here 
! -------------------------------------------------------------------------------------------------
module mo_gas_optics_ssm 
  use mo_rte_kind,           only: wp, wl
  use mo_rte_config,         only: check_extents, check_values
  use mo_rte_util_array,     only: zero_array
  use mo_rte_util_array_validation, &
                             only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props,      only: ty_optical_props
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_optics_util_string, only: lower_case, string_in_array, string_loc_in_array
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_optics,         only: ty_gas_optics
  implicit none
  private
  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_gas_optics), public :: ty_gas_optics_ssm
    private
    ! list of parameters - ideally as small as possible, and ideally extensible 
    !   maybe? a derived type with a gas name, N slopes/peak locations, 1 continuum? 
    !
    ! Need to be able to distinguish LW and SW data
    ! Need to provide spectral discretization (number of bands and gpoints, which will be the same; starting and ending band wavenumbers)
    contains 
      procedure, public :: configure
      procedure, public :: source_is_internal
      procedure, public :: source_is_external
      procedure, public :: gas_optics_int
      procedure, public :: gas_optics_ext
      procedure, public :: get_press_min
      procedure, public :: get_press_max
      procedure, public :: get_temp_min
      procedure, public :: get_temp_max

  end type ty_gas_optics_ssm

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Configure the simple spectral model parameters 
  !
  function configure(this) result(error_msg)
    !
    ! All the parameters for the SSM need to get added to the argument list 
    !   At least one of the arguments needs to distinguish LW from SW
    !   Probably need to specify the spectral discretization during configuration 
    ! 
    class(ty_gas_optics_ssm), intent(inout) :: this
    character(len=128)                      :: error_msg     !! Empty if succssful
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    error_msg = ""
  end function configure
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Compute gas optical depth and Planck source functions,
  !>  given temperature, pressure, and composition
  !
  function gas_optics_int(this,                             &
                          play, plev, tlay, tsfc, gas_desc, &
                          optical_props, sources,           &
                          col_dry, tlev) result(error_msg)
    ! inputs
    class(ty_gas_optics_ssm), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in   ) :: tsfc      !! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props !! Optical properties
    class(ty_source_func_lw    ),  &
                              intent(inout) :: sources       !! Planck sources
    character(len=128)                      :: error_msg     !! Empty if succssful
    ! Optional inputs
    real(wp), dimension(:,:),   intent(in   ), &
                           optional, target :: col_dry, &  !! Column dry amount; dim(ncol,nlay)
                                               tlev        !! level temperatures [K]; (ncol,nlay+1)
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    error_msg = ""
    
    ! How much data validation to implement? 

    
  end function gas_optics_int

  !------------------------------------------------------------------------------------------
  !
  !> Compute gas optical depth given temperature, pressure, and composition
  !>    Top-of-atmosphere stellar insolation is also reported
  !
  function gas_optics_ext(this,                         &
                          play, plev, tlay, gas_desc,   & ! mandatory inputs
                          optical_props, toa_src,       & ! mandatory outputs
                          col_dry) result(error_msg)      ! optional input

    class(ty_gas_optics_ssm), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props
    real(wp), dimension(:,:), intent(  out) :: toa_src     !! Incoming solar irradiance(ncol,ngpt)
    character(len=128)                      :: error_msg   !! Empty if successful

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    error_msg = ""
  end function gas_optics_ext

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for internal sources/longwave, false otherwise
  !
  pure function source_is_internal(this)
    class(ty_gas_optics_ssm), intent(in) :: this
    logical                              :: source_is_internal
    source_is_internal = .true.
  end function source_is_internal
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if configured for external sources/shortwave, false otherwise
  !
  pure function source_is_external(this)
    class(ty_gas_optics_ssm), intent(in) :: this
    logical                              :: source_is_external
    source_is_external = .false. 
  end function source_is_external
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Boilerplate functions - not relevant for SSM
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Minimum valid pressure
  !
  pure function get_press_min(this)
    class(ty_gas_optics_ssm), intent(in) :: this
    real(wp)                             :: get_press_min !! minimum valid pressure

    get_press_min = 0._wp
  end function get_press_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Maximum valid pressure
  !
  pure function get_press_max(this)
    class(ty_gas_optics_ssm), intent(in) :: this
    real(wp)                                :: get_press_max !! maximum valid pressure

    get_press_max = huge(1._wp)
  end function get_press_max

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Minimum valid temperature 
  !
  pure function get_temp_min(this)
    class(ty_gas_optics_ssm), intent(in) :: this
    real(wp)                             :: get_temp_min !! minimum valid temperature 
    get_temp_min = 0._wp
  end function get_temp_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Maximum valid pressure
  !
  pure function get_temp_max(this)
    class(ty_gas_optics_ssm), intent(in) :: this
    real(wp)                             :: get_temp_max !! maximum valid pressure

    get_temp_max = huge(1._wp)
  end function get_temp_max
  !--------------------------------------------------------------------------------------------------------------------

end module mo_gas_optics_ssm