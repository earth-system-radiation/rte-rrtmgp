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
module mo_optics_ssm
  use mo_rte_kind,           only: wp, wl
  use mo_rte_config,         only: check_extents, check_values
  use mo_rte_util_array,     only: zero_array
  use mo_rte_util_array_validation, &
                             only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props,      only: ty_optical_props
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_optics_util_string, only: lower_case, string_in_array, string_loc_in_array
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_constants,   only: avogad, m_dry, m_h2o, grav
  use mo_optics_ssm_kernels, only: compute_tau, compute_Planck_source

  implicit none
  private

  real(wp), parameter, public :: Tsun_ssm = 5760._wp ! Default sun temperature for SSM
  real(wp), parameter, public :: tsi = 1360._wp ! Default total solar irradiance
  real(wp), parameter, public :: mw_h2o = 0.018_wp ! Molecular weight h2o
  real(wp), parameter, public :: mw_co2 = 0.044_wp ! Molecular weight co2
  real(wp), parameter, public :: mw_o3  = 0.048_wp ! Molecular weight o3
  !
  ! Do the other SSM defaults - absorption parameters, spectral dicretization -
  !   get declared here as public entities? Or do we add general introscpection?
  !   (General introspection won't let us recover the triangles)
  ! If used before configured do we just call this%configure() using defaults?
  !
  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_gas_optics), public :: ty_optics_ssm
    private
    character(32), dimension(:),    allocatable :: gas_names 
    real(wp),      dimension(:),    allocatable :: mol_weights
    real(wp),      dimension(:, :), allocatable :: absorption_coeffs
    real(wp),      dimension(:),    allocatable :: nus, dnus
      ! total absorption coefficients at spectral points (nnus, ngases)
    real(wp) :: Tstar  = 0._wp, &
                pref   = 500._wp * 500._wp, & ! 500 hPa
                m_dry  = 0.029_wp ! molecular weight of dry air [kg/mol]
    contains
      procedure, public :: configure
      procedure, public :: gas_optics_int
      procedure, public :: gas_optics_ext
      procedure, public :: source_is_internal
      procedure, public :: source_is_external
      procedure, public :: get_press_min
      procedure, public :: get_press_max
      procedure, public :: get_temp_min
      procedure, public :: get_temp_max

  end type ty_optics_ssm

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Configure the simple spectral model parameters
  !> Gas optics for simple spectral models are described as N triangles of ln(kappa) for
  !>   each of M gases.
  !> Each triangle id defined by
  !>    a central wavenumber nu0
  !>    the absorption coefficient at that wavenumber
  !>    the slope d ln(kappa)/d nu
  !> By default there are two gases (h2o and co2).
  !>   h2o has three triangles (rotational lines, and a continuum with infinite slope)
  !>   co2 has a single triangle
  !> Absorption coefficents are defined at pref and pressure-broadened as p/pref
  !> Spectral discretization is defined as
  !>   max and min wavenumbers and
  !>   a set of wavenumbers at which to evaluate absorption and Planck function
  !> If Tstar is provided the gas optics are for insolation
  !
  function configure(this,                       &
                     gas_names, triangle_params, &
                     nus, nu_min, nu_max,        &
                     Tstar) result(error_msg)
    !
    ! All the parameters for the SSM need to get added to the argument list
    !   At least one of the arguments needs to distinguish LW from SW
    !   Probably need to specify the spectral discretization during configuration
    !
    class(ty_optics_ssm),      intent(inout) :: this
    character(32), dimension(:),   intent(in   ) :: gas_names
    real(wp),      dimension(:,:), intent(in   ) :: triangle_params
      !! (ntriangles, 4) where the second dimension holds [gas_index, kappa0, nu0, l]
    real(wp),                      intent(in   ) ::  kappa_cld_lw
      !! cloud optical properties in longwave
    real(wp),                      intent(in   ) ::  kappa_cld_sw
    real(wp),                      intent(in   ) ::  g_cld_sw
    real(wp),                      intent(in   ) ::  ssa_cld_sw
      !! cloud optical properties in shortwave
    real(wp),      dimension(:),   intent(in   ) :: nus
      !! Wavenumbers at which to evaluate Planck function and absorption coefficient
    real(wp),                      intent(in   ) :: nu_min, nu_max
      !! Upper and lower bounds of spectrum
    real(wp),      optional,       intent(in   ) :: Tstar
      !! Temperature for stellar insolation
    real(wp),      optional,       intent(in   ) :: tsi
      !! Total solar irradiance
    character(len=128)                      :: error_msg     !! Empty if successful
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    integer :: nnu, ngas
    integer :: inu, itri, igas

    error_msg = ""
    ngas = size(gas_names)
    nnu  = size(nus)

    ! Input sanitizing?
    ! triangle params: index <= ngases, kappa0 >= 0; nu_min < nu0s < nu_max; l > 0
    ! nus > 0; ascending? nu_min <= nus <= max_nu
    ! Tstar > 0 if specified
    if (.not. all(nus > 0.0_wp)) then
      error_msg = "ssm_gas_optics(): all wavenumbers must be > 0"
    end if

    if (present(Tstar)) then
      if (.not. all(Tstar > 0.0_wp)) then
        error_msg = "ssm_gas_optics(): Tstar must be > 0"
      end if
    end if
    
    if (.not. all(triangle_params(:, 2) >= 0.0_wp)) then
      error_msg = "ssm_gas_optics(): kappa0 needs to be >=0"
    end if
    
    if (.not. all(nu_min < triangle_params(:, 3) < nu_max)) then
      error_msg = "ssm_gas_optics(): nu0 needs to be less than nu_max and greater than nu_min"
    end if

    if (.not. all(triangle_params(:, 4) > 0.0_wp)) then
      error_msg = "ssm_gas_optics(): l needs to be > 0"
    end if
    if (error_msg /= '') return
    
    if(allocated(this%gas_names)) &
      deallocate(this%gas_names, this%mol_weights, this%nus, this%dnus, this%absorption_coeffs)

    allocate(this%gas_names(ngas), &
             this%mol_weights(ngas), &
             this%nus (nnu), &
             this%dnus(nnu), &
             this%absorption_coeffs(ngas,nnu) )

    this%nus(1:nnu) = nus(1:nnu)
    this%gas_names(:) = gas_names(:)

    ! Set molar masses based on gas names
    !   Maybe this is better as module data...
    do igas = 1, ngas
      select case (trim(lower_case(gas_names(igas))))
        case ('h2o')
          this%mol_weights(igas) = mw_h2o
        case ('co2')
          this%mol_weights(igas) = mw_co2
        case ('o3')
          this%mol_weights(igas) = mw_o3
        case default
          error_msg = "Unknown gas: " // trim(gas_names(igas))
          return
      end select
    end do
    if (error_msg /= '') return

    ! Spectral discretization - edges of "bands"
    ! err_message = this%ty_optical_props%init(band_lims_wavenum, band2gpt)
    ! where band2gpt = [(i, i = 1, nnu)]

    ! Compute absorption coefficients by summing exponentials at each nu
    ! Initialize absorption coefficients to zero
    this%absorption_coeffs(:,:) = 0._wp

    ! robert, can we parallelize this?
    do itri = 1, size(triangle_params, 1)
      do inu = 1, nnu
        this%absorption_coeffs(nint(triangle_params(itri, 1)), inu) = &
          this%absorption_coeffs(nint(triangle_params(itri, 1)), inu) + &
          triangle_params(itri, 2) * exp(-abs(this%nus(inu) - triangle_params(itri, 3)) / triangle_params(itri, 4))
      end do
    end do

    if(present(Tstar)) this%Tstar = Tstar
    if(present(tsi))   this%tsi = tsi
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
    class(ty_optics_ssm),     intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa]; (ncol,nlay+1)
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
                           optional, target :: col_dry, &  !! Column dry amount; dim(ncol,nlay), will be ignored
                                               tlev        !! level temperatures [K]; (ncol,nlay+1)
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    integer :: ncol, nlay, nnu, ngas
    integer :: igas, icol, ilay, idx_gas
    real(wp), dimension(size(this%gas_names), size(play,1), size(play,2)) :: layer_mass
    error_msg = ""

    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = size(this%nus)
    ngas = size(this%gas_names)

    ! How much data validation to implement?

    call compute_layer_mass(ncol, nlay, ngas, &
                            this, plev, gas_desc, layer_mass, &
                            error_msg)
    if (error_msg /= '') return

    !
    ! Absorption optical depth
    !
    call compute_tau(ncol, nlay, nnu, ngas,   &
                     this%absorption_coeffs, play, this%pref, layer_mass, &
                     optical_props%tau)

    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    !
    ! Planck function sources
    ! Longwave: sources%sfc_source, sources%lay_source, sources%lev_source, &
    !            sources%sfc_source_Jac
    ! How to compute 1D/sfc - another kernel subroutine?
    ! How to compute Tlev if not provided - make generic?
    call compute_Planck_source_2D(ncol, nlay,   nnu, &
                               this%nus, this%dnus, tlay,   &
                               sources%lay_source)
    call compute_Planck_source_2D(ncol, nlay+1, nnu, &
                               this%nus, this%dnus, tlev,   &
                               sources%lev_source)
    call compute_Planck_source_1D(ncol, nnu, &
                               this%nus, this%dnus, tsfc,   &
                               sources%sfc_source)

    call zero_array(ncol, nnu, sources%sfc_source_Jac)
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

    class(ty_optics_ssm),     intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props
    real(wp), dimension(:),   intent(  out) :: toa_src     !! Incoming solar irradiance(ncol)
    character(len=128)                      :: error_msg   !! Empty if successful

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    integer :: ncol, nlay, nnu, ngas
    integer :: igas, icol, ilay, idx_gas
    real(wp), dimension(size(this%gas_names), size(play,1), size(play,2)) :: layer_mass
    error_msg = ""

    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = size(this%nus)
    ngas = size(this%gas_names)

    ! How much data validation to implement?

    call compute_layer_mass(ncol, nlay, ngas, &
                            this, plev, gas_desc, layer_mass, &
                            error_msg)
    if (error_msg /= '') return

    !
    ! Absorption optical depth
    !
    call compute_tau(ncol, nlay, nnu, ngas,   &
                     this%absorption_coeffs, play, this%pref, layer_mass, &
                     optical_props%tau)
                     
    !
    ! Planck function sources
    ! Shortwave: incoming solar irradiance
    call compute_Planck_source_1D(ncol, nlay,   nnu, &
                               this%nus, this%dnus, this%Tstar,   &
                               toa_src)
    
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
    class(ty_optics_ssm), intent(in) :: this
    logical                          :: source_is_internal
    source_is_internal = (this%Tstar <= 0._wp)
  end function source_is_internal
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if configured for external sources/shortwave, false otherwise
  !
  pure function source_is_external(this)
    class(ty_optics_ssm), intent(in) :: this
    logical                              :: source_is_external
    source_is_external = (this%Tstar > 0._wp)
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
    class(ty_optics_ssm), intent(in) :: this
    real(wp)                         :: get_press_min !! minimum valid pressure

    get_press_min = 0._wp
  end function get_press_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Maximum valid pressure
  !
  pure function get_press_max(this)
    class(ty_optics_ssm), intent(in) :: this
    real(wp)                         :: get_press_max !! maximum valid pressure

    get_press_max = huge(1._wp)
  end function get_press_max

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Minimum valid temperature
  !
  pure function get_temp_min(this)
    class(ty_optics_ssm), intent(in) :: this
    real(wp)                         :: get_temp_min !! minimum valid temperature
    get_temp_min = 0._wp
  end function get_temp_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Maximum valid pressure
  !
  pure function get_temp_max(this)
    class(ty_optics_ssm), intent(in) :: this
    real(wp)                         :: get_temp_max !! maximum valid pressure

    get_temp_max = huge(1._wp)
  end function get_temp_max
  
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Compute layer masses from gas concentrations and pressure levels 
  !
  subroutine compute_layer_mass(ncol, nlay, ngas, this, plev, gas_desc, layer_mass, error_msg)
    integer,  intent(in ) :: ncol, nlay, ngas
    class(ty_optics_ssm),     intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: plev       !! level pressures [Pa]; (ncol,nlay+1)
    type(ty_gas_concs),       intent(in   ) :: gas_desc   !! Gas volume mixing ratios
    real(wp), dimension(:,:,:), intent(  out) :: layer_mass !! (ngas, ncol, nlay)
    character(len=128),       intent(  out) :: error_msg

    integer :: igas, icol, ilay
    real(wp), dimension(size(layer_mass,1), size(layer_mass,2), size(layer_mass,3)) :: vmr

    error_msg = ""

    vmr(:,:,:) = 0._wp

    ! Get vmr if gas is provided in ty_gas_concs
    do igas = 1, ngas
      if (any(lower_case(this%gas_names(igas)) == gas_desc%get_gas_names())) then
        error_msg = gas_desc%get_vmr(this%gas_names(igas), vmr(igas,:,:))
        if (error_msg /= '') return
      endif
    end do

    ! Convert pressures and vmr to layer masses (ngas, ncol, nlay)
    ! mmr = vmr * (Mgas/Mair)
    ! layer_mass = mmr * dp / g
    do ilay = 1, nlay
      do icol = 1, ncol
        do igas = 1, ngas
          layer_mass(igas, icol, ilay) = vmr(igas, icol, ilay) * &
            (this%mol_weights(igas) / this%m_dry) * &
            abs(plev(icol, ilay+1) - plev(icol, ilay)) / grav
        end do
      end do
    end do

  end subroutine compute_layer_mass
end module mo_optics_ssm
