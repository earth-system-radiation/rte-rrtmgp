! See documentation in other modules...
!
! Contacts: Andrew Williams and Robert Pincus
! email:  andrewwilliams@ucsd.edu
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
  use mo_rte_util_array,     only: zero_array, set_to_scalar
  use mo_rte_util_array_validation, &
                             only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props,      only: ty_optical_props
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_optics_util_string, only: lower_case
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_constants,   only: grav
  use mo_optics_ssm_kernels, only: compute_tau, compute_Planck_source, compute_layer_mass

  implicit none
  interface configure
    module procedure configure_with_values, configure_with_defaults
  end interface

  private
  public :: configure
  real(wp), parameter, public :: Tsun_ssm = 5760._wp ! Default sun temperature for SSM
  real(wp), parameter, public :: tsi = 1360._wp ! Default total solar irradiance

  real(wp), parameter, public :: mw_h2o = 0.018_wp ! Molecular weight h2o
  real(wp), parameter, public :: mw_co2 = 0.044_wp ! Molecular weight co2
  real(wp), parameter, public :: mw_o3  = 0.048_wp ! Molecular weight o3

  real(wp), parameter, public :: kappa_cld_lw = 50._wp ! Default for lw cloud absorption coefficient (m2/kg)
  real(wp), parameter, public :: kappa_cld_sw = 0.0001_wp ! Default for sw cloud absorption coefficient (m2/kg)

  real(wp), parameter, public :: ssa_cld_lw = 0._wp ! Default for lw cloud single scattering albedo
  real(wp), parameter, public :: ssa_cld_sw = 0.9999_wp ! Default for sw cloud single scattering albedo

  real(wp), parameter, public :: g_cld_lw = 0._wp ! Default for lw cloud asymmetry
  real(wp), parameter, public :: g_cld_sw = 0.85_wp ! Default for sw cloud asymmetry

  ! Default wavenumber arrays
  integer :: i
  integer,  parameter :: nnu_def = 41           ! Default nnu
  real(wp), parameter :: nu_min_lw_def = 0._wp    ! cm⁻¹
  real(wp), parameter :: nu_max_lw_def = 3500._wp  ! cm⁻¹

  real(wp), parameter :: nu_min_sw_def = 0._wp    ! cm⁻¹
  real(wp), parameter :: nu_max_sw_def = 50000._wp  ! cm⁻¹

  real(wp), dimension(nnu_def), parameter :: nus_lw_def = &
    [(50._wp + (i-1) * (3000._wp - 50._wp) / (nnu_def - 1), i = 1, nnu_def)] ! Default wavenumber array, LW

  real(wp), dimension(nnu_def), parameter :: nus_sw_def = &
    [(1000._wp + (i-1) * (45000._wp - 1000._wp) / (nnu_def - 1), i = 1, nnu_def)] ! Default wavenumber array, SW

  ! Default spectoscopic params
  real(wp), dimension(4,4), parameter :: triangle_params_def_lw = reshape( &
    [1._wp, 298._wp,    0._wp, 64._wp,  &
     1._wp,  12._wp, 1600._wp, 36._wp,  &
     1._wp,  16._wp, 1600._wp, 54._wp,  &
     2._wp, 110._wp,  667._wp, 12._wp], &
    shape = [4, 4], order = [2, 1])

  character(len=32), dimension(2), parameter :: gas_names_def_lw = [character(32) :: "h2o", "co2"]

  real(wp), dimension(2,4), parameter :: triangle_params_def_sw = reshape( &
    [1._wp,   1._wp,  0._wp,    1200._wp,  &
     2._wp,   0._wp,  0._wp, 1000000._wp], &  ! No O3 triangle for now, todo
    shape = [2, 4], order = [2, 1])

  character(len=32), dimension(2), parameter :: gas_names_def_sw = [character(32) :: "h2o", "o3"]

  !
  !
  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_gas_optics), public :: ty_optics_ssm
    private
    character(32), dimension(:),    allocatable :: gas_names
    real(wp),      dimension(:),    allocatable :: mol_weights
    real(wp),      dimension(:, :), allocatable :: absorption_coeffs
      ! total absorption coefficients at spectral points (nnus, ngases)
    real(wp),      dimension(:),    allocatable :: nus, dnus
      ! spectral discretization
    real(wp),      dimension(:),    allocatable :: toa_src
      ! determined at configuration time
    real(wp) :: Tstar     = 0._wp, &
                pref      = 500._wp * 100._wp, & ! reference pressure, assumes 500 hPa by default
                m_dry     = 0.029_wp, & ! molecular weight of dry air [kg/mol]
                tsi       = 0._wp, &           ! total solar irradiance
                kappa_cld = 0._wp, &              ! cloud absorption coefficient [m2/kg]
                g_cld     = 0._wp, &              ! cloud asymmetry factor
                ssa_cld   = 0._wp                 ! cloud single scattering albedo
    contains
      procedure, private :: configure_with_values
      procedure, private :: configure_with_defaults
      generic,   public  :: configure => configure_with_values, configure_with_defaults
      procedure, public  :: gas_optics_int
      procedure, public  :: gas_optics_ext
      procedure, public  :: cloud_optics
      procedure, public  :: source_is_internal
      procedure, public  :: source_is_external
      procedure, public  :: get_press_min
      procedure, public  :: get_press_max
      procedure, public  :: get_temp_min
      procedure, public  :: get_temp_max

  end type ty_optics_ssm

contains
  !--------------------------------------------------------------------------------------------------------------------
  function configure_with_defaults(this, do_sw) result(error_msg)
    class(ty_optics_ssm), intent(inout) :: this
    logical, optional,    intent(in)    :: do_sw ! Configure with SW defaults? Else LW
    character(len=128)                  :: error_msg

    logical :: do_sw_local

    do_sw_local = .false.
    if (present(do_sw)) do_sw_local = do_sw

    if (.not. do_sw_local) then
      error_msg = this%configure_with_values(gas_names_def_lw, triangle_params_def_lw, &
                                             nus_lw_def, nu_min_lw_def, nu_max_lw_def, &
                                             kappa_cld=kappa_cld_lw, g_cld=g_cld_lw, ssa_cld=ssa_cld_lw)
    else
      error_msg = this%configure_with_values(gas_names_def_sw, triangle_params_def_sw, &
                                             nus_sw_def, nu_min_sw_def, nu_max_sw_def, &
                                             Tstar=Tsun_ssm, tsi=tsi,                  &
                                             kappa_cld=kappa_cld_sw, g_cld=g_cld_sw, ssa_cld=ssa_cld_sw)
    end if
  end function configure_with_defaults

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Configure the simple spectral model parameters
  !> Gas optics for simple spectral models are described as N triangles of ln(kappa) for
  !>   each of M gases.
  !> Each triangle id defined by
  !>    the absorption coefficient at central wavenumber, nu0
  !>    a central wavenumber nu0
  !>    the slope d ln(kappa)/d nu
  !> By default there are two gases (h2o and co2) (for the LW).
  !>   h2o has three triangles (rotational lines, and a continuum with infinite slope)
  !>   co2 has a single triangle
  !> Absorption coefficents are defined at pref and pressure-broadened as p/pref
  !> Spectral discretization is defined as
  !>   max and min wavenumbers and
  !>   a set of wavenumbers at which to evaluate absorption and Planck function
  !> If Tstar is provided the gas optics are for insolation
  !
  function configure_with_values(this,           &
                     gas_names, triangle_params, &
                     nus, nu_min, nu_max,        &
                     Tstar, tsi,                 &
                     kappa_cld, g_cld, ssa_cld) result(error_msg)
    !
    !
    class(ty_optics_ssm),          intent(inout) :: this
    character(32), dimension(:),   intent(in   ) :: gas_names
    real(wp),      dimension(:,:), intent(in   ) :: triangle_params
      !! (ntriangles, 4) where the second dimension holds [gas_index, kappa0, nu0, l]
    real(wp),      dimension(:),   intent(in   ) :: nus
      !! Wavenumbers at which to evaluate Planck function and absorption coefficient
    real(wp),                      intent(in   ) :: nu_min, nu_max
      !! Upper and lower bounds of spectrum
    real(wp),      optional,       intent(in   ) :: Tstar
      !! Temperature for stellar insolation
    real(wp),      optional,       intent(in   ) :: tsi
      !! Total solar irradiance
    real(wp),      optional,       intent(in   ) :: kappa_cld
    real(wp),      optional,       intent(in   ) :: g_cld
    real(wp),      optional,       intent(in   ) :: ssa_cld
      !! cloud optical properties
    character(len=128)                      :: error_msg     !! Empty if successful
    ! ----------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------
    integer :: nnu, ngas
    integer :: inu, itri, igas
    real(wp), dimension(2, size(nus)) :: band_lims_wavenum
    real(wp), dimension(1, size(nus)) :: toa_src_1col

    error_msg = ""
    ngas = size(gas_names)
    nnu  = size(nus)

    !
    ! Input sanitizing
    ! triangle params: index <= ngases, kappa0 >= 0; nu_min < nu0s < nu_max; l > 0
    ! nu_min <= nus <= max_nu
    ! Tstar > 0 if specified
    !

    if (.not. all(nus > nu_min .and. nus < nu_max)) then
      error_msg = "ssm_gas_optics(): nu must be less than nu_max and greater than nu_min"
    end if

    if (present(Tstar)) then
      if (Tstar <= 0.0_wp) then
        error_msg = "ssm_gas_optics(): if specified Tstar must be > 0"
      end if
    end if

    if (present(tsi)) then
      if (tsi <= 0.0_wp) then
        error_msg = "ssm_gas_optics(): if specified tsi must be > 0"
      end if
    end if

    ! Validate gas indices in triangle_params
    if (.not. all(triangle_params(:, 1) >= 1._wp .and. &
                  triangle_params(:, 1) <= real(ngas, wp) .and. &
                  triangle_params(:, 1) == floor(triangle_params(:, 1)))) then
      error_msg = "ssm_gas_optics(): gas index in triangle_params must be integer >= 1 and <= ngas"
      return
    end if
    if (.not. all(triangle_params(:, 2) >= 0.0_wp)) then
      error_msg = "ssm_gas_optics(): kappa0 needs to be >=0"
    end if

    if (.not. all(triangle_params(:, 4) > 0.0_wp)) then
      error_msg = "ssm_gas_optics(): l needs to be > 0"
    end if
    if (error_msg /= '') return

    if(allocated(this%gas_names)) &
      deallocate(this%gas_names, this%mol_weights, this%nus, this%dnus, &
                 this%absorption_coeffs, this%toa_src)

    allocate(this%gas_names(ngas), &
             this%mol_weights(ngas), &
             this%nus (nnu), &
             this%dnus(nnu), &
             this%absorption_coeffs(ngas,nnu), &
             this%toa_src(nnu))

    this%nus(1:nnu)        = nus(1:nnu)
    this%gas_names(1:ngas) = gas_names(1:ngas)

    !
    ! Spectral discretization - edges of "bands"
    ! Then construct the band limits (place edges at midpoints between nus values):
    !
    ! First band: starts at nu_min
    band_lims_wavenum(1, 1) = nu_min
    band_lims_wavenum(2, 1) = (nus(1) + nus(2)) * 0.5_wp

    ! Middle bands: edges at midpoints
    do inu = 2, nnu - 1
      band_lims_wavenum(1, inu) = (nus(inu-1) + nus(inu))   * 0.5_wp
      band_lims_wavenum(2, inu) = (nus(inu)   + nus(inu+1)) * 0.5_wp
    end do

    ! Last band: ends at nu_max
    band_lims_wavenum(1, nnu) = (nus(nnu-1) + nus(nnu)) * 0.5_wp
    band_lims_wavenum(2, nnu) = nu_max

    !
    ! Initialize the parent class
    error_msg = this%ty_optical_props%init(band_lims_wavenum, name="ssm")
    if (error_msg /= '') return

    ! Now you can get dnus from the initialized structure:
    this%dnus = band_lims_wavenum(2, :) - band_lims_wavenum(1, :)

    ! Set molar masses based on gas names
    ! Maybe this is better as module data...
    do igas = 1, ngas
      select case (trim(lower_case(gas_names(igas))))
        case ('h2o')
          this%mol_weights(igas) = mw_h2o
        case ('co2')
          this%mol_weights(igas) = mw_co2
        case ('o3')
          this%mol_weights(igas) = mw_o3
        case default
          error_msg = "Don't know the molecular weight for gas: " // trim(gas_names(igas))
          return
      end select
    end do
    if (error_msg /= '') return

    ! Compute absorption coefficients by summing exponentials at each nu
    ! Initialize absorption coefficients to zero
    this%absorption_coeffs(:,:) = 0._wp

    do itri = 1, size(triangle_params, 1)  ! Loop over triangles
      igas = int(triangle_params(itri, 1)) ! Identify which gas this triangle corresponds to
      do inu = 1, nnu                      ! Loop over wavenumber
        this%absorption_coeffs(igas, inu) = &
          this%absorption_coeffs(igas, inu) + &
          triangle_params(itri, 2) * exp(-abs(this%nus(inu) - triangle_params(itri, 3)) / triangle_params(itri, 4))
      end do
    end do

    if(present(Tstar)) this%Tstar = Tstar
    if(present(tsi))   this%tsi = tsi

    if(this%Tstar > 0) then
      !
      ! Shortwave: incoming solar irradiance
      !   (Need to have a 1xnnu array for compute_Planck_source)
      !
      call compute_Planck_source(1, nnu, &
                                 this%nus, this%dnus, [this%Tstar],   &
                                 toa_src_1col)

      ! Tstar sets spectral distribution; normalize to top-of-atmosphere
      !   total solar irradiance (flux)
      this%toa_src(:) = toa_src_1col(1,:) * this%tsi/sum(toa_src_1col(1,:))
    else
      call zero_array(nnu, this%toa_src)
    end if

    if(present(kappa_cld)) then
      this%kappa_cld = kappa_cld
    else
      this%kappa_cld = 0._wp
    end if

    if(present(g_cld)) then
      this%g_cld = g_cld
    else
      this%g_cld = 0._wp
    end if

    if(present(ssa_cld)) then
      this%ssa_cld = ssa_cld
    else
      this%ssa_cld = 0._wp
    end if

    !$acc        enter data copyin(this%gas_names, this%mol_weights, this%absorption_coeffs, &
    !$acc                          this%nus, this%dnus, this%toa_src)
    !$omp target enter data map(to:this%gas_names, this%mol_weights, this%absorption_coeffs, &
    !$omp                          this%nus, this%dnus, this%toa_src)

  end function configure_with_values

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
    !
    integer :: ncol, nlay, nnu, ngas
    integer :: igas, icol, ilay, idx_gas
    real(wp), dimension(size(this%gas_names), size(play,1), size(play,2)) :: layer_mass
    ! ----------------------------------------------------------
    error_msg = ""

    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt()
    ngas = size(this%gas_names)

    ! Ideally some data sanitization goes here - size of optical props agrees with size of play etc.
    !   use extents_are() function

    call get_layer_mass(ncol, nlay, ngas, &
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
    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    ! Planck function sources
    !
    call compute_Planck_source(ncol, nlay,   nnu, &
                               this%nus, this%dnus, tlay,   &
                               sources%lay_source)
    ! This will fail if Tlev isn't provided
    !   There's interpolation code in RRTMGP gas optics -
    !   should we make this generic and package it with the gas optics type?
    call compute_Planck_source(ncol, nlay+1, nnu, &
                               this%nus, this%dnus, tlev,   &
                               sources%lev_source)
    call compute_Planck_source(ncol, nnu, &
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
    real(wp), dimension(:,:), intent(  out) :: toa_src     !! Incoming solar irradiance(ncol, nnu)
    character(len=128)                      :: error_msg   !! Empty if successful

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    !
    integer :: ncol, nlay, nnu, ngas
    integer :: igas, icol, ilay, inu, idx_gas
    real(wp), dimension(size(this%gas_names), size(play,1), size(play,2)) :: layer_mass
    error_msg = ""

    ! ----------------------------------------------------------
    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt()
    ngas = size(this%gas_names)

    ! Ideally some data sanitization goes here - size of optical props agrees with size of play etc.
    !   use extents_are() function

    call get_layer_mass(ncol, nlay, ngas, &
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
    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    ! Not doing scattering of gases
    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    !$acc parallel loop collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do inu = 1, this%get_ngpt()
      do icol = 1, optical_props%get_ncol()
        toa_src(icol, inu) = this%toa_src(inu)
      end do
    end do
  end function gas_optics_ext

  !------------------------------------------------------------------------------------------
  !
  !> Derive cloud optical properties from provided cloud physical properties
  !
  function cloud_optics(this,                     &
                        clwp, ciwp, reliq, reice, &
                        optical_props) result(error_msg)
    class(ty_optics_ssm), &
              intent(in   ) :: this
    real(wp), intent(in   ) :: clwp  (:,:), &   ! cloud liquid water path (g/m2)
                               ciwp  (:,:), &   ! cloud ice water path    (g/m2)
                               reliq (:,:), &   ! cloud liquid particle effective size (microns)
                               reice (:,:)      ! cloud ice particle effective radius  (microns)
    class(ty_optical_props_arry), &
              intent(inout) :: optical_props
    character(len=128)      :: error_msg
    ! ----------------------------------------------------------
    ! Local variables
    integer :: icol, ilay, inu
    ! ----------------------------------------------------------
    error_msg = ""

    ! Get cloud optical depth by multiplying
    ! [kg/m2] of cloud by [m2/kg] absorption coeff

    !$acc parallel loop collapse(3) copyout(optical_props%tau)
    !$omp target teams distribute parallel do simd collapse(3) map(from:optical_props%tau)
    do inu = 1, this%get_ngpt()
      do ilay = 1, optical_props%get_nlay()
        do icol = 1, optical_props%get_ncol()
          optical_props%tau(icol, ilay, inu) = 1000._wp * (clwp(icol, ilay) + ciwp(icol, ilay)) * this%kappa_cld
        end do
      end do
    end do

    select type(optical_props)
      type is (ty_optical_props_2str)
          call set_to_scalar(optical_props%get_ncol(), &
                             optical_props%get_nlay(), &
                             optical_props%get_ngpt(), &
                         optical_props%ssa, this%ssa_cld)
          call set_to_scalar(optical_props%get_ncol(), &
                             optical_props%get_nlay(), &
                             optical_props%get_ngpt(), &
                         optical_props%g,    this%g_cld)
      type is (ty_optical_props_nstr)
        ! Toss an error here? No one uses n-streams
    end select

  end function cloud_optics

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Compute layer masses from gas concentrations and pressure levels
  !
  subroutine get_layer_mass(ncol, nlay, ngas, this, plev, gas_desc, layer_mass, error_msg)
    integer,  intent(in ) :: ncol, nlay, ngas
    class(ty_optics_ssm),     intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: plev       !! level pressures [Pa]; (ncol,nlay+1)
    type(ty_gas_concs),       intent(in   ) :: gas_desc   !! Gas volume mixing ratios
    real(wp), dimension(:,:,:), intent(  out) :: layer_mass !! (ngas, ncol, nlay)
    character(len=128),       intent(  out) :: error_msg

    integer :: igas
    real(wp), dimension(size(layer_mass,1), size(layer_mass,2), size(layer_mass,3)) :: vmr
    ! --------------
    error_msg = ""

    !$acc        data create(   vmr)
    !$omp target data map(alloc:vmr)
    call zero_array(size(vmr,1), size(vmr,2), size(vmr,3), vmr)

    ! Get vmr if gas is provided in ty_gas_concs
    do igas = 1, ngas
      if (any(lower_case(this%gas_names(igas)) == gas_desc%get_gas_names())) then
        error_msg = gas_desc%get_vmr(this%gas_names(igas), vmr(igas,:,:))
        if (error_msg /= '') return
      endif
    end do

    ! Convert pressures and vmr to layer masses (ngas, ncol, nlay)
    call compute_layer_mass(ncol, nlay, ngas, vmr, plev, this%mol_weights, this%m_dry, layer_mass)

    !$acc end data
    !$omp end target data

  end subroutine get_layer_mass

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

end module mo_optics_ssm
