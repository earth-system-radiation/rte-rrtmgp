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
!
! The gas optics class used by RRMTGP needs to be initialized with data stored in a netCDF file.
!    RRTMGP itself doesn't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward implementation of reading the data
!    and calling gas_optics%load().
!
! -------------------------------------------------------------------------------------------------
module mo_load_coefficients
  !
  ! Modules for working with rte and rrtmgp
  !
  use mo_rte_kind,           only: wp, wl
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp

  use fms_mod, only: error_mesg, fatal
  use fms2_io_mod, only: close_file, FmsNetcdfFile_t, get_dimension_size, open_file, &
                         read_data, variable_exists
  implicit none
  private
  public :: load_and_init

contains


  subroutine stop_on_err(msg)

    character(len=*), intent(in) :: msg

    if (trim(msg) .ne. "") then
      call error_mesg("mo_load_coefficients", trim(msg), fatal)
    endif
  end subroutine


  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_and_init(kdist, filename, available_gases)
    class(ty_gas_optics_rrtmgp), intent(inout) :: kdist
    character(len=*),     intent(in   ) :: filename
    class(ty_gas_concs),  intent(in   ) :: available_gases ! Which gases does the host model have available?
    ! --------------------------------------------------
    !
    ! Variables that will be passed to gas_optics%load()
    !
    type(FmsNetcdfFile_t) :: dataset
    character(len=32), dimension(:), allocatable :: gas_names
    integer,  dimension(:,:,:),      allocatable :: key_species
    integer,  dimension(:,:  ),      allocatable :: band2gpt
    real(wp), dimension(:,:  ),      allocatable :: band_lims
    real(wp)                                     :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:      ),    allocatable :: press_ref
    real(wp), dimension(:      ),    allocatable :: temp_ref
    real(wp), dimension(:,:,:  ),    allocatable :: vmr_ref
    real(wp), dimension(:,:,:,:),    allocatable :: kmajor

    character(len=32), dimension(:),  allocatable :: gas_minor, identifier_minor
    character(len=32), dimension(:),  allocatable :: minor_gases_lower,               minor_gases_upper
    integer, dimension(:,:),          allocatable :: minor_limits_gpt_lower,          minor_limits_gpt_upper
    logical(wl), dimension(:),        allocatable :: minor_scales_with_density_lower, minor_scales_with_density_upper
    character(len=32), dimension(:),  allocatable :: scaling_gas_lower,               scaling_gas_upper
    logical(wl), dimension(:),        allocatable :: scale_by_complement_lower,       scale_by_complement_upper
    integer, dimension(:),            allocatable :: kminor_start_lower,              kminor_start_upper
    real(wp), dimension(:,:,:),       allocatable :: kminor_lower,                    kminor_upper

    real(wp), dimension(:,:,:  ), allocatable :: rayl_lower, rayl_upper
    real(wp), dimension(:      ), allocatable :: solar_quiet, solar_facular, solar_sunspot
    real(wp)                                  :: tsi_default, mg_default, sb_default
    real(wp), dimension(:,:    ), allocatable :: totplnk
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac
    real(wp), dimension(:,:)    , allocatable :: optimal_angle_fit

    ! -----------------
    !
    ! Book-keeping variables
    !
    integer, dimension(:), allocatable :: buffer
    integer :: ntemps,          &
               npress,          &
               nabsorbers,      &
               nextabsorbers,   &
               nminorabsorbers, &
               nmixingfracs,    &
               nlayers,         &
               nbnds,           &
               ngpts,           &
               npairs,          &
               nminor_absorber_intervals_lower, &
               nminor_absorber_intervals_upper, &
               ncontributors_lower, &
               ncontributors_upper, &
               ninternalSourcetemps, &
               nfit_coeffs
    ! --------------------------------------------------
    !
    ! How big are the various arrays?
    !
    if (.not. open_file(dataset, trim(fileName), "read")) then
      call stop_on_err("load_and_init(): can't open file "//trim(fileName)//".")
    endif
    call get_dimension_size(dataset, "temperature", ntemps)
    call get_dimension_size(dataset, "pressure", npress)
    call get_dimension_size(dataset, "absorber", nabsorbers)
    call get_dimension_size(dataset, "minor_absorber", nminorabsorbers)
    call get_dimension_size(dataset, "absorber_ext", nextabsorbers)
    call get_dimension_size(dataset, "mixing_fraction", nmixingfracs)
    call get_dimension_size(dataset, "atmos_layer", nlayers)
    call get_dimension_size(dataset, "bnd", nbnds)
    call get_dimension_size(dataset, "gpt", ngpts)
    call get_dimension_size(dataset, "pair", npairs)
    call get_dimension_size(dataset, "minor_absorber_intervals_lower", &
                            nminor_absorber_intervals_lower)
    call get_dimension_size(dataset, "minor_absorber_intervals_upper", &
                            nminor_absorber_intervals_upper)
    call get_dimension_size(dataset, "temperature_Planck", ninternalSourcetemps)
    call get_dimension_size(dataset, "contributors_lower", ncontributors_lower)
    call get_dimension_size(dataset, "contributors_upper", ncontributors_upper)
    call get_dimension_size(dataset, "fit_coeffs", nfit_coeffs) ! Will be 0 for SW

    ! -----------------
    !
    ! Read the many arrays
    !
    allocate(gas_names(nabsorbers))
    call read_data(dataset, "gas_names", gas_names)
    allocate(key_species(2, nlayers, nbnds))
    call read_data(dataset, "key_species", key_species)
    allocate(band_lims(2, nbnds))
    call read_data(dataset, "bnd_limits_wavenumber", band_lims)
    allocate(band2gpt(2, nbnds))
    call read_data(dataset, "bnd_limits_gpt", band2gpt)
    allocate(press_ref(npress))
    call read_data(dataset, "press_ref", press_ref)
    allocate(temp_ref(ntemps))
    call read_data(dataset, "temp_ref", temp_ref)
    call read_data(dataset, "absorption_coefficient_ref_P", temp_ref_p)
    call read_data(dataset, "absorption_coefficient_ref_T", temp_ref_t)
    call read_data(dataset, "press_ref_trop", press_ref_trop)
    allocate(kminor_lower(ncontributors_lower, nmixingfracs, ntemps))
    call read_data(dataset, "kminor_lower", kminor_lower)
    allocate(kminor_upper(ncontributors_upper, nmixingfracs, ntemps))
    call read_data(dataset, "kminor_upper", kminor_upper)
    allocate(gas_minor(nminorabsorbers))
    call read_data(dataset, "gas_minor", gas_minor)
    allocate(identifier_minor(nminorabsorbers))
    call read_data(dataset, "identifier_minor", identifier_minor)
    allocate(minor_gases_lower(nminor_absorber_intervals_lower))
    call read_data(dataset, "minor_gases_lower", minor_gases_lower)
    allocate(minor_gases_upper(nminor_absorber_intervals_upper))
    call read_data(dataset, "minor_gases_upper", minor_gases_upper)
    allocate(minor_limits_gpt_lower(npairs, nminor_absorber_intervals_lower))
    call read_data(dataset, "minor_limits_gpt_lower", minor_limits_gpt_lower)
    allocate(minor_limits_gpt_upper(npairs, nminor_absorber_intervals_upper))
    call read_data(dataset, "minor_limits_gpt_upper", minor_limits_gpt_upper)
    allocate(buffer(nminor_absorber_intervals_lower))
    allocate(minor_scales_with_density_lower(nminor_absorber_intervals_lower))
    call read_data(dataset, "minor_scales_with_density_lower", buffer)
    minor_scales_with_density_lower(:) = buffer(:) .ne. 0
    allocate(scale_by_complement_lower(nminor_absorber_intervals_lower))
    call read_data(dataset, "scale_by_complement_lower", buffer)
    scale_by_complement_lower(:) = buffer(:) .ne. 0
    deallocate(buffer)
    allocate(buffer(nminor_absorber_intervals_upper))
    allocate(minor_scales_with_density_upper(nminor_absorber_intervals_upper))
    call read_data(dataset, "minor_scales_with_density_upper", buffer)
    minor_scales_with_density_upper(:) = buffer(:) .ne. 0
    allocate(scale_by_complement_upper(nminor_absorber_intervals_upper))
    call read_data(dataset, "scale_by_complement_upper", buffer)
    scale_by_complement_upper(:) = buffer(:) .ne. 0
    deallocate(buffer)
    allocate(scaling_gas_lower(nminor_absorber_intervals_lower))
    call read_data(dataset, "scaling_gas_lower", scaling_gas_lower)
    allocate(scaling_gas_upper(nminor_absorber_intervals_upper))
    call read_data(dataset, "scaling_gas_upper", scaling_gas_upper)
    allocate(kminor_start_lower(nminor_absorber_intervals_lower))
    call read_data(dataset, "kminor_start_lower", kminor_start_lower)
    allocate(kminor_start_upper(nminor_absorber_intervals_upper))
    call read_data(dataset, "kminor_start_upper", kminor_start_upper)
    allocate(vmr_ref(nlayers, nextabsorbers, ntemps))
    call read_data(dataset, "vmr_ref", vmr_ref)
    allocate(kmajor(ngpts, nmixingfracs, npress+1, ntemps))
    call read_data(dataset, "kmajor", kmajor)
    if (variable_exists(dataset, "rayl_lower")) then
      allocate(rayl_lower(ngpts, nmixingfracs, ntemps))
      call read_data(dataset, "rayl_lower", rayl_lower)
      allocate(rayl_upper(ngpts, nmixingfracs, ntemps))
      call read_data(dataset, "rayl_upper", rayl_upper)
    end if
    ! --------------------------------------------------
    !
    ! Initialize the gas optics class with data. The calls look slightly different depending
    !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
    ! gas_optics%load() returns a string; a non-empty string indicates an error.
    !
    if (variable_exists(dataset, "totplnk")) then
      !
      ! If there's a totplnk variable in the file it's a longwave (internal sources) type
      !
      allocate(totplnk(ninternalsourcetemps, nbnds))
      call read_data(dataset, "totplnk", totplnk)
      allocate(planck_frac(ngpts, nmixingfracs, npress+1, ntemps))
      call read_data(dataset, "plank_fraction", planck_frac)
      allocate(optimal_angle_fit(nfit_coeffs, nbnds))
      call read_data(dataset, "optimal_angle_fit", optimal_angle_fit)
      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor, &
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  totplnk, planck_frac,       &
                                  rayl_lower, rayl_upper, &
                                  optimal_angle_fit))
    else
      !
      ! Solar source doesn't have an dependencies yet
      !
      allocate(solar_quiet(ngpts))
      call read_data(dataset, "solar_source_quiet", solar_quiet)
      allocate(solar_facular(ngpts))
      call read_data(dataset, "solar_source_facular", solar_facular)
      allocate(solar_sunspot(ngpts))
      call read_field(dataset, "solar_source_sunspot", solar_sunspot)
      call read_field(dataset, "tsi_default", tsi_default)
      call read_field(dataset, "mg_default", mg_default)
      call read_field(dataset, "sb_default", sb_default)
      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor,&
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  solar_quiet, solar_facular, solar_sunspot, &
                                  tsi_default, mg_default, sb_default, &
                                  rayl_lower, rayl_upper))
    end if
    ! --------------------------------------------------
    call close_file(dataset)
    if (allocated(gas_names)) deallocate(gas_names)
    if (allocated(key_species)) deallocate(key_species)
    if (allocated(band2gpt)) deallocate(band2gpt)
    if (allocated(band_lims)) deallocate(band_lims)
    if (allocated(press_ref)) deallocate(press_ref)
    if (allocated(temp_ref)) deallocate(temp_ref)
    if (allocated(vmr_ref)) deallocate(vmr_ref)
    if (allocated(kmajor)) deallocate(kmajor)
    if (allocated(gas_minor)) deallocate(gas_minor)
    if (allocated(identifier_minor)) deallocate(identifier_minor)
    if (allocated(minor_gases_lower)) deallocate(minor_gases_lower)
    if (allocated(minor_gases_upper)) deallocate(minor_gases_upper)
    if (allocated(minor_limits_gpt_lower)) deallocate(minor_limits_gpt_lower)
    if (allocated(minor_limits_gpt_upper)) deallocate(minor_limits_gpt_upper)
    if (allocated(minor_scales_with_density_lower)) deallocate(minor_scales_with_density_lower)
    if (allocated(minor_scales_with_density_upper)) deallocate(minor_scales_with_density_upper)
    if (allocated(scaling_gas_lower)) deallocate(scaling_gas_lower)
    if (allocated(scaling_gas_upper)) deallocate(scaling_gas_upper)
    if (allocated(scale_by_complement_lower)) deallocate(scale_by_complement_lower)
    if (allocated(scale_by_complement_upper)) deallocate(scale_by_complement_upper)
    if (allocated(kminor_start_lower)) deallocate(kminor_start_lower)
    if (allocated(kminor_start_upper)) deallocate(kminor_start_upper)
    if (allocated(kminor_lower)) deallocate(kminor_lower)
    if (allocated(kminor_upper)) deallocate(kminor_upper)
    if (allocated(rayl_lower)) deallocate(rayl_lower)
    if (allocated(rayl_upper)) deallocate(rayl_upper)
    if (allocated(solar_quiet)) deallocate(solar_quiet)
    if (allocated(solar_facular)) deallocate(solar_facular)
    if (allocated(solar_sunspot)) deallocate(solar_sunspot)
    if (allocated(totplnk)) deallocate(totplnk)
    if (allocated(planck_frac)) deallocate(planck_frac)
    if (allocated(optimal_angle_fit)) deallocate(optimal_angle_fit)
  end subroutine load_and_init
end module
