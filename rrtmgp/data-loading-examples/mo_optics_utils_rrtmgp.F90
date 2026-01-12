! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2024-  Atmospheric and Environmental Research,
!    Regents of the University of Colorado,
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ----------------------------------------------------------------------------
!!
! Gas, cloud, and aerosol optics classes need to be initialized with data; for RRTMGP data comes from a netCDF file.
!    The gas optics classes themselves don't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward serial implementation of reading the data
!    and calling gas_optics%load().
!
!
module mo_optics_utils_rrtmgp
  use mo_rte_kind,           only: wp, wl
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics_rrtmgp, &
                             only: ty_cloud_optics_rrtmgp
  use mo_aerosol_optics_rrtmgp_merra,  &
                             only: ty_aerosol_optics_rrtmgp_merra
  use mo_testing_utils,      only: stop_on_err
  ! --------------------------------------------------
  use mo_simple_netcdf, only: read_field, read_char_vec, read_logical_vec, dim_exists, var_exists, get_dim_size
  use netcdf
  implicit none

  private
  public :: load_gas_optics, load_cloud_optics, load_aerosol_optics

contains
  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_gas_optics(kdist, filename, available_gases)
    class(ty_gas_optics_rrtmgp), intent(inout) :: kdist
    character(len=*),            intent(in   ) :: filename
    class(ty_gas_concs),         intent(in   ) :: available_gases ! Which gases does the host model have available?
    ! --------------------------------------------------
    !
    ! Variables that will be passed to gas_optics%load()
    !
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
    integer :: ncid
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
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("load_gas_optics(): can't open file " // trim(fileName))
    ntemps            = get_dim_size(ncid,'temperature')
    npress            = get_dim_size(ncid,'pressure')
    nabsorbers        = get_dim_size(ncid,'absorber')
    nminorabsorbers   = get_dim_size(ncid,'minor_absorber')
    nextabsorbers     = get_dim_size(ncid,'absorber_ext')
    nmixingfracs      = get_dim_size(ncid,'mixing_fraction')
    nlayers           = get_dim_size(ncid,'atmos_layer')
    nbnds             = get_dim_size(ncid,'bnd')
    ngpts             = get_dim_size(ncid,'gpt')
    npairs            = get_dim_size(ncid,'pair')
    nminor_absorber_intervals_lower &
                      = get_dim_size(ncid,'minor_absorber_intervals_lower')
    nminor_absorber_intervals_upper  &
                      = get_dim_size(ncid,'minor_absorber_intervals_upper')
    ninternalSourcetemps &
                      = get_dim_size(ncid,'temperature_Planck')
    ncontributors_lower = get_dim_size(ncid,'contributors_lower')
    ncontributors_upper = get_dim_size(ncid,'contributors_upper')
    nfit_coeffs         = get_dim_size(ncid,'fit_coeffs') ! Will be 0 for SW

    ! -----------------
    !
    ! Read the many arrays
    !
    gas_names         = read_char_vec(ncid, 'gas_names', nabsorbers)
    key_species       = int(read_field(ncid, 'key_species',  2, nlayers, nbnds))
    band_lims         = read_field(ncid, 'bnd_limits_wavenumber', 2, nbnds)
    band2gpt          = int(read_field(ncid, 'bnd_limits_gpt', 2, nbnds))
    press_ref         = read_field(ncid, 'press_ref', npress)
    temp_ref          = read_field(ncid, 'temp_ref',  ntemps)
    temp_ref_p        = read_field(ncid, 'absorption_coefficient_ref_P')
    temp_ref_t        = read_field(ncid, 'absorption_coefficient_ref_T')
    press_ref_trop    = read_field(ncid, 'press_ref_trop')
    kminor_lower      = read_field(ncid, 'kminor_lower', &
        ncontributors_lower, nmixingfracs, ntemps)
    kminor_upper      = read_field(ncid, 'kminor_upper', &
        ncontributors_upper, nmixingfracs, ntemps)
    gas_minor = read_char_vec(ncid, 'gas_minor', nminorabsorbers)
    identifier_minor = read_char_vec(ncid, 'identifier_minor', nminorabsorbers)
    minor_gases_lower = read_char_vec(ncid, 'minor_gases_lower', nminor_absorber_intervals_lower)
    minor_gases_upper = read_char_vec(ncid, 'minor_gases_upper', nminor_absorber_intervals_upper)
    minor_limits_gpt_lower &
                      = int(read_field(ncid, 'minor_limits_gpt_lower', npairs,nminor_absorber_intervals_lower))
    minor_limits_gpt_upper &
                      = int(read_field(ncid, 'minor_limits_gpt_upper', npairs,nminor_absorber_intervals_upper))
    minor_scales_with_density_lower &
                      = read_logical_vec(ncid, 'minor_scales_with_density_lower', nminor_absorber_intervals_lower)
    minor_scales_with_density_upper &
                      = read_logical_vec(ncid, 'minor_scales_with_density_upper', nminor_absorber_intervals_upper)
    scale_by_complement_lower &
                      = read_logical_vec(ncid, 'scale_by_complement_lower', nminor_absorber_intervals_lower)
    scale_by_complement_upper &
                      = read_logical_vec(ncid, 'scale_by_complement_upper', nminor_absorber_intervals_upper)
    scaling_gas_lower &
                      = read_char_vec(ncid, 'scaling_gas_lower', nminor_absorber_intervals_lower)
    scaling_gas_upper &
                      = read_char_vec(ncid, 'scaling_gas_upper', nminor_absorber_intervals_upper)
    kminor_start_lower &
                      = int(read_field(ncid, 'kminor_start_lower', nminor_absorber_intervals_lower))
    kminor_start_upper &
                      = int(read_field(ncid, 'kminor_start_upper', nminor_absorber_intervals_upper))
    vmr_ref           = read_field(ncid, 'vmr_ref', nlayers, nextabsorbers, ntemps)

    kmajor            = read_field(ncid, 'kmajor',  ngpts, nmixingfracs,  npress+1, ntemps)
    if(var_exists(ncid, 'rayl_lower')) then
      rayl_lower = read_field(ncid, 'rayl_lower',   ngpts, nmixingfracs,            ntemps)
      rayl_upper = read_field(ncid, 'rayl_upper',   ngpts, nmixingfracs,            ntemps)
    end if
    ! --------------------------------------------------
    !
    ! Initialize the gas optics class with data. The calls look slightly different depending
    !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
    ! gas_optics%load() returns a string; a non-empty string indicates an error.
    !
    if(var_exists(ncid, 'totplnk')) then
      !
      ! If there's a totplnk variable in the file it's a longwave (internal sources) type
      !
      totplnk     = read_field(ncid, 'totplnk', ninternalSourcetemps, nbnds)
      planck_frac = read_field(ncid, 'plank_fraction', ngpts, nmixingfracs, npress+1, ntemps)
      optimal_angle_fit = read_field(ncid, 'optimal_angle_fit', nfit_coeffs, nbnds)
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
      solar_quiet   = read_field(ncid, 'solar_source_quiet', ngpts)
      solar_facular = read_field(ncid, 'solar_source_facular', ngpts)
      solar_sunspot = read_field(ncid, 'solar_source_sunspot', ngpts)
      tsi_default   = read_field(ncid, 'tsi_default')
      mg_default    = read_field(ncid, 'mg_default')
      sb_default    = read_field(ncid, 'sb_default')
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
    ncid = nf90_close(ncid)
  end subroutine load_gas_optics
  !--------------------------------------------------------------------------------------------------------------------
!
  ! read cloud optical property LUT coefficients from NetCDF file
  !
  subroutine load_cloud_optics(cloud_spec, cld_coeff_file)
    class(ty_cloud_optics_rrtmgp),         intent(inout) :: cloud_spec
    character(len=*),                      intent(in   ) :: cld_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrghice, nsize_liq, nsize_ice, nspec

    ! Lookup table interpolation constants
    real(wp) :: radliq_lwr          ! liquid particle size lower bound for interpolation
    real(wp) :: radliq_upr          ! liquid particle size upper bound for interpolation
    real(wp) :: diamice_lwr         ! ice particle size lower bound for interpolation
    real(wp) :: diamice_upr         ! ice particle size upper bound for interpolation
    ! LUT coefficients
    real(wp), dimension(:,:),   allocatable :: extliq   ! extinction: liquid
    real(wp), dimension(:,:),   allocatable :: ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:),   allocatable :: asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:), allocatable :: extice   ! extinction: ice
    real(wp), dimension(:,:,:), allocatable :: ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:), allocatable :: asyice   ! asymmetry parameter: ice

    logical(wl) :: defined_on_gpts
    integer,  dimension(:,:), allocatable :: band_lims_gpt
    real(wp), dimension(:,:), allocatable :: band_lims_wvn
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cloud_optics(): can't open file " // trim(cld_coeff_file))

    defined_on_gpts = dim_exists(ncid, 'ngpt')
    ! Read LUT coefficient dimensions
    nband = get_dim_size(ncid,'nband')
    if (defined_on_gpts) then
       nspec = get_dim_size(ncid,'ngpt')
    else
       nspec = nband
    endif
    nrghice   = get_dim_size(ncid,'nrghice')
    nsize_liq = get_dim_size(ncid,'nsize_liq')
    nsize_ice = get_dim_size(ncid,'nsize_ice')

    ! Read LUT constants
    radliq_lwr  = read_field(ncid, 'radliq_lwr')
    radliq_upr  = read_field(ncid, 'radliq_upr')
    diamice_lwr = read_field(ncid, 'diamice_lwr')
    diamice_upr = read_field(ncid, 'diamice_upr')

    ! Allocate cloud property lookup table input arrays
    allocate(extliq(nsize_liq, nspec), &
             ssaliq(nsize_liq, nspec), &
             asyliq(nsize_liq, nspec), &
             extice(nsize_ice, nspec, nrghice), &
             ssaice(nsize_ice, nspec, nrghice), &
             asyice(nsize_ice, nspec, nrghice))

    ! Read LUT coefficients
     extliq = read_field(ncid, 'extliq',  nsize_liq, nspec)
     ssaliq = read_field(ncid, 'ssaliq',  nsize_liq, nspec)
     asyliq = read_field(ncid, 'asyliq',  nsize_liq, nspec)
     extice = read_field(ncid, 'extice',  nsize_ice, nspec, nrghice)
     ssaice = read_field(ncid, 'ssaice',  nsize_ice, nspec, nrghice)
     asyice = read_field(ncid, 'asyice',  nsize_ice, nspec, nrghice)

    ! Read band wavenumber limits
    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)
    if(defined_on_gpts) then
      allocate(band_lims_gpt(2, nband))
      band_lims_gpt = read_field(ncid, 'bnd_limits_gpt', 2, nband)
      call stop_on_err(cloud_spec%load(band_lims_wvn, &
                                       radliq_lwr, radliq_upr, &
                                       diamice_lwr, diamice_upr, &
                                       extliq, ssaliq, asyliq, &
                                       extice, ssaice, asyice, band_lims_gpt))
    else
      call stop_on_err(cloud_spec%load(band_lims_wvn, &
                                       radliq_lwr, radliq_upr, &
                                       diamice_lwr, diamice_upr, &
                                       extliq, ssaliq, asyliq, &
                                       extice, ssaice, asyice))
    end if

    ncid = nf90_close(ncid)
  end subroutine load_cloud_optics
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read aerosol optical property LUT coefficients from NetCDF file
  !
  subroutine load_aerosol_optics(aerosol_spec, aero_coeff_file)
    class(ty_aerosol_optics_rrtmgp_merra),   &
                                intent(inout) :: aerosol_spec
    character(len=*),           intent(in   ) :: aero_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrh, nbin, nval, npair

    real(wp), dimension(:,:),   allocatable :: band_lims_wvn    ! spectral band wavenumber limits (npair,nband)
    ! LUT interpolation arrays
    real(wp), dimension(:,:),   allocatable :: merra_aero_bin_lims ! aerosol size bin limits (npair, nbin)
    real(wp), dimension(:),     allocatable :: aero_rh             ! aerosol relative humidity values (nrh)
    ! LUT coefficients (extinction: m2/kg; ssa, g: unitless)
    real(wp), dimension(:,:,:),   allocatable :: aero_dust_tbl    ! extinction, ssa, g: dust (nval,nbin,nband)
    real(wp), dimension(:,:,:,:), allocatable :: aero_salt_tbl    ! extinction, ssa, g: sea salt (nval,nrh,nbin,nband)
    real(wp), dimension(:,:,:),   allocatable :: aero_sulf_tbl    ! extinction, ssa, g: sulfate (nval,nrh,nband)
    real(wp), dimension(:,:),     allocatable :: aero_bcar_tbl    ! extinction, ssa, g: black carbon, hydrophobic (nval,nband)
    real(wp), dimension(:,:,:),   allocatable :: aero_bcar_rh_tbl ! extinction, ssa, g: black carbon, hydrophilic (nval,nrh,nband)
    real(wp), dimension(:,:),     allocatable :: aero_ocar_tbl    ! extinction, ssa, g: organic carbon, hydrophobic (nval,nband)
    real(wp), dimension(:,:,:),   allocatable :: aero_ocar_rh_tbl ! extinction, ssa, g: organic carbon, hydrophilic (nval,nrh,nband)
    ! -----------------
    ! Open aerosol optical property coefficient file
    if(nf90_open(trim(aero_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_aersol_optics(): can't open file " // trim(aero_coeff_file))

    ! Read LUT coefficient dimensions
    nband = get_dim_size(ncid,'nband')
    nrh   = get_dim_size(ncid,'nrh')
    nbin  = get_dim_size(ncid,'nbin')
    nval  = get_dim_size(ncid,'nval')
    npair  = get_dim_size(ncid,'pair')

    allocate(band_lims_wvn(npair, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)

    ! Allocate aerosol property lookup table interpolation arrays
    allocate(merra_aero_bin_lims(npair,nbin), &
             aero_rh(nrh))

    ! Read LUT interpolation arryas
    merra_aero_bin_lims = read_field(ncid, 'merra_aero_bin_lims',  npair, nbin)
    aero_rh= read_field(ncid, 'aero_rh',  nrh)

    ! Allocate aerosol property lookup table input arrays
    allocate(aero_dust_tbl(nval, nbin, nband), &
             aero_salt_tbl(nval, nrh, nbin, nband), &
             aero_sulf_tbl(nval, nrh, nband), &
             aero_bcar_tbl(nval, nband), &
             aero_bcar_rh_tbl(nval,nrh, nband), &
             aero_ocar_tbl(nval, nband), &
             aero_ocar_rh_tbl(nval, nrh, nband))

    ! Read LUT coefficients
    aero_dust_tbl = read_field(ncid, 'aero_dust_tbl',  nval, nbin, nband)
    aero_salt_tbl = read_field(ncid, 'aero_salt_tbl',  nval, nrh, nbin, nband)
    aero_sulf_tbl = read_field(ncid, 'aero_sulf_tbl',  nval, nrh, nband)
    aero_bcar_tbl = read_field(ncid, 'aero_bcar_tbl',  nval, nband)
    aero_bcar_rh_tbl = read_field(ncid, 'aero_bcar_rh_tbl',  nval, nrh, nband)
    aero_ocar_tbl = read_field(ncid, 'aero_ocar_tbl',  nval, nband)
    aero_ocar_rh_tbl = read_field(ncid, 'aero_ocar_rh_tbl',  nval, nrh, nband)

    ncid = nf90_close(ncid)

    call stop_on_err(aerosol_spec%load(band_lims_wvn, &
                                  merra_aero_bin_lims, aero_rh, &
                                  aero_dust_tbl, aero_salt_tbl, aero_sulf_tbl, &
                                  aero_bcar_tbl, aero_bcar_rh_tbl, &
                                  aero_ocar_tbl, aero_ocar_rh_tbl))

  end subroutine load_aerosol_optics
  !--------------------------------------------------------------------------------------------------------------------

end module mo_optics_utils_rrtmgp
