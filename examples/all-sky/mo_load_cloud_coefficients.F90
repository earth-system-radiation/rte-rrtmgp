module mo_load_cloud_coefficients
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_cloud_optics_rrtmgp, &
                        only: ty_cloud_optics_rrtmgp
  use mo_simple_netcdf, only: read_field, read_string, dim_exists, var_exists, &
                              get_dim_size, write_field, create_dim, create_var
  use netcdf

  implicit none
  private
  public :: load_cld_lutcoeff
  ! ----------------------------------------------------------------------------------

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read cloud optical property LUT coefficients from NetCDF file
  !
  subroutine load_cld_lutcoeff(cloud_spec, cld_coeff_file, ngpnt, nspec, &
                               band_lims_gpt, band_lims_wvn)
    class(ty_cloud_optics_rrtmgp),         intent(inout) :: cloud_spec
    character(len=*),                      intent(in   ) :: cld_coeff_file
    integer,                               intent(in   ) :: ngpnt
    integer,                               intent(  out) :: nspec
    integer,  dimension(:,:), allocatable, intent(  out) :: band_lims_gpt
    real(wp), dimension(:,:), allocatable, intent(  out) :: band_lims_wvn
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrghice, nsize_liq, nsize_ice

    ! Lookup table interpolation constants
    real(wp) :: radliq_lwr          ! liquid particle size lower bound for interpolation
    real(wp) :: radliq_upr          ! liquid particle size upper bound for interpolation
    real(wp) :: radice_lwr          ! ice particle size lower bound for interpolation
    real(wp) :: radice_upr          ! ice particle size upper bound for interpolation
    ! LUT coefficients
    real(wp), dimension(:,:),   allocatable :: extliq   ! extinction: liquid
    real(wp), dimension(:,:),   allocatable :: ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:),   allocatable :: asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:), allocatable :: extice   ! extinction: ice
    real(wp), dimension(:,:,:), allocatable :: ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:), allocatable :: asyice   ! asymmetry parameter: ice

    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cld_lutcoeff(): can't open file " // trim(cld_coeff_file))

    ! Read band LUT coefficient dimensions
    nband = get_dim_size(ncid,'nband')
    if (dim_exists(ncid, 'ngpt')) then 
       nspec = get_dim_size(ncid,'ngpt')
    else
       nspec = nband
    endif
    nrghice   = get_dim_size(ncid,'nrghice')
    nsize_liq = get_dim_size(ncid,'nsize_liq')
    nsize_ice = get_dim_size(ncid,'nsize_ice')

    ! Read band wavenumber limits
    allocate(band_lims_gpt(2, nband), &
             band_lims_wvn(2, nband))
    band_lims_gpt = read_field(ncid, 'bnd_limits_gpt', 2, nband)
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)

    ! Read LUT constants
    radliq_lwr = read_field(ncid, 'radliq_lwr')
    radliq_upr = read_field(ncid, 'radliq_upr')
    radice_lwr = read_field(ncid, 'radice_lwr')
    radice_upr = read_field(ncid, 'radice_upr')

    ! Allocate cloud property band lookup table input arrays
    allocate(extliq(nsize_liq, nspec), &
             ssaliq(nsize_liq, nspec), &
             asyliq(nsize_liq, nspec), &
             extice(nsize_ice, nspec, nrghice), &
             ssaice(nsize_ice, nspec, nrghice), &
             asyice(nsize_ice, nspec, nrghice))

    if (dim_exists(ncid, 'ngpt')) then 
    ! Read g-point LUT coefficients
       extliq = read_field(ncid, 'gpt_extliq',  nsize_liq, nspec)
       ssaliq = read_field(ncid, 'gpt_ssaliq',  nsize_liq, nspec)
       asyliq = read_field(ncid, 'gpt_asyliq',  nsize_liq, nspec)
       extice = read_field(ncid, 'gpt_extice',  nsize_ice, nspec, nrghice)
       ssaice = read_field(ncid, 'gpt_ssaice',  nsize_ice, nspec, nrghice)
       asyice = read_field(ncid, 'gpt_asyice',  nsize_ice, nspec, nrghice)
    else
    ! Read band LUT coefficients
       extliq = read_field(ncid, 'bnd_extliq',  nsize_liq, nspec)
       ssaliq = read_field(ncid, 'bnd_ssaliq',  nsize_liq, nspec)
       asyliq = read_field(ncid, 'bnd_asyliq',  nsize_liq, nspec)
       extice = read_field(ncid, 'bnd_extice',  nsize_ice, nspec, nrghice)
       ssaice = read_field(ncid, 'bnd_ssaice',  nsize_ice, nspec, nrghice)
       asyice = read_field(ncid, 'bnd_asyice',  nsize_ice, nspec, nrghice)
    endif

    ncid = nf90_close(ncid)
    call stop_on_err(cloud_spec%load_cloud_optics(ngpnt, nband, nspec, &
                                                  band_lims_gpt, band_lims_wvn, &
                                                  radliq_lwr, radliq_upr, &
                                                  radice_lwr, radice_upr, &
                                                  extliq, ssaliq, asyliq, &
                                                  extice, ssaice, asyice))
  end subroutine load_cld_lutcoeff

  ! -----------------------------------------------------------------------------------
    subroutine stop_on_err(msg)
      !
      ! Print error message and stop
      !
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      if(len_trim(msg) > 0) then
        write (error_unit,*) trim(msg)
        error stop 1
      end if
    end subroutine

end module
