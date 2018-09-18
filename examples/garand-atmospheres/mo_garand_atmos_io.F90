! This code is part of Radiative Transfer for Energetics (RTE) and
!   RRTM for GCM Applications - Parallel (RRTMGP)
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
! This module reads and writes to netCDF files (nominally for the Garand atmospheres) using
!   serial I/O. The files follow arbitrary conventions adopted by RTE+RRTMGP developers.
!   It may be useful as an example for other formats/conventions.
!
! Reading routines use "allocation on assignment," a feature of Fortran 2003 that may require
!   particular compilation flags.
!
module mo_garand_atmos_io
  !
  ! RTE+RRTMGP modules
  !
  use mo_rte_kind,           only: wp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_util_reorder,       only: reorder123x312
  use mo_optical_props,      only: ty_optical_props
  !
  ! NetCDF I/O routines, shared with other RTE+RRTMGP examples
  !
  use mo_simple_netcdf,      only: read_field, read_string, var_exists, get_dim_size, &
                                   write_field, create_dim, create_var

  use netcdf
  implicit none
  private

  public :: read_atmos, is_lw, is_sw,           &
            read_lw_bc, read_sw_bc, read_lw_rt, &
            write_spectral_disc,                &
            write_fluxes, write_dir_fluxes, write_heating_rates
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read profiles for all columns  -- T, p, and gas concentrations
  !   Allocation occurs on assignments (says the F2003 standard)
  !
  subroutine read_atmos(fileName,                          &
                        p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)
    character(len=*),   intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable,                 &
                        intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
    type(ty_gas_concs), intent(inout) :: gas_concs
    ! -------------------
    integer :: ncid, ncol, nlay, nlev
    ! -------------------

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    nlev = get_dim_size(ncid, 'lev')
    if(nlev /= nlay+1) call stop_on_err("read_atmos: nlev should be nlay+1")

    !
    ! These lines assume that compilers follow the Fortran 2003 standard for
    !   allocating on assignment. This may require explicit compiler support
    !   e.g. -assume realloc_lhs flag for Intel
    !
    p_lay = read_field(ncid, 'p_lay', ncol, nlay)
    t_lay = read_field(ncid, 't_lay', ncol, nlay)
    p_lev = read_field(ncid, 'p_lev', ncol, nlev)
    t_lev = read_field(ncid, 't_lev', ncol, nlev)

    if(var_exists(ncid, 'vmr_h2o')) call stop_on_err(gas_concs%set_vmr('h2o', read_field(ncid, 'vmr_h2o', ncol, nlay)))
    if(var_exists(ncid, 'vmr_co2')) call stop_on_err(gas_concs%set_vmr('co2', read_field(ncid, 'vmr_co2', ncol, nlay)))
    if(var_exists(ncid, 'vmr_o3' )) call stop_on_err(gas_concs%set_vmr('o3' , read_field(ncid, 'vmr_o3' , ncol, nlay)))
    if(var_exists(ncid, 'vmr_n2o')) call stop_on_err(gas_concs%set_vmr('n2o', read_field(ncid, 'vmr_n2o', ncol, nlay)))
    if(var_exists(ncid, 'vmr_co' )) call stop_on_err(gas_concs%set_vmr('co' , read_field(ncid, 'vmr_co' , ncol, nlay)))
    if(var_exists(ncid, 'vmr_ch4')) call stop_on_err(gas_concs%set_vmr('ch4', read_field(ncid, 'vmr_ch4', ncol, nlay)))
    if(var_exists(ncid, 'vmr_o2' )) call stop_on_err(gas_concs%set_vmr('o2' , read_field(ncid, 'vmr_o2' , ncol, nlay)))
    if(var_exists(ncid, 'vmr_n2' )) call stop_on_err(gas_concs%set_vmr('n2' , read_field(ncid, 'vmr_n2' , ncol, nlay)))

    ! col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
    if(var_exists(ncid, 'col_dry')) col_dry = read_field(ncid, 'col_dry', ncol, nlay)

    ncid = nf90_close(ncid)

  end subroutine read_atmos
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Does this file contain variables needed to do SW calculations ?
  !
  function is_sw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_sw
    ! -------------------
    integer :: ncid
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("is_sw: can't find file " // trim(fileName))

    is_sw = var_exists(ncid, 'solar_zenith_angle')
    ncid = nf90_close(ncid)
  end function is_sw

  ! ----------------------
  function is_lw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw

    is_lw = .not. is_sw(fileName)
  end function is_lw
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read LW boundary conditions for all columns
  !
  subroutine read_lw_bc(fileName, t_sfc, emis_sfc)
    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:),   allocatable, intent(inout) :: t_sfc
    real(wp), dimension(:,:), allocatable, intent(inout) :: emis_sfc
    ! -------------------
    integer :: ncid
    integer :: ncol, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_lw_bc: can't find file " // trim(fileName))

    ncol  = get_dim_size(ncid, 'col')
    nband = get_dim_size(ncid, 'band')

    t_sfc    =  read_field(ncid, 't_sfc',           ncol)
    emis_sfc =  read_field(ncid, 'emis_sfc', nband, ncol)

    ncid = nf90_close(ncid)
  end subroutine read_lw_bc
   !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read LW radiative transfer parameters (number of quadrature angle for no-scattering calculations)
  !
  subroutine read_lw_rt(fileName, n_quad_angles)
    character(len=*), intent(in   ) :: fileName
    integer,          intent(  out) :: n_quad_angles
    ! -------------------
    integer :: ncid, ncol, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_lw_bc: can't find file " // trim(fileName))
    n_quad_angles  = get_dim_size(ncid, 'angle')
    ncid = nf90_close(ncid)
  end subroutine read_lw_rt
 !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read SW boundary conditions for all columns
  !
  subroutine read_sw_bc(fileName, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:),   allocatable, intent(inout) :: sza, tsi
    real(wp), dimension(:,:), allocatable, intent(inout) :: sfc_alb_dir, sfc_alb_dif
    real(wp),                              intent(inout) :: tsi_scaling
    ! -------------------
    integer :: ncid, ncol, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_sw_bc: can't find file " // trim(fileName))

    ncol  = get_dim_size(ncid, 'col')
    nband = get_dim_size(ncid, 'band')

    sza         =  read_field(ncid, 'solar_zenith_angle',        ncol)
    tsi         =  read_field(ncid, 'total_solar_irradiance',    ncol)
    sfc_alb_dir =  read_field(ncid, 'sfc_alb_direct',  nband, ncol)
    sfc_alb_dif =  read_field(ncid, 'sfc_alb_diffuse', nband, ncol)

    ! read tsi_scaling only if variable is present in the netCDF file
    if(var_exists(ncid, 'tsi_scaling')) tsi_scaling = read_field(ncid, 'tsi_scaling' )

    ncid = nf90_close(ncid)
  end subroutine read_sw_bc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write broadband and by-band fluxes
  !
  subroutine write_fluxes(fileName, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) ::     flux_up,     flux_dn,     flux_net
    real(wp), dimension(:,:,:), optional, &
                                intent(in) :: bnd_flux_up, bnd_flux_dn, bnd_flux_net
    ! -------------------
    integer :: ncid, ncol, nlev, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlev  = get_dim_size(ncid, 'lev')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid,      "flux_up",          ["col",  "lev"],         [ncol, nlev])
    call create_var(ncid,      "flux_dn",          ["col",  "lev"],         [ncol, nlev])
    call create_var(ncid,      "flux_net",         ["col",  "lev"],         [ncol, nlev])
    if(present(bnd_flux_up )) call create_var(ncid, "band_flux_up",  ["band", "col ", "lev "], [nband, ncol, nlev])
    if(present(bnd_flux_dn )) call create_var(ncid, "band_flux_dn",  ["band", "col ", "lev "], [nband, ncol, nlev])
    if(present(bnd_flux_net)) call create_var(ncid, "band_flux_net", ["band", "col ", "lev "], [nband, ncol, nlev])

    call stop_on_err(write_field(ncid, "flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "flux_dn",  flux_dn ))
    call stop_on_err(write_field(ncid, "flux_net", flux_net))
    ! col,lay,bnd -> bnd,col,lay
    if(present(bnd_flux_up )) call stop_on_err(write_field(ncid, "band_flux_up",  reorder123x312(bnd_flux_up )))
    if(present(bnd_flux_dn )) call stop_on_err(write_field(ncid, "band_flux_dn",  reorder123x312(bnd_flux_dn )))
    if(present(bnd_flux_net)) call stop_on_err(write_field(ncid, "band_flux_net", reorder123x312(bnd_flux_net)))

    ncid = nf90_close(ncid)
  end subroutine write_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write direct-beam fluxes
  !
  subroutine write_dir_fluxes(fileName, flux_dir, bnd_flux_dir)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) ::     flux_dir
    real(wp), dimension(:,:,:), optional, &
                                intent(in) :: bnd_flux_dir
    ! -------------------
    integer :: ncid, ncol, nlay, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_dir_fluxes: can't open file " // trim(fileName))

    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid,      "flux_dir_dn",         ["col",  "lev"],         [ncol, nlay+1])
    if(present(bnd_flux_dir)) call create_var(ncid, "band_flux_dir_dn", ["band", "col ", "lev "], [nband, ncol, nlay+1])

    call stop_on_err(write_field(ncid, "flux_dir_dn",  flux_dir))
    if(present(bnd_flux_dir)) call stop_on_err(write_field(ncid, "band_flux_dir_dn",  reorder123x312(bnd_flux_dir)))

    ncid = nf90_close(ncid)
  end subroutine write_dir_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write heating rates (broadband, by-band)
  !
  subroutine write_heating_rates(fileName, heating_rate, bnd_heating_rate)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) ::     heating_rate
    real(wp), dimension(:,:,:), intent(in) :: bnd_heating_rate
    ! -------------------
    integer :: ncid, ncol, nlay, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_heating_rates: can't open file " // trim(fileName))

    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid,      "heating_rate",          ["col", "lay"],         [ncol, nlay])
    call create_var(ncid, "band_heating_rate", ["band", "col ", "lay "], [nband, ncol, nlay])

    call stop_on_err(write_field(ncid,     "heating_rate",                     heating_rate))
    call stop_on_err(write_field(ncid, "band_heating_rate", reorder123x312(bnd_heating_rate)))

    ncid = nf90_close(ncid)
  end subroutine write_heating_rates
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write spectral discretization
  !
  subroutine write_spectral_disc(fileName, spectral_disc)
    character(len=*),        intent(in) :: fileName
    class(ty_optical_props), intent(in) :: spectral_disc
    ! -------------------
    integer :: ncid, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_spectral_disc: can't open file " // trim(fileName))

    nband = spectral_disc%get_nband()
    call create_dim(ncid, 'band', nband)
    call create_dim(ncid, "pair", 2)

    call create_var(ncid, "band_lims_wvn", ["pair", "band"], [2, nband])
    call stop_on_err(write_field(ncid, "band_lims_wvn", spectral_disc%get_band_lims_wavenumber()))
    call create_var(ncid, "band_lims_gpt", ["pair", "band"], [2, nband], NF90_INT)
    call stop_on_err(write_field(ncid, "band_lims_gpt", spectral_disc%get_band_lims_gpoint()))

    ncid = nf90_close(ncid)
  end subroutine write_spectral_disc
  !--------------------------------------------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop
    !
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then
      write(error_unit,*) trim(msg)
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
end module mo_garand_atmos_io
