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
  use mo_rte_kind,           only: wp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics,         only: ty_gas_optics
  ! --------------------------------------------------
  use netcdf
  implicit none
  private
  public :: load_and_init

  interface read_field
    module procedure read_scalar, read_1d_field, read_2d_field, read_3d_field, read_4d_field
  end interface
contains
  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_and_init(kdist, filename, available_gases)
    class(ty_gas_optics), intent(inout) :: kdist
    character(len=*),     intent(in   ) :: filename
    class(ty_gas_concs),  intent(in   ) :: available_gases ! Which gases does the host model have available?
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

    character(len=256), dimension(:), allocatable :: gas_minor, identifier_minor
    character(len=256), dimension(:), allocatable :: minor_gases_lower,               minor_gases_upper
    integer, dimension(:,:),          allocatable :: minor_limits_gpt_lower,          minor_limits_gpt_upper
    logical, dimension(:),            allocatable :: minor_scales_with_density_lower, minor_scales_with_density_upper
    character(len=256), dimension(:), allocatable :: scaling_gas_lower,               scaling_gas_upper
    logical, dimension(:),            allocatable :: scale_by_complement_lower,       scale_by_complement_upper
    integer, dimension(:),            allocatable :: kminor_start_lower,              kminor_start_upper
    real(wp), dimension(:,:,:),       allocatable :: kminor_lower,                    kminor_upper

    real(wp), dimension(:,:,:  ), allocatable :: rayl_lower, rayl_upper
    real(wp), dimension(:      ), allocatable :: solar_src
    real(wp), dimension(:,:    ), allocatable :: totplnk
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac
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
               ninternalSourcetemps
    ! --------------------------------------------------
    !
    ! How big are the various arrays?
    !
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("load_and_init(): can't open file " // trim(fileName))
    ntemps            = get_dim_length(ncid,'temperature')
    npress            = get_dim_length(ncid,'pressure')
    nabsorbers        = get_dim_length(ncid,'absorber')
    nminorabsorbers   = get_dim_length(ncid,'minor_absorber')
    nextabsorbers     = get_dim_length(ncid,'absorber_ext')
    nmixingfracs      = get_dim_length(ncid,'mixing_fraction')
    nlayers           = get_dim_length(ncid,'atmos_layer')
    nbnds             = get_dim_length(ncid,'bnd')
    ngpts             = get_dim_length(ncid,'gpt')
    npairs            = get_dim_length(ncid,'pair')
    nminor_absorber_intervals_lower &
                      = get_dim_length(ncid,'minor_absorber_intervals_lower')
    nminor_absorber_intervals_upper  &
                      = get_dim_length(ncid,'minor_absorber_intervals_upper')
    ninternalSourcetemps &
                      = get_dim_length(ncid,'temperature_Planck')
    ncontributors_lower = get_dim_length(ncid,'contributors_lower')
    ncontributors_upper = get_dim_length(ncid,'contributors_upper')
    ! -----------------
    !
    ! Read the many arrays
    !
    gas_names         = read_char_vec(ncid, 'gas_names', nabsorbers)
    key_species       = read_field(ncid, 'key_species',  2, nlayers, nbnds)
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
                      = read_field(ncid, 'kminor_start_lower', nminor_absorber_intervals_lower)
    kminor_start_upper &
                      = read_field(ncid, 'kminor_start_upper', nminor_absorber_intervals_upper)
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
                                  rayl_lower, rayl_upper))
    else
      !
      ! Solar source doesn't have an dependencies yet
      !
      solar_src = read_field(ncid, 'solar_source', ngpts)
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
                                  solar_src, &
                                  rayl_lower, rayl_upper))
    end if
    ! --------------------------------------------------
    ncid = nf90_close(ncid)
  end subroutine load_and_init
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function read_scalar(ncid, varName)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    real(wp)                     :: read_scalar

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_scalar)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_scalar
  !--------------------------------------------------------------------------------------------------------------------
  function read_1d_field(ncid, varName, nx)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx
    real(wp), dimension(nx)      :: read_1d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 1) /= [nx])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent.")
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_1d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_1d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_2d_field(ncid, varName, nx, ny)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny
    real(wp), dimension(nx, ny)  :: read_2d_field

    integer :: varid
    if(any(get_data_size(ncid, varName, 2) /= [nx, ny])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent.")
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_2d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_2d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_3d_field(ncid, varName, nx, ny, nz)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx, ny, nz
    real(wp), dimension(nx, ny, nz)  :: read_3d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 3) /= [nx, ny, nz])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent.")
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_3d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_3d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_4d_field(ncid, varName, nw, nx, ny, nz)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nw, nx, ny, nz
    real(wp), dimension(nw, nx, ny, nz)  :: read_4d_field

    integer :: varid

    if(any(get_data_size(ncid, varName, 4) /= [nw, nx, ny, nz])) &
      call stop_on_err("read_field: variable " // trim(varName) // " size is inconsistent." )
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_field: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_4d_field)  /= NF90_NOERR) &
      call stop_on_err("read_field: can't read variable " // trim(varName))

  end function read_4d_field
  !--------------------------------------------------------------------------------------------------------------------
  function read_logical_vec(ncid, varName, nx)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx
    integer,      dimension(nx) :: read_logical_tmp
    logical,      dimension(nx) :: read_logical_vec

    integer :: varid
    integer :: ix

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_logical_vec: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_logical_tmp)  /= NF90_NOERR) &
      call stop_on_err("read_logical_vec: can't read variable " // trim(varName))
    do ix = 1, nx
      if (read_logical_tmp(ix) .eq. 0) then
        read_logical_vec(ix) = .false.
      else
        read_logical_vec(ix) = .true.
      endif
    enddo

  end function read_logical_vec
  !--------------------------------------------------------------------------------------------------------------------
  function read_char_vec(ncid, varName, nx)
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: nx
    character(len=32), dimension(nx) :: read_char_vec

    integer :: varid

    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("read_char_vec: can't find variable " // trim(varName))
    if(nf90_get_var(ncid, varid, read_char_vec)  /= NF90_NOERR) &
      call stop_on_err("read_char_vec: can't read variable " // trim(varName))

  end function read_char_vec
  !--------------------------------------------------------------------------------------------------------------------
  function var_exists(ncid, varName)
    !
    ! Does this variable exist (have a valid var_id) in the open netCDF file?
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    logical :: var_exists

    integer :: varId
    var_exists = nf90_inq_varid(ncid, trim(varName), varid) == NF90_NOERR
  end function var_exists
  !--------------------------------------------------------------------------------------------------------------------
  function get_dim_length(ncid, dimname)
    !
    ! Get the length of a dimension from an open netCDF file
    !  This is unfortunately a two-step process
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: get_dim_length

    integer :: dimid

    if(nf90_inq_dimid(ncid, trim(dimname), dimid) == NF90_NOERR) then
      if(nf90_inquire_dimension(ncid, dimid, len=get_dim_length) /= NF90_NOERR) get_dim_length = 0
    else
      get_dim_length = 0
    end if

  end function get_dim_length
  !--------------------------------------------------------------------------------------------------------------------
  function get_data_size(ncid, varName, n)
    !
    ! Returns the extents of a netcdf variable on disk
    !
    integer,          intent(in) :: ncid
    character(len=*), intent(in) :: varName
    integer,          intent(in) :: n
    integer                      :: get_data_size(n)

    integer :: i
    integer :: varid, ndims, dimids(n)

    get_data_size(n) = -1
    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't find variable " // trim(varName))
    if(nf90_inquire_variable(ncid, varid, ndims = ndims) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't get information for variable " // trim(varName))
    if(ndims /= n) &
      call stop_on_err("get_data_size:  variable " // trim(varName) // " has the wrong number of dimensions" )
    if(nf90_inquire_variable(ncid, varid, dimids = dimids) /= NF90_NOERR) &
      call stop_on_err("get_data_size: can't read dimension ids for variable " // trim(varName))
    do i = 1, n
      if(nf90_inquire_dimension(ncid, dimids(i), len = get_data_size(i)) /= NF90_NOERR) &
        call stop_on_err("get_data_size: can't get dimension lengths for variable " // trim(varName))
    end do

  end function get_data_size
  !--------------------------------------------------------------------------------------------------------------------
end module
