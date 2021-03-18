subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rte_rrtmgp_clouds stopping"
    stop
  end if
end subroutine stop_on_err

subroutine vmr_2d_to_1d(gas_concs, gas_concs_garand, name, sz1, sz2)
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_rte_kind,           only: wp

  type(ty_gas_concs), intent(in)    :: gas_concs_garand
  type(ty_gas_concs), intent(inout) :: gas_concs
  character(len=*),   intent(in)    :: name
  integer,            intent(in)    :: sz1, sz2

  real(wp) :: tmp(sz1, sz2), tmp_col(sz2)

  call stop_on_err(gas_concs_garand%get_vmr(name, tmp))
  tmp_col(:) = tmp(1, :)

  call stop_on_err(gas_concs%set_vmr       (name, tmp_col))
end subroutine vmr_2d_to_1d
! ----------------------------------------------------------------------------------
program rte_rrtmgp_clouds
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_load_cloud_coefficients, &
                             only: load_cld_lutcoeff, load_cld_padecoeff
  use mo_garand_atmos_io,    only: read_atmos
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  real(wp), dimension(:,:),   allocatable :: col_dry
  real(wp), dimension(:,:),   allocatable :: temp_array

  real(wp), dimension(:,:),   allocatable :: t_lev
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  !
  ! Source functions
  !
  type(ty_source_func_lw), save               :: lw_sources
  !
  ! Clouds
  !
  real(wp), allocatable, dimension(:,:) :: iwp, rei, lwp, rel
  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_optics_rrtmgp) :: k_dist
  type(ty_cloud_optics)      :: cloud_optics
  type(ty_gas_concs)         :: gas_concs, gas_concs_garand, gas_concs_1col
  class(ty_optical_props_arry), &
                 allocatable :: atmos, clouds
  type(ty_fluxes_broadband)  :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: top_at_1, is_sw, is_lw

  integer  :: nlay, nbnd, ngpt
  integer  :: icol, ilay, ibnd, igas
  real(wp) :: rei_val

  character(len=8) :: char_input
  integer :: nUserArgs=0, cloud_lay
  logical :: use_luts = .true.
  integer, parameter :: ngas = 8
  character(len=3), dimension(ngas) &
                     :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']

  character(len=256) :: input_file, k_dist_file, cloud_optics_file

  integer, parameter :: ncol = 5
  real(wp), dimension(ncol), parameter :: cloud_taus = [0.01_wp, 0.1_wp, 1._wp, 10._wp, 100._wp]
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line for any file names, block size
  !
  nUserArgs = command_argument_count()
  if (nUserArgs <  3) call stop_on_err("Need to supply input_file k_distribution_file ncol.")
  if (nUserArgs >  3) print *, "Ignoring command line arguments beyond the first three..."
  if (nUserArgs >= 1) call get_command_argument(1,input_file)
  if (nUserArgs >= 2) call get_command_argument(2,k_dist_file)
  if (nUserArgs >= 3) call get_command_argument(3,cloud_optics_file)
  if(trim(input_file) == '-h' .or. trim(input_file) == "--help") then
    call stop_on_err("cloud_scattering input_file absorption_coefficients_file cloud_optics_file")
  end if
  !
  ! Read temperature, pressure, gas concentrations.
  !   Arrays are allocated as they are read
  !
  call read_atmos(input_file,                 &
                  p_lay, t_lay, p_lev, t_lev, &
                  gas_concs_garand, col_dry)
  deallocate(col_dry)
  nlay = size(p_lay, 2)
  ! For clouds we'll use the first column, repeated over and over
  call stop_on_err(gas_concs%init(gas_names))
  do igas = 1, ngas
    call vmr_2d_to_1d(gas_concs, gas_concs_garand, gas_names(igas), size(p_lay, 1), nlay)
  end do
  !  If we trusted in Fortran allocate-on-assign we could skip the temp_array here
  allocate(temp_array(ncol, nlay))
  temp_array = spread(p_lay(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, p_lay)
  allocate(temp_array(ncol, nlay))
  temp_array = spread(t_lay(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, t_lay)
  allocate(temp_array(ncol, nlay+1))
  temp_array = spread(p_lev(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, p_lev)
  allocate(temp_array(ncol, nlay+1))
  temp_array = spread(t_lev(1,:), dim = 1, ncopies=ncol)
  call move_alloc(temp_array, t_lev)
  ! This puts pressure and temperature arrays on the GPU
  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)
  is_sw = k_dist%source_is_external()
  if(is_sw) call stop_on_err("Only checking out LW cloud scattering")
  is_lw = .not. is_sw
  !
  ! Should also try with Pade calculations
  !  call load_cld_padecoeff(cloud_optics, cloud_optics_file)
  !
  if(use_luts) then
    call load_cld_lutcoeff (cloud_optics, cloud_optics_file)
  else
    call load_cld_padecoeff(cloud_optics, cloud_optics_file)
  end if
  call stop_on_err(cloud_optics%set_ice_roughness(2))
  ! ----------------------------------------------------------------------------
  !
  ! Problem sizes
  !
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

  ! ----------------------------------------------------------------------------
  !  Boundary conditions
  call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist))
  allocate(t_sfc(ncol), emis_sfc(nbnd, ncol))
  ! Surface temperature
  t_sfc = t_lev(1, merge(nlay+1, 1, top_at_1))
  emis_sfc = 0.98_wp
  ! ----------------------------------------------------------------------------
  !
  ! Cloud physical descrpition
  !
  allocate(lwp(ncol,nlay), iwp(ncol,nlay), &
           rel(ncol,nlay), rei(ncol,nlay))
  rel = 0._wp; lwp = 0._wp
  rei = 0._wp; iwp = 0._wp

  cloud_lay = minloc(abs(p_lay(1,:) - 25000._wp), dim=1) ! Put the clouds near 250 hPa
  print *, "cloud layer is ", cloud_lay, ", p = ", p_lay(1,cloud_lay)
  iwp(:, cloud_lay) = 67._wp ! Gives an optical near of order 1
  rei(:, cloud_lay) = 0.5 * (cloud_optics%get_min_radius_ice() + cloud_optics%get_max_radius_ice())
  ! ----------------------------------------------------------------------------
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1))
  fluxes%flux_up => flux_up(:,:)
  fluxes%flux_dn => flux_dn(:,:)

  !
  ! No-scattering solutions
  !
  call make_optical_props_1scl
  call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                     t_lay, t_sfc, &
                                     gas_concs,    &
                                     atmos,        &
                                     lw_sources,   &
                                     tlev = t_lev))
  call stop_on_err(rte_lw(atmos, top_at_1, &
                          lw_sources,      &
                          emis_sfc,        &
                          fluxes))
  call write_fluxes("cloud_scattering.nc", flux_up, flux_dn, "clear_1scl")

  call stop_on_err(cloud_optics%cloud_optics(lwp, iwp, rel, rei, clouds))
  print *, "cloud taus ", clouds%tau(1, cloud_lay, :)
  do ibnd = 1, ncol
      clouds%tau(:, cloud_lay, ibnd) = cloud_taus * clouds%tau(:, cloud_lay, ibnd)
  end do

  call stop_on_err(clouds%increment(atmos))
  call stop_on_err(rte_lw(atmos, top_at_1, &
                          lw_sources,      &
                          emis_sfc,        &
                          fluxes))
  call write_fluxes("cloud_scattering.nc", flux_up, flux_dn, "cld_1scl")


  !
  ! No-scattering solutions
  !
  call make_optical_props_2str
  call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                     t_lay, t_sfc, &
                                     gas_concs,    &
                                     atmos,        &
                                     lw_sources,   &
                                     tlev = t_lev))
  call stop_on_err(rte_lw(atmos, top_at_1, &
                          lw_sources,      &
                          emis_sfc,        &
                          fluxes))
  call write_fluxes("cloud_scattering.nc", flux_up, flux_dn, "clear_2str")

  call stop_on_err(cloud_optics%cloud_optics(lwp, iwp, rel, rei, clouds))
  do ibnd = 1, ncol
      clouds%tau(:, cloud_lay, ibnd) = cloud_taus * clouds%tau(:, cloud_lay, ibnd)
  end do

  call stop_on_err(clouds%increment(atmos))
  call stop_on_err(rte_lw(atmos, top_at_1, &
                          lw_sources,      &
                          emis_sfc,        &
                          fluxes))
  call write_fluxes("cloud_scattering.nc", flux_up, flux_dn, "cld_tang")
  call stop_on_err(rte_lw(atmos, top_at_1, &
                          lw_sources,      &
                          emis_sfc,        &
                          fluxes, use_2stream = .true.))
  call write_fluxes("cloud_scattering.nc", flux_up, flux_dn, "cld_2str")
contains
  ! ----------------------------------------------------------------------------
  subroutine make_optical_props_1scl
    if(allocated(atmos )) deallocate(atmos)
    allocate(ty_optical_props_1scl::atmos)
    if(allocated(clouds )) deallocate(clouds)
    allocate(ty_optical_props_1scl::clouds)
    ! Clouds optical props are defined by band
    !
    ! Allocate arrays for the optical properties themselves.
    !
    select type(atmos)
      class is (ty_optical_props_1scl)
        call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
    select type(clouds)
      class is (ty_optical_props_1scl)
        call stop_on_err(clouds%alloc_1scl(ncol, nlay, cloud_optics))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_1scl
  ! ----------------------------------------------------------------------------
  subroutine make_optical_props_2str
    if(allocated(atmos )) deallocate(atmos)
    allocate(ty_optical_props_2str::atmos)
    if(allocated(clouds )) deallocate(clouds)
    allocate(ty_optical_props_2str::clouds)
    ! Clouds optical props are defined by band
    !
    ! Allocate arrays for the optical properties themselves.
    !
    select type(atmos)
      class is (ty_optical_props_2str)
        call stop_on_err(atmos%alloc_2str(ncol, nlay, k_dist))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
    select type(clouds)
      class is (ty_optical_props_2str)
        call stop_on_err(clouds%alloc_2str(ncol, nlay, cloud_optics))
      class default
        call stop_on_err("rte_rrtmgp_atmos: Don't recognize the kind of optical properties ")
    end select
  end subroutine make_optical_props_2str
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_fluxes(fileName, flux_up, flux_dn, id)
    use netcdf
    use mo_simple_netcdf, only: get_dim_size, create_dim, create_var, write_field
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn
    character(len=*),         intent(in) :: id
    ! -------------------
    integer :: ncid, ncol, nlev
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(flux_up, dim=1)
    nlev  = get_dim_size(ncid, 'lev')
    call create_dim(ncid, "col_flx", ncol)

    call create_var(ncid, trim(id) // "_up",  ["col_flx",  "lev    "], [ncol, nlev])
    call create_var(ncid, trim(id) // "_dn",  ["col_flx",  "lev    "], [ncol, nlev])

    call stop_on_err(write_field(ncid, trim(id) // "_up",  flux_up ))
    call stop_on_err(write_field(ncid, trim(id) // "_dn",  flux_dn ))

    ncid = nf90_close(ncid)
  end subroutine write_fluxes
  !--------------------------------------------------------------------------------------------------------------------
end program rte_rrtmgp_clouds
