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

  !$acc data create(tmp, tmp_col)
  call stop_on_err(gas_concs_garand%get_vmr(name, tmp))
  !$acc kernels
  tmp_col(:) = tmp(1, :)
  !$acc end kernels

  call stop_on_err(gas_concs%set_vmr       (name, tmp_col))
  !$acc end data
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
  use mo_garand_atmos_io,    only: read_atmos, write_lw_fluxes, write_sw_fluxes
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  real(wp), dimension(:,:),   allocatable :: col_dry
  real(wp), dimension(:,:),   allocatable :: temp_array

  !
  ! Longwave only
  !
  real(wp), dimension(:,:),   allocatable :: t_lev
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  !
  ! Shortwave only
  !
  real(wp), dimension(:),     allocatable :: mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  !
  ! Source functions
  !
  !   Longwave
  type(ty_source_func_lw), save               :: lw_sources
  !   Shortwave
  real(wp), dimension(:,:), allocatable, save :: toa_flux
  !
  ! Clouds
  !
  real(wp), allocatable, dimension(:,:) :: lwp, iwp, rel, rei
  logical,  allocatable, dimension(:,:) :: cloud_mask
  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir
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

  integer  :: ncol, nlay, nbnd, ngpt
  integer  :: icol, ilay, ibnd, iloop, igas
  real(wp) :: rel_val, rei_val

  character(len=8) :: char_input
  integer  :: nUserArgs=0, nloops
  logical :: use_luts = .true., write_fluxes = .true.
  integer, parameter :: ngas = 8
  character(len=3), dimension(ngas) &
                     :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']

  character(len=256) :: input_file, k_dist_file, cloud_optics_file
  !
  ! Timing variables
  !
  integer(kind=8)              :: start, finish, start_all, finish_all, clock_rate
  real(wp)                     :: avg
  integer(kind=8), allocatable :: elapsed(:)
  !$omp threadprivate( lw_sources, toa_flux, flux_up, flux_dn, flux_dir )
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line for any file names, block size
  !
  ! rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc  128 1
  ! rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc  128 1
  nUserArgs = command_argument_count()
  nloops = 1
  if (nUserArgs <  4) call stop_on_err("Need to supply input_file k_distribution_file ncol.")
  if (nUserArgs >= 1) call get_command_argument(1,input_file)
  if (nUserArgs >= 2) call get_command_argument(2,k_dist_file)
  if (nUserArgs >= 3) call get_command_argument(3,cloud_optics_file)
  if (nUserArgs >= 4) then
    call get_command_argument(4, char_input)
    read(char_input, '(i8)') ncol
    if(ncol <= 0) call stop_on_err("Specify positive ncol.")
  end if
  if (nUserArgs >= 5) then
    call get_command_argument(5, char_input)
    read(char_input, '(i8)') nloops
    if(nloops <= 0) call stop_on_err("Specify positive nloops.")
  end if
  if (nUserArgs >  6) print *, "Ignoring command line arguments beyond the first five..."
  if(trim(input_file) == '-h' .or. trim(input_file) == "--help") then
    call stop_on_err("rrtmgp_clouds input_file absorption_coefficients_file cloud_optics_file ncol")
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
  !$acc enter data copyin(p_lay, p_lev, t_lay, t_lev)
  ! ----------------------------------------------------------------------------
  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)
  is_sw = k_dist%source_is_external()
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
  ! LW calculations neglect scattering; SW calculations use the 2-stream approximation
  !   Here we choose the right variant of optical_props.
  !
  if(is_sw) then
    allocate(ty_optical_props_2str::atmos)
    allocate(ty_optical_props_2str::clouds)
  else
    allocate(ty_optical_props_1scl::atmos)
    allocate(ty_optical_props_1scl::clouds)
  end if
  ! Clouds optical props are defined by band
  call stop_on_err(clouds%init(k_dist%get_band_lims_wavenumber()))
  !
  ! Allocate arrays for the optical properties themselves.
  !
  select type(atmos)
    class is (ty_optical_props_1scl)
      !$acc enter data copyin(atmos)
      call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
      !$acc enter data copyin(atmos) create(atmos%tau)
    class is (ty_optical_props_2str)
      call stop_on_err(atmos%alloc_2str( ncol, nlay, k_dist))
      !$acc enter data copyin(atmos) create(atmos%tau, atmos%ssa, atmos%g)
    class default
      call stop_on_err("rte_rrtmgp_clouds: Don't recognize the kind of optical properties ")
  end select
  select type(clouds)
    class is (ty_optical_props_1scl)
      call stop_on_err(clouds%alloc_1scl(ncol, nlay))
      !$acc enter data copyin(clouds) create(clouds%tau)
    class is (ty_optical_props_2str)
      call stop_on_err(clouds%alloc_2str(ncol, nlay))
      !$acc enter data copyin(clouds) create(clouds%tau, clouds%ssa, clouds%g)
    class default
      call stop_on_err("rte_rrtmgp_clouds: Don't recognize the kind of optical properties ")
  end select
  ! ----------------------------------------------------------------------------
  !  Boundary conditions depending on whether the k-distribution being supplied
  !   is LW or SW
  if(is_sw) then
    ! toa_flux is threadprivate
    !$omp parallel
    allocate(toa_flux(ncol, ngpt))
    !$omp end parallel
    !
    allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), mu0(ncol))
    !$acc enter data create(sfc_alb_dir, sfc_alb_dif, mu0)
    ! Ocean-ish values for no particular reason
    !$acc kernels
    sfc_alb_dir = 0.06_wp
    sfc_alb_dif = 0.06_wp
    mu0 = .86_wp
    !$acc end kernels
  else
    ! lw_sorces is threadprivate
    !$omp parallel
    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist))
    !$omp end parallel

    allocate(t_sfc(ncol), emis_sfc(nbnd, ncol))
    !$acc enter data create(t_sfc, emis_sfc)
    ! Surface temperature
    !$acc kernels
    t_sfc = t_lev(1, merge(nlay+1, 1, top_at_1))
    emis_sfc = 0.98_wp
    !$acc end kernels
  end if
  ! ----------------------------------------------------------------------------
  !
  ! Fluxes
  !
  !$omp parallel
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1))
  !$omp end parallel

  !$acc enter data create(flux_up, flux_dn)
  if(is_sw) then
    allocate(flux_dir(ncol,nlay+1))
    !$acc enter data create(flux_dir)
  end if
  !
  ! Clouds
  !
  allocate(lwp(ncol,nlay), iwp(ncol,nlay), &
           rel(ncol,nlay), rei(ncol,nlay), cloud_mask(ncol,nlay))
  !$acc enter data create(cloud_mask, lwp, iwp, rel, rei)

  ! Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
  !   and not very close to the ground (< 900 hPa), and
  !   put them in 2/3 of the columns since that's roughly the
  !   total cloudiness of earth
  rel_val = 0.5 * (cloud_optics%get_min_radius_liq() + cloud_optics%get_max_radius_liq())
  rei_val = 0.5 * (cloud_optics%get_min_radius_ice() + cloud_optics%get_max_radius_ice())
  !$acc parallel loop collapse(2) copyin(t_lay) copyout(lwp, iwp, rel, rei)
  do ilay=1,nlay
    do icol=1,ncol
      cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100._wp * 100._wp .and. &
                              p_lay(icol,ilay) < 900._wp * 100._wp .and. &
                              mod(icol, 3) /= 0
      !
      ! Ice and liquid will overlap in a few layers
      !
      lwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) .and. t_lay(icol,ilay) > 263._wp)
      iwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) .and. t_lay(icol,ilay) < 273._wp)
      rel(icol,ilay) = merge(rel_val, 0._wp, lwp(icol,ilay) > 0._wp)
      rei(icol,ilay) = merge(rei_val, 0._wp, iwp(icol,ilay) > 0._wp)
    end do
  end do
  !$acc exit data delete(cloud_mask)
  ! ----------------------------------------------------------------------------
  !
  ! Multiple iterations for big problem sizes, and to help identify data movement
  !   For CPUs we can introduce OpenMP threading over loop iterations
  !
  allocate(elapsed(nloops))
  !
  call system_clock(start_all)
  !
  !$omp parallel do firstprivate(fluxes)
  do iloop = 1, nloops
    call system_clock(start)
    call stop_on_err(                                      &
      cloud_optics%cloud_optics(lwp, iwp, rel, rei, clouds))
    !
    ! Solvers
    !
    fluxes%flux_up => flux_up(:,:)
    fluxes%flux_dn => flux_dn(:,:)
    if(is_lw) then
      !$acc enter data create(lw_sources, lw_sources%lay_source, lw_sources%lev_source_inc, lw_sources%lev_source_dec, lw_sources%sfc_source)
      call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                         t_lay, t_sfc, &
                                         gas_concs,    &
                                         atmos,        &
                                         lw_sources,   &
                                         tlev = t_lev))
      call stop_on_err(clouds%increment(atmos))
      call stop_on_err(rte_lw(atmos, top_at_1, &
                              lw_sources,      &
                              emis_sfc,        &
                              fluxes))
      !$acc exit data delete(lw_sources%lay_source, lw_sources%lev_source_inc, lw_sources%lev_source_dec, lw_sources%sfc_source, lw_sources)
    else
      !$acc enter data create(toa_flux)
      fluxes%flux_dn_dir => flux_dir(:,:)

      call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                         t_lay,        &
                                         gas_concs,    &
                                         atmos,        &
                                         toa_flux))
      call stop_on_err(clouds%delta_scale())
      call stop_on_err(clouds%increment(atmos))
      call stop_on_err(rte_sw(atmos, top_at_1, &
                              mu0,   toa_flux, &
                              sfc_alb_dir, sfc_alb_dif, &
                              fluxes))
      !$acc exit data delete(toa_flux)
    end if
    !print *, "******************************************************************"
    call system_clock(finish, clock_rate)
    elapsed(iloop) = finish - start
  end do
  !
  call system_clock(finish_all, clock_rate)
  !
  !$acc exit data delete(lwp, iwp, rel, rei)
  !$acc exit data delete(p_lay, p_lev, t_lay, t_lev)

#ifdef _OPENACC
  avg = sum( elapsed(merge(2,1,nloops>1):) ) / real(merge(nloops-1,nloops,nloops>1))

  print *, "Execution times - min(s)        :", minval(elapsed) / real(clock_rate)
  print *, "                - avg(s)        :", avg / real(clock_rate)
  print *, "                - per column(ms):", avg / real(ncol) / (1.0e-3*clock_rate)
#else
  print *, "Execution times - total(s)      :", (finish_all-start_all) / real(clock_rate)
  print *, "                - per column(ms):", (finish_all-start_all) / real(ncol*nloops) / (1.0e-3*clock_rate)
#endif

  if(is_lw) then
    !$acc exit data copyout(flux_up, flux_dn)
    if(write_fluxes) call write_lw_fluxes(input_file, flux_up, flux_dn)
    !$acc exit data delete(t_sfc, emis_sfc)
  else
    !$acc exit data copyout(flux_up, flux_dn, flux_dir)
    if(write_fluxes) call write_sw_fluxes(input_file, flux_up, flux_dn, flux_dir)
    !$acc exit data delete(sfc_alb_dir, sfc_alb_dif, mu0)
  end if
  !$acc enter data create(lwp, iwp, rel, rei)
end program rte_rrtmgp_clouds
