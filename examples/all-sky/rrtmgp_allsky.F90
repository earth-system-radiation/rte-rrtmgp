program rte_rrtmgp_allsky
  use, intrinsic :: iso_fortran_env, & 
                             only: output_unit
  use mo_rte_kind,           only: wp, i8, wl
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics_rrtmgp,only: ty_cloud_optics_rrtmgp
  use mo_aerosol_optics_rrtmgp_merra ! Includes aerosol type integers
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_load_cloud_coefficients, &
                             only: load_cld_lutcoeff, load_cld_padecoeff
  use mo_load_aerosol_coefficients, &
                             only: load_aero_lutcoeff
  use mo_rte_config,         only: rte_config_checks
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev ! t_lev is only needed for LW
  real(wp), dimension(:,:),   allocatable :: q, o3, col_dry
  real(wp), dimension(:,:),   allocatable :: temp_array

  !
  ! Longwave only
  !
  real(wp), dimension(:),     allocatable :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  !
  ! Shortwave only
  !
  real(wp), dimension(:),     allocatable :: mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  !
  ! Gas concentrations 
  !
  character(len=3), dimension(8), parameter :: &
                   gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']
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
  ! Aerosols
  !
  logical :: cell_has_aerosols
  integer,  dimension(:,:), allocatable :: aero_type 
                                           ! MERRA2/GOCART aerosol type
  real(wp), dimension(:,:), allocatable :: aero_size
                                           ! Aerosol size for dust and sea salt
  real(wp), dimension(:,:), allocatable :: aero_mass
                                           ! Aerosol mass column (kg/m2)
  real(wp), dimension(:,:), allocatable :: relhum
                                           ! Relative humidity (fraction)
  logical, dimension(:,:), allocatable  :: aero_mask
                                           ! Aerosol mask

  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_optics_rrtmgp)   :: k_dist
  type(ty_cloud_optics_rrtmgp) :: cloud_optics
  type(ty_aerosol_optics_rrtmgp_merra)   & 
                               :: aerosol_optics
  type(ty_gas_concs)           :: gas_concs, gas_concs_garand, gas_concs_1col
  class(ty_optical_props_arry), &
                 allocatable   :: atmos, clouds, aerosols
  type(ty_fluxes_broadband)    :: fluxes

  !
  ! Inputs to RRTMGP
  !
  logical :: top_at_1, is_sw, is_lw

  integer  :: nbnd, ngpt
  integer  :: icol, ilay, ibnd, iloop, igas

  character(len=8) :: char_input
  integer :: nUserArgs, nloops, ncol, nlay
  ! logical :: write_fluxes = .false.
  logical :: do_clouds = .false., use_luts = .true. 
  logical :: do_aerosols = .false.
  integer, parameter :: ngas = 8

  character(len=256) :: output_file, k_dist_file, cloud_optics_file, aerosol_optics_file
  !
  ! Timing variables
  !
  integer(kind=i8)              :: start, finish, start_all, finish_all, clock_rate
  real(wp)                      :: avg, mint
  integer(kind=i8), allocatable :: elapsed(:)
  ! NAR OpenMP CPU directives in compatible with OpenMP GPU directives
  !!$omp threadprivate( lw_sources, toa_flux, flux_up, flux_dn, flux_dir )
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !
  ! Parse command line: rrtmgp_allsky ncol nlay nreps kdist [clouds [aerosols]] 

  !
  nUserArgs = command_argument_count()
  if (nUserArgs <  5) call stop_on_err("Usage: rrtmgp_allsky ncol nlay nreps output_file gas-optics [cloud-optics [aerosol-optics]]")

  call get_command_argument(1, char_input)
  read(char_input, '(i8)') ncol
  if(ncol <= 0) call stop_on_err("Specify positive ncol.")

  call get_command_argument(2, char_input)
  read(char_input, '(i8)') nlay
  if(nlay <= 0) call stop_on_err("Specify positive nlay.")

  call get_command_argument(3, char_input)
  read(char_input, '(i8)') nloops
  if(nloops <= 0) call stop_on_err("Specify positive nreps (number of times to do ncol examples.")

  call get_command_argument(4,output_file)
  call get_command_argument(5,k_dist_file)

  if (nUserArgs >= 6) then 
    call get_command_argument(6,cloud_optics_file)
    do_clouds = .true. 
  end if 
  if (nUserArgs >= 7) then 
    call get_command_argument(7,aerosol_optics_file)
    do_aerosols = .true. 
  end if 
  if (nUserArgs >  7) print *, "Ignoring command line arguments beyond the first seven..."
  ! -----------------------------------------------------------------------------------
  allocate(p_lay(ncol, nlay), t_lay(ncol, nlay), p_lev(ncol, nlay+1), t_lev(ncol, nlay+1))
  allocate(q    (ncol, nlay),    o3(ncol, nlay))
  !$acc        data create(   p_lay, t_lay, p_lev, t_lev, q, o3)
  !$omp target data map(alloc:p_lay, t_lay, p_lev, t_lev, q, o3)
  call compute_profiles(300._wp, ncol, nlay, p_lay, t_lay, p_lev, t_lev, q, o3)

  call stop_on_err(gas_concs%init(gas_names))
  call stop_on_err(gas_concs%set_vmr("h2o", q )) 
  call stop_on_err(gas_concs%set_vmr("o3",  o3)) 
  call stop_on_err(gas_concs%set_vmr("co2", 348.e-6_wp)) 
  call stop_on_err(gas_concs%set_vmr("ch4", 1650.e-9_wp)) 
  call stop_on_err(gas_concs%set_vmr("n2o", 306.e-9_wp)) 
  call stop_on_err(gas_concs%set_vmr("n2",  0.7808_wp)) 
  call stop_on_err(gas_concs%set_vmr("o2",  0.2095_wp)) 
  call stop_on_err(gas_concs%set_vmr("co",  0._wp)) 
  ! ----------------------------------------------------------------------------
  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)
  is_sw = k_dist%source_is_external()
  is_lw = .not. is_sw
  if (do_clouds) then 
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
  end if

  if (do_aerosols) then 
    !
    ! Load aerosol optics coefficients from lookup tables
    !
    call load_aero_lutcoeff (aerosol_optics, aerosol_optics_file)
  end if 

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
  else
    allocate(ty_optical_props_1scl::atmos)
  end if
  !
  ! Allocate arrays for the optical properties themselves.
  !
  select type(atmos)
    class is (ty_optical_props_1scl)
      !$acc enter data copyin(atmos)
      call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist))
      !$acc enter data copyin(atmos) create(atmos%tau)
      !$omp target enter data map(alloc:atmos%tau)
    class is (ty_optical_props_2str)
      call stop_on_err(atmos%alloc_2str( ncol, nlay, k_dist))
      !$acc enter data copyin(atmos) create(atmos%tau, atmos%ssa, atmos%g)
      !$omp target enter data map(alloc:atmos%tau, atmos%ssa, atmos%g)
    class default
      call stop_on_err("rte_rrtmgp_allsky: Don't recognize the kind of optical properties ")
  end select

  ! ----------------------------------------------------------------------------
  !  Boundary conditions depending on whether the k-distribution being supplied
  !   is LW or SW
  if(is_sw) then
    ! toa_flux is threadprivate
    !!$omp parallel
    allocate(toa_flux(ncol, ngpt))
    !!$omp end parallel
    !
    allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), mu0(ncol))
    !$acc         enter data create(   sfc_alb_dir, sfc_alb_dif, mu0)
    !$omp target  enter data map(alloc:sfc_alb_dir, sfc_alb_dif, mu0)
    ! Ocean-ish values for no particular reason
    !$acc kernels
    !$omp target
    sfc_alb_dir = 0.06_wp
    sfc_alb_dif = 0.06_wp
    mu0 = .86_wp
    !$acc end kernels
    !$omp end target
  else
    ! lw_sorces is threadprivate
    !!$omp parallel
    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist))
    !!$omp end parallel

    allocate(t_sfc(ncol), emis_sfc(nbnd, ncol))
    !$acc         enter data create   (t_sfc, emis_sfc)
    !$omp target  enter data map(alloc:t_sfc, emis_sfc)
    ! Surface temperature
    !$acc kernels
    !$omp target
    t_sfc = t_lev(1, merge(nlay+1, 1, top_at_1))
    emis_sfc = 0.98_wp
    !$acc end kernels
    !$omp end target
  end if
  ! ----------------------------------------------------------------------------
  !
  ! Fluxes
  !
  !!$omp parallel
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1))
  !!$omp end parallel

  !$acc         data create(   flux_up, flux_dn)
  !$omp target  data map(alloc:flux_up, flux_dn)
  if(is_sw) then
    allocate(flux_dir(ncol,nlay+1))
    !$acc enter data create(flux_dir)
    !$omp target enter data map(alloc:flux_dir)
  end if

  if (do_clouds)   call compute_clouds
  if (do_aerosols) call compute_aerosols

  ! ----------------------------------------------------------------------------
  !
  ! Multiple iterations for big problem sizes, and to help identify data movement
  !   For CPUs we can introduce OpenMP threading over loop iterations
  !
  allocate(elapsed(nloops))
  !
  call system_clock(start_all)
  !
  !!$omp parallel do firstprivate(fluxes)
  do iloop = 1, nloops
    ! Omit the checks starting with the second iteration
    if (iloop > 1) call rte_config_checks(logical(.false., wl))

    call system_clock(start)
    !
    ! Cloud optics
    !
    if(do_clouds) & 
      call stop_on_err(cloud_optics%cloud_optics(lwp, iwp, rel, rei, clouds))
    !
    ! Aerosol optics
    !
    if(do_aerosols) & 
      call stop_on_err(aerosol_optics%aerosol_optics(aero_type, aero_size,  &
                                                     aero_mass, relhum, aerosols))
    !
    ! Solvers
    !
    fluxes%flux_up => flux_up(:,:)
    fluxes%flux_dn => flux_dn(:,:)
    if(is_lw) then
      !
      ! Should we allocate these once, rather than once per loop? They're big. 
      ! 
      !$acc        data create(   lw_sources, lw_sources%lay_source,     lw_sources%lev_source_inc) &
      !$acc             create(               lw_sources%lev_source_dec, lw_sources%sfc_source,     lw_sources%sfc_source_Jac)
      !$omp target data map(alloc:            lw_sources%lay_source,     lw_sources%lev_source_inc) &
      !$omp             map(alloc:            lw_sources%lev_source_dec, lw_sources%sfc_source,     lw_sources%sfc_source_Jac)
      call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                         t_lay, t_sfc, &
                                         gas_concs,    &
                                         atmos,        &
                                         lw_sources,   &
                                         tlev = t_lev))
      if(do_clouds)   call stop_on_err(clouds%increment(atmos))
      if(do_aerosols) call stop_on_err(aerosols%increment(atmos))
      call stop_on_err(rte_lw(atmos, top_at_1, &
                              lw_sources,      &
                              emis_sfc,        &
                              fluxes))
      !$acc        end data
      !$omp end target data

    else
      !$acc         data create(   toa_flux)
      !$omp target  data map(alloc:toa_flux)
      fluxes%flux_dn_dir => flux_dir(:,:)

      call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                         t_lay,        &
                                         gas_concs,    &
                                         atmos,        &
                                         toa_flux))
      if(do_clouds) then 
        call stop_on_err(clouds%delta_scale())
        call stop_on_err(clouds%increment(atmos))
      end if 
      if(do_aerosols) then 
        call stop_on_err(aerosols%delta_scale())
        call stop_on_err(aerosols%increment(atmos))
      end if 
      call stop_on_err(rte_sw(atmos, top_at_1, &
                              mu0,   toa_flux, &
                              sfc_alb_dir, sfc_alb_dif, &
                              fluxes))
      !$acc        end data   
      !$omp end target data
    end if
    call system_clock(finish, clock_rate)
    elapsed(iloop) = finish - start
  end do
  !
  call system_clock(finish_all, clock_rate)

  avg  = sum( elapsed(merge(2,1,nloops>1):) ) / real(merge(nloops-1,nloops,nloops>1))
  mint = minval(elapsed) 

  ! What to print? 
  !   ncol, nlay, ngpt; are clouds used, are aerosols used; time per column, total, min; 
  print *, " ncol   nlay   ngpt  clouds aerosols time_per_col_ms nloops time_total_s time_min_s"
  write(output_unit, '(3(i6, 1x), 6x, 2(i1, 8x), 1x, f7.3, 1x, i6, 2x, 2(4x,f7.3))') & 
    ncol, nlay, ngpt, merge(1,0,do_clouds), merge(1,0,do_aerosols),  & 
    avg/(real(ncol) * (1.0e-3*clock_rate)),  nloops,  sum(elapsed) / real(clock_rate),  mint / real(clock_rate)

  call write_fluxes

  ! 
  ! Memory for bounday conditions on the GPU was allocated with unstructured data dataments 
  !   (acc enter data). Deallocate it expliicity 
  !
  if(is_lw) then
    !$acc        exit data delete(     t_sfc, emis_sfc)
    !$omp target exit data map(release:t_sfc, emis_sfc)
  else
    !$acc        exit data delete(     sfc_alb_dir, sfc_alb_dif, mu0)
    !$omp target exit data map(release:sfc_alb_dir, sfc_alb_dif, mu0)
  end if
  
  !
  ! Clouds and aerosols also used enter data  
  !
  if(do_clouds) then
    !$acc        exit data delete(     cloud_mask, lwp, iwp, rel, rei)
    !$omp target exit data map(release:cloud_mask, lwp, iwp, rel, rei)
    select type(clouds)
      class is (ty_optical_props_1scl)
        !$acc        exit data delete     (clouds%tau, clouds)
        !$omp target exit data map(release:clouds%tau)
      class is (ty_optical_props_2str)
        !$acc        exit data delete     (clouds%tau, clouds%ssa, clouds%g, clouds)
        !$omp target exit data map(release:clouds%tau, clouds%ssa, clouds%g)
    end select
    ! 
    ! Explicit finalization of cloud optical properties - not really necessary since memory 
    !   will be freed when the program ends, but useful for testing 
    !
    call clouds%finalize
  end if
  if(do_aerosols) then
    !$acc        exit data delete(     aero_type, aero_size, aero_mass, relhum)
    !$omp target exit data map(release:aero_type, aero_size, aero_mass, relhum)
    select type(aerosols)
      class is (ty_optical_props_1scl)
        !$acc        exit data delete     (aerosols%tau, aerosols)
        !$omp target exit data map(release:aerosols%tau)
      class is (ty_optical_props_2str)
        !$acc        exit data delete     (aerosols%tau, aerosols%ssa, aerosols%g, aerosols)
        !$omp target exit data map(release:aerosols%tau, aerosols%ssa, aerosols%g)
    end select
    ! 
    ! Explicit finalization of aerosol optical properties - not really necessary since memory 
    !   will be freed when the program ends, but useful for testing 
    !
    call aerosols%finalize
  end if
  !
  ! k-distribution
  !
  call k_dist%finalize
  
  if(.not. is_lw) then
    !$acc        exit data delete(     flux_dir)
    !$omp target exit data map(release:flux_dir)
  end if

  ! fluxes - but not flux_dir, which used enter data 
  !$acc end        data 
  !$omp end target data
  ! p_lay etc
  !$acc end        data 
  !$omp end target data
contains
  ! ----------------------------------------------------------------------------------
  subroutine compute_profiles(SST, ncol, nlay, p_lay, t_lay, p_lev, t_lev, q_lay, o3)
    !
    ! Construct profiles of pressure, temperature, humidity, and ozone 
    !   more or less following the RCEMIP protocol for a surface temperature of 300K
    !   more or less follows a Python implementation by Chiel van Heerwardeen
    ! Extensions for future - variable SST and T profile, variable RH, lapse rate in stratosphere 
    !   will all access absorption coefficient data more realistically 
    !
    real(wp),                          intent(in ) :: SST 
    integer,                           intent(in ) :: ncol, nlay
    real(wp), dimension(ncol, nlay  ), intent(out) :: p_lay, t_lay, q_lay, o3
    real(wp), dimension(ncol, nlay+1), intent(out) :: p_lev, t_lev

    real(wp) :: z_lay(nlay), z_lev(nlay+1)
    real(wp) :: z, q, T, p
    real(wp) :: Tv, Tv0, p_hpa
    integer  :: icol, ilay, i

    real(wp), parameter :: z_trop = 15000._wp, z_top = 70.e3_wp
    ! Ozone profile - maybe only a single profile? 
    real(wp), parameter :: g1 = 3.6478_wp, g2 = 0.83209_wp, g3 = 11.3515_wp, o3_min = 1e-13_wp 
    ! According to CvH RRTMGP in Single Precision will fail with lower ozone concentrations

    real(wp), parameter :: g = 9.79764, Rd = 287.04, p0 = 101480. ! Surface pressure 
    real(wp), parameter :: z_q1 = 4.0e3, z_q2 = 7.5e3,  q_t = 1.e-8
    real(wp), parameter :: gamma = 6.7e-3
    
    real(wp), parameter :: q_0 = 0.01864 ! for 300 K SST.
    ! -------------------
    Tv0 = (1. + 0.608*q_0) * SST
    !
    ! Split resolution above and below RCE tropopause (15 km or about 125 hPa)
    !
    z_lev(:) = [0._wp,  2._wp*           z_trop /nlay * [(i, i=1, nlay/2)],  & 
               z_trop + 2._wp * (z_top - z_trop)/nlay * [(i, i=1, nlay/2)]]
    z_lay(:) = 0.5_wp * (z_lev(1:nlay)  + z_lev(2:nlay+1))
    
    !$acc        data copyin(z_lev, z_lay) 
    !$omp target data map(to:z_lev, z_lay)

    !
    ! The two loops are the same, except applied to layers and levels 
    !   but nvfortran doesn't seems to support elemental procedures in OpenACC loops
    !
    !$acc                         parallel loop    collapse(2) 
    !$omp target teams distribute parallel do simd collapse(2) 
    do ilay = 1, nlay 
      do icol = 1, ncol 
        z = z_lay(ilay) 
        if (z > z_trop) then 
          q = q_t
          T = SST - gamma*z_trop/(1. + 0.608*q_0)
          Tv  = (1. + 0.608*q  ) *   T
          p = p0 * (Tv/Tv0)**(g/(Rd*gamma)) * exp( -((g*(z-z_trop))/(Rd*Tv)) )
        else 
          q = q_0 * exp(-z/z_q1) * exp(-(z/z_q2)**2)
          T = SST - gamma*z / (1. + 0.608*q)    
          Tv  = (1. + 0.608*q  ) *   T
          p = p0 * (Tv/Tv0)**(g/(Rd*gamma))
        end if 
        p_lay(icol,ilay) = p 
        t_lay(icol,ilay) = T
        q_lay(icol,ilay) = q
        p_hpa = p_lay(icol,ilay) / 100._wp
        o3(icol, ilay) = max(o3_min, & 
                             g1 * p_hpa**g2 * exp(-p_hpa/g3) * 1.e-6_wp)
      end do
    end do 

    !$acc                         parallel loop    collapse(2) 
    !$omp target teams distribute parallel do simd collapse(2) 
    do ilay = 1, nlay+1
      do icol = 1, ncol 
        z = z_lev(ilay) 
        if (z > z_trop) then 
          q = q_t
          T = SST - gamma*z_trop/(1. + 0.608*q_0)
          Tv  = (1. + 0.608*q  ) *   T
          p = p0 * (Tv/Tv0)**(g/(Rd*gamma)) * exp( -((g*(z-z_trop))/(Rd*Tv)) )
        else 
          q = q_0 * exp(-z/z_q1) * exp(-(z/z_q2)**2)
          T = SST - gamma*z / (1. + 0.608*q)    
          Tv  = (1. + 0.608*q  ) *   T
          p = p0 * (Tv/Tv0)**(g/(Rd*gamma))
        end if 
        p_lev(icol,ilay) = p 
        t_lev(icol,ilay) = T
      end do 
    end do 
    !$acc end        data
    !$omp end target data
  end subroutine compute_profiles
  ! ----------------------------------------------------------------------------------
  subroutine stop_on_err(error_msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_msg

    if(error_msg /= "") then
      write (error_unit,*) trim(error_msg)
      write (error_unit,*) "rrtmgp_allsky stopping"
      error stop 1
    end if
  end subroutine stop_on_err
  ! --------------------------------------------------------------------------------------
  !
  subroutine compute_clouds 
    real(wp) :: rel_val, rei_val
    ! 
    ! Variable and memory allocation 
    !
    if(is_sw) then
      allocate(ty_optical_props_2str::clouds)
    else
      allocate(ty_optical_props_1scl::clouds)
    end if
    ! Clouds optical props are defined by band
    call stop_on_err(clouds%init(k_dist%get_band_lims_wavenumber()))

    select type(clouds)
      class is (ty_optical_props_1scl)
        call stop_on_err(clouds%alloc_1scl(ncol, nlay))
        !$acc enter data copyin(clouds) create(clouds%tau)
        !$omp target enter data map(alloc:clouds%tau)
      class is (ty_optical_props_2str)
        call stop_on_err(clouds%alloc_2str(ncol, nlay))
        !$acc enter data copyin(clouds) create(clouds%tau, clouds%ssa, clouds%g)
        !$omp target enter data map(alloc:clouds%tau, clouds%ssa, clouds%g)
      class default
        call stop_on_err("rte_rrtmgp_allsky: Don't recognize the kind of optical properties ")
    end select
    !
    ! Cloud physical properties 
    !
    allocate(lwp(ncol,nlay), iwp(ncol,nlay), &
             rel(ncol,nlay), rei(ncol,nlay), cloud_mask(ncol,nlay))
    !$acc enter        data create(   cloud_mask, lwp, iwp, rel, rei)
    !$omp target enter data map(alloc:cloud_mask, lwp, iwp, rel, rei)

    ! Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
    !   and not very close to the ground (< 900 hPa), and
    !   put them in 2/3 of the columns since that's roughly the
    !   total cloudiness of earth
    rel_val = 0.5 * (cloud_optics%get_min_radius_liq() + cloud_optics%get_max_radius_liq())
    rei_val = 0.5 * (cloud_optics%get_min_radius_ice() + cloud_optics%get_max_radius_ice())
    !$acc                         parallel loop    collapse(2) copyin(t_lay) copyout( lwp, iwp, rel, rei)
    !$omp target teams distribute parallel do simd collapse(2) map(to:t_lay) map(from:lwp, iwp, rel, rei)
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
    !$omp target exit data map(release:cloud_mask)
   
  end subroutine compute_clouds
  !
  ! --------------------------------------------------------------------------------------
  !
  subroutine compute_aerosols
    real(wp), dimension(ncol,nlay) :: vmr_h2o ! h2o vmr
    logical :: is_sulfate, is_dust, is_even_column 
    ! 
    ! Variable and memory allocation 
    !
    if(is_sw) then
      allocate(ty_optical_props_2str::aerosols)
    else
      allocate(ty_optical_props_1scl::aerosols)
    end if
    call stop_on_err(aerosols%init(k_dist%get_band_lims_wavenumber()))
    select type(aerosols)
      class is (ty_optical_props_1scl)
        call stop_on_err(aerosols%alloc_1scl(ncol, nlay))
        !$acc        enter data copyin(aerosols) create(aerosols%tau)
        !$omp target enter data map              (alloc:aerosols%tau)
      class is (ty_optical_props_2str)
        call stop_on_err(aerosols%alloc_2str(ncol, nlay))
        !$acc        enter data copyin(aerosols) create(aerosols%tau, aerosols%ssa, aerosols%g)
        !$omp target enter data map              (alloc:aerosols%tau, aerosols%ssa, aerosols%g)
      class default
        call stop_on_err("rte_rrtmgp_allsky: Don't recognize the kind of optical properties ")
    end select
    !
    ! Derive relative humidity from profile
    !   Keep vmr_h2o on the GPU
    !
    !$acc        data create(   vmr_h2o)
    !$omp target data map(alloc:vmr_h2o)
    call stop_on_err(gas_concs%get_vmr("h2o",vmr_h2o))
    !
    ! Aerosol properties 
    ! 
    allocate(aero_type(ncol,nlay), aero_size(ncol,nlay), &
             aero_mass(ncol,nlay), relhum   (ncol,nlay))
    !$acc        enter data create(   aero_type, aero_size, aero_mass, relhum)
    !$omp target enter data map(alloc:aero_type, aero_size, aero_mass, relhum)
    call get_relhum(ncol, nlay, p_lay, t_lay, vmr_h2o, relhum)
    !$acc end data
    !$omp end target data

    ! Restrict sulfate aerosols to lower stratosphere (> 50 hPa = 50*100 Pa; < 100 hPa = 100*100 Pa)
    !   and dust aerosols to the lower troposphere (> 700 hPa; < 900 hPa), and
    !   put them in 1/2 of the columns
    !
    !
    !$acc                         parallel loop    collapse(2) copyin(p_lay) 
    !$omp target teams distribute parallel do simd collapse(2) map(to:p_lay) 
    do ilay=1,nlay
      do icol=1,ncol
        is_sulfate = (p_lay(icol,ilay) >  50._wp * 100._wp .and. & 
                      p_lay(icol,ilay) < 100._wp * 100._wp)
        is_dust    = (p_lay(icol,ilay) > 700._wp * 100._wp .and. & 
                      p_lay(icol,ilay) < 900._wp * 100._wp)
        is_even_column = mod(icol, 2) /= 0
        if      (is_even_column .and. is_sulfate) then 
          aero_type(icol,ilay) = merra_aero_sulf
          aero_size(icol,ilay) = 0.2_wp
          aero_mass(icol,ilay) = 1.e-6_wp
        else if(is_even_column .and. is_dust) then 
          ! Dust aerosol
          aero_type(icol,ilay) = merra_aero_dust
          aero_size(icol,ilay) = 0.5_wp
          aero_mass(icol,ilay) = 3.e-5_wp
        else
          aero_type(icol,ilay) = 0
          aero_size(icol,ilay) = 0._wp
          aero_mass(icol,ilay) = 0._wp
        end if
      end do
    end do

    end subroutine compute_aerosols
  ! --------------------------------------------------------------------------------------
  !
  ! Calculate layer relative humidity for aerosol optics calculations
  !
  subroutine get_relhum(ncol, nlay, p_lay, t_lay, vmr_h2o, relhum)
    use mo_rte_kind,           only: wp
    use mo_gas_optics_constants,   only: m_h2o, m_dry

    integer,  intent(in) :: ncol, nlay
    real(wp), intent(in) :: p_lay(ncol,nlay)    ! layer pressure (Pa)
    real(wp), intent(in) :: t_lay(ncol,nlay)    ! layer temperature (K)
    real(wp), intent(in) :: vmr_h2o(ncol,nlay)  ! water volume mixing ratio

    real(wp), intent(inout) :: relhum(ncol,nlay) ! relative humidity (fraction, 0-1)

    ! Local variables 
    integer :: i, k

    real(wp) :: mmr_h2o          ! water mass mixing ratio
    real(wp) :: q_lay            ! water specific humidity
    real(wp) :: q_lay_min, q_tmp, es_tmp
    real(wp) :: mwd, t_ref, rh

    ! Set constants
    mwd       = m_h2o/m_dry      ! ratio of water to dry air molecular weights
    t_ref     = 273.16_wp        ! reference temperature (K)
    q_lay_min = 1.e-7_wp         ! minimum water mass mixing ratio
    ! -------------------

    ! Derive layer virtual temperature
    !$acc                         parallel loop    collapse(2) copyin(p_lay, vmr_h2o, t_lay) copyout( relhum)
    !$omp target teams distribute parallel do simd collapse(2) map(to:p_lay, vmr_h2o, t_lay) map(from:relhum) 
    do i = 1, ncol 
       do k = 1, nlay
          ! Convert h2o vmr to mmr
          mmr_h2o = vmr_h2o(i,k) * mwd
          q_lay = mmr_h2o / (1 + mmr_h2o)
          q_tmp = max(q_lay_min, q_lay)
          es_tmp = exp( (17.67_wp * (t_lay(i,k)-t_ref)) / (t_lay(i,k)-29.65_wp) )
          rh = (0.263_wp * p_lay(i,k) * q_tmp) / es_tmp
          ! Convert rh from percent to fraction
          relhum(i,k) = 0.01_wp * rh
       enddo
    enddo
  end subroutine get_relhum
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_fluxes 
    use netcdf
    use mo_simple_netcdf, only: write_field
    integer :: ncid, i, col_dim, lay_dim, lev_dim, varid
    real(wp) :: vmr(ncol, nlay)
    character(len=3) :: flux_prefix
    !
    ! Write fluxes - make this optional? 
    !

    !
    ! Define dimensions 
    !
    if(nf90_create(trim(output_file),  NF90_CLOBBER, ncid) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't create file " // trim(output_file))

    if(nf90_def_dim(ncid, "col", ncol, col_dim) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't define col dimension")
    if(nf90_def_dim(ncid, "lay", nlay, lay_dim) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't define lay dimension")
    if(nf90_def_dim(ncid, "lev", nlay+1, lev_dim) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't define lev dimension")
    
    !
    ! Define variables 
    !
    ! State 
    !
    call create_var("p_lev", ncid, [col_dim, lev_dim])
    call create_var("t_lev", ncid, [col_dim, lev_dim])
    call create_var("p_lay", ncid, [col_dim, lay_dim])
    call create_var("t_lay", ncid, [col_dim, lay_dim])
    call create_var("h2o",   ncid, [col_dim, lay_dim])
    call create_var("o3",    ncid, [col_dim, lay_dim])

    ! All the gases except h2o, o3 - write as attributes? Or not bother? 

    if(do_clouds) then 
      call create_var("lwp", ncid, [col_dim, lay_dim])
      call create_var("iwp", ncid, [col_dim, lay_dim])
      call create_var("rel", ncid, [col_dim, lay_dim])
      call create_var("rei", ncid, [col_dim, lay_dim])
    end if 
    if(do_aerosols) then 
      if(nf90_def_var(ncid, "aero_type", NF90_SHORT, [col_dim, lay_dim], varid) /= NF90_NOERR) &
        call stop_on_err("create_var: can't define variable aero_type")
      call create_var("aero_size", ncid, [col_dim, lay_dim])
      call create_var("aero_mass", ncid, [col_dim, lay_dim])
   end if 
    !
    ! Fluxes - definitions 
    !
    if(is_sw) then 
      flux_prefix = "sw_"
      call create_var(flux_prefix // "flux_dir", ncid, [col_dim, lev_dim])
    else
      flux_prefix = "lw_"  
    end if 
    call create_var(flux_prefix // "flux_up", ncid, [col_dim, lev_dim])
    call create_var(flux_prefix //"flux_dn", ncid, [col_dim, lev_dim])
    if(nf90_enddef(ncid) /= NF90_NOERR) &
      call stop_on_err("rrtmgp_allsky: can't end file definition??")

    !
    ! Write variables
    !
    ! State - writing 
    !$acc        update host(p_lev, t_lev, p_lay, t_lay)
    !$omp target update from(p_lev, t_lev, p_lay, t_lay)
    call stop_on_err(write_field(ncid, "p_lev",  p_lev))
    call stop_on_err(write_field(ncid, "t_lev",  t_lev))
    call stop_on_err(write_field(ncid, "p_lay",  p_lay))
    call stop_on_err(write_field(ncid, "t_lay",  t_lay))
    ! Array vmr is on the host, not the device, but is copied-out
    call stop_on_err(gas_concs%get_vmr("h2o", vmr))
    call stop_on_err(write_field(ncid, "h2o",    vmr))
    call stop_on_err(gas_concs%get_vmr("o3",  vmr))
    call stop_on_err(write_field(ncid, "o3",     vmr))

    if(do_clouds) then 
      !$acc        update host(lwp, iwp, rel, rei)
      !$omp target update from(lwp, iwp, rel, rei)
      call stop_on_err(write_field(ncid, "lwp",  lwp))
      call stop_on_err(write_field(ncid, "iwp",  iwp))
      call stop_on_err(write_field(ncid, "rel",  rel))
      call stop_on_err(write_field(ncid, "rei",  rei))
    end if 

    if(do_aerosols) then 
      !$acc        update host(aero_size, aero_mass, aero_type)
      !$omp target update from(aero_size, aero_mass, aero_type)
      call stop_on_err(write_field(ncid, "aero_size",  aero_size))
      call stop_on_err(write_field(ncid, "aero_mass",  aero_mass))
      call stop_on_err(write_field(ncid, "aero_type",  aero_type))
    end if 

    ! Fluxes - writing 
    !$acc        update host(flux_up, flux_dn)
    !$omp target update from(flux_up, flux_dn)
    call stop_on_err(write_field(ncid, flux_prefix // "flux_up",  flux_up))
    call stop_on_err(write_field(ncid, flux_prefix // "flux_dn",  flux_dn))
    if(.not. is_lw) then 
      !$acc        update host(flux_dir)
      !$omp target update from(flux_dir)
      call stop_on_err(write_field(ncid, flux_prefix // "flux_dir",  flux_dir))
    end if 

    ! Close netCDF 
    if(nf90_close(ncid) /= NF90_NOERR) call stop_on_err("rrtmgp_allsky: error closing file??")
  end subroutine write_fluxes
  ! ---------------------------------------------------------
  subroutine create_var(name, ncid, dim_ids)
    use netcdf
    character(len=*),      intent(in) :: name 
    integer,               intent(in) :: ncid
    integer, dimension(:), intent(in) :: dim_ids 

    integer :: varid

    if(nf90_def_var(ncid, trim(name), NF90_DOUBLE, dim_ids, varid) /= NF90_NOERR) &
      call stop_on_err("create_var: can't define " // trim(name) // " variable")
  end subroutine create_var
  ! ---------------------------------------------------------
end program rte_rrtmgp_allsky
