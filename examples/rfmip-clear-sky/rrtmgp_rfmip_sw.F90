! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Example program to demonstrate the calculation of shortwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_sw [block_size input_file  coefficient_file upflux_file downflux_file]
!   All arguments are optional but need to be specified in order.
! Or as ssm_rfmip_sw [block size] to test the simple spectral model
!
! -------------------------------------------------------------------------------------------------
program rrtmgp_rfmip_sw
  ! --------------------------------------------------
  !
  ! Modules for working with rte and rrtmgp
  !
  ! Working precision for real variables
  !
  use mo_rte_kind,           only: wp
  !
  ! Array utilities
  !
  use mo_rte_util_array,     only: zero_array
  !
  ! Optical properties of the atmosphere as array of values
  !   In the longwave we include only absorption optical depth (_1scl)
  !   Shortwave calculations use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  !
  use mo_optical_props,      only: ty_optical_props_2str
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly
  !
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! RTE shortwave driver
  !
  use mo_rte_sw,             only: rte_sw
  !
  ! RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  !   Here we're just reporting broadband fluxes
  !
  use mo_fluxes,             only: ty_fluxes_broadband
  ! --------------------------------------------------
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !    This example can use either a k-distribution from RRTMGP or a simple spectral model
  !    The optics that gets used is chosen at run time from the program name
  !
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_optics_ssm,         only: ty_optics_ssm
  !
  ! modules for reading and writing files
  !
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, &
                                   read_and_block_sw_bc, determine_gas_names
  use mo_testing_utils,      only: stop_on_err
  !
  ! RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
  !
  use mo_optics_utils_rrtmgp,only: load_gas_optics
  implicit none
  ! --------------------------------------------------
  !
  ! Local variables
  !
  character(len=512) :: invoked_name
  character(len=512) :: rfmip_file = 'multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc'
  character(len=512) :: kdist_file = 'coefficients_sw.nc'
  character(len=132) :: flxdn_file, flxup_file
  integer            :: nargs, ncol, nlay, nbnd, ngpt, nexp, nblocks, block_size, forcing_index
  integer            :: b, icol, ibnd, igpt
  character(len=4)   :: block_size_char, forcing_index_char = '1'
  logical            :: do_rrtmgp, do_ssm

  character(len=32 ), &
            dimension(:),             allocatable :: kdist_gas_names, rfmip_gas_games
  real(wp), dimension(:,:,:),         allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
  real(wp), dimension(:,:  ),         allocatable :: surface_albedo, total_solar_irradiance, solar_zenith_angle
                                                     ! block_size, nblocks
  real(wp), dimension(:,:  ),         allocatable :: sfc_alb_spec ! nbnd, block_size; spectrally-resolved surface albedo
  !
  ! Classes used by rte+rrtmgp
  !
  ! mo_optics_utils_rrtmgp declares a variable
  type(ty_optical_props_2str)                    :: optical_props
  type(ty_fluxes_broadband)                      :: fluxes

  real(wp), dimension(:,:), allocatable          :: toa_flux ! block_size, ngpt
  real(wp), dimension(:  ), allocatable          :: def_tsi, mu0    ! block_size
  logical , dimension(:,:), allocatable          :: usecol ! block_size, nblocks
  !
  ! Optics is determined at run time
  !
  class(ty_gas_optics), allocatable :: gas_optics
  !
  ! ty_gas_concentration holds multiple columns; we make an array of these objects to
  !   leverage what we know about the input file
  !
  type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array
  real(wp), parameter :: deg_to_rad = acos(-1._wp)/180._wp
  ! -------------------------------------------------------------------------------------------------
  !
  ! Code starts
  !
  ! Determine which gas optics to use based on the name by which the program is evoked
  !   (possibly fragile)
  ! Based on the possibilities: rrtmgp_rfmip_lw, ssm_rfmip_lw
  call get_command_argument(0, invoked_name)
  do_rrtmgp = (invoked_name(len_trim(invoked_name)-14:len_trim(invoked_name)-8) == "rrtmgp_")
  do_ssm    = (invoked_name(len_trim(invoked_name)-11:len_trim(invoked_name)-8) == "ssm_")
  if (.not. (do_rrtmgp .or. do_ssm)) call stop_on_err("Don't recogize which optics to use")

  nargs = command_argument_count()
  if(nargs >= 2) call get_command_argument(2, rfmip_file)
  call read_size(rfmip_file, ncol, nlay, nexp)
  if(nargs >= 1) then
    call get_command_argument(1, block_size_char)
    read(block_size_char, '(i4)') block_size
  else
    block_size = ncol
  end if

  if(do_rrtmgp) then
    allocate(ty_gas_optics_rrtmgp::gas_optics)
    print *, "Usage: rrtmgp_rfmip_sw [block_size] [rfmip_file] [k-distribution_file] [forcing_index (1,2,3)]"
    !
    ! Read arguments out of order
    !
    if(nargs >= 3) call get_command_argument(3, kdist_file)
    if(nargs >= 4) then
      call get_command_argument(4, forcing_index_char)
    end if
    read(forcing_index_char, '(i4)') forcing_index
    if(forcing_index < 1 .or. forcing_index > 3) &
      stop "Forcing index is invalid (must be 1,2 or 3)"

    flxdn_file = 'rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f' // trim(forcing_index_char) // '_gn.nc'
    flxup_file = 'rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f' // trim(forcing_index_char) // '_gn.nc'

    !
    ! Identify the set of gases used in the calculation based on the forcing index
    !   A gas might have a different name in the k-distribution than in the files
    !   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
    !
    call determine_gas_names(rfmip_file, kdist_file, forcing_index, kdist_gas_names, rfmip_gas_games)
    print *, "Calculation uses RFMIP gases: ", (trim(rfmip_gas_games(b)) // " ", b = 1, size(rfmip_gas_games))
  else if (do_ssm) then
    allocate(ty_optics_ssm::gas_optics)
    print *, "Usage: ssm_rfmip_sw [block_size]"
    flxdn_file = 'rsd_ssm_rfmip-rad-irf.nc'
    flxup_file = 'rsu_ssm_rfmip-rad-irf.nc'
    !
    ! These variables are needed for the fragile RFMIP IO
    !
    kdist_gas_names = ["co2"]
    rfmip_gas_games = ["carbon_dioxide"]
  end if
  !
  ! How big is the problem? Does it fit into blocks of the size we've specified?
  !
  if(mod(ncol*nexp, block_size) /= 0 ) call stop_on_err("rrtmgp_rfmip_sw: number of columns doesn't fit evenly into blocks.")
  nblocks = (ncol*nexp)/block_size
  print *, "Doing ",  nblocks, "blocks of size ", block_size

  ! --------------------------------------------------
  !
  ! Prepare data for use in rte+rrtmgp
  !
  !
  ! Allocation on assignment within reading routines
  !
  call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, t_lay, t_lev)

  !
  ! Read the gas concentrations and surface properties
  !
  call read_and_block_gases_ty(rfmip_file, block_size, kdist_gas_names, rfmip_gas_games, gas_conc_array)
  call read_and_block_sw_bc(rfmip_file, block_size, surface_albedo, total_solar_irradiance, solar_zenith_angle)

  select type (gas_optics)
    type is (ty_gas_optics_rrtmgp)
      !
      ! Read k-distribution information. load_gas_optics() reads data from netCDF and calls
      !   gas_optics%init(); users might want to use their own reading methods
      !
      call load_gas_optics(gas_optics, trim(kdist_file), gas_conc_array(1))
      if(.not. gas_optics%source_is_external()) &
        stop "rrtmgp_rfmip_sw: k-distribution file isn't SW"
      ngpt = gas_optics%get_ngpt()

      !
      ! RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
      !   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
      !   This introduces an error but shows input sanitizing.
      !
      !
      ! Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
      !
      if(p_lay(1, 1, 1) < p_lay(1, nlay, 1)) then
        p_lev(:,1,:) = gas_optics%get_press_min() + epsilon(gas_optics%get_press_min())
      else
        p_lev(:,nlay+1,:) &
                     = gas_optics%get_press_min() + epsilon(gas_optics%get_press_min())
      end if
    type is (ty_optics_ssm)
      call stop_on_err(gas_optics%configure(do_sw = .true.))
  end select
  nbnd = gas_optics%get_nband()


  allocate(toa_flux(block_size, gas_optics%get_ngpt()), &
           def_tsi(block_size), usecol(block_size,nblocks))
  !
  ! RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
  !   nighttime columns with a default solar zenith angle. We'll mask these out later, of
  !   course, but this gives us more work and so a better measure of timing.
  !
  do b = 1, nblocks
    usecol(1:block_size,b)  = solar_zenith_angle(1:block_size,b) < 90._wp - 2._wp * spacing(90._wp)
  end do

  !
  ! Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
  !   gas optical properties, and source functions. The %alloc() routines carry along
  !   the spectral discretization from the k-distribution.
  !
  allocate(flux_up(block_size, nlay+1, nblocks), &
           flux_dn(block_size, nlay+1, nblocks))
  allocate(mu0(block_size), sfc_alb_spec(nbnd,block_size))
  call stop_on_err(optical_props%alloc_2str(block_size, nlay, gas_optics))
  !$acc enter data create(optical_props, optical_props%tau, optical_props%ssa, optical_props%g)
  !$omp target enter data map(alloc:optical_props%tau, optical_props%ssa, optical_props%g)
  !$acc enter data create (toa_flux, def_tsi)
  !$omp target enter data map(alloc:toa_flux, def_tsi)
  !$acc enter data create (sfc_alb_spec, mu0)
  !$omp target enter data map(alloc:sfc_alb_spec, mu0)
  ! --------------------------------------------------
  !
  ! Loop over blocks
  !
  do b = 1, nblocks
    fluxes%flux_up => flux_up(:,:,b)
    fluxes%flux_dn => flux_dn(:,:,b)
    !
    ! Compute the optical properties of the atmosphere and the Planck source functions
    !    from pressures, temperatures, and gas concentrations...
    !
    call stop_on_err(gas_optics%gas_optics(p_lay(:,:,b), &
                                           p_lev(:,:,b),       &
                                           t_lay(:,:,b),       &
                                           gas_conc_array(b),  &
                                           optical_props,      &
                                           toa_flux))
    ! Boundary conditions
    !   (This is partly to show how to keep work on GPUs using OpenACC in a host application)
    ! What's the total solar irradiance assumed by RRTMGP?
    !
#if defined(_OPENACC) || defined(_OPENMP)
    call zero_array(block_size, def_tsi)
    !$acc parallel loop collapse(2) copy(def_tsi) copyin(toa_flux)
    !$omp target teams distribute parallel do simd collapse(2) map(tofrom:def_tsi) map(to:toa_flux)
    do igpt = 1, ngpt
      do icol = 1, block_size
        !$acc atomic update
        !$omp atomic update
        def_tsi(icol) = def_tsi(icol) + toa_flux(icol, igpt)
      end do
    end do
#else
    !
    ! More compactly...
    !
    def_tsi(1:block_size) = sum(toa_flux, dim=2)
#endif
    !
    ! Normalize incoming solar flux to match RFMIP specification
    !
    !$acc parallel loop collapse(2) copyin(total_solar_irradiance, def_tsi) copy(toa_flux)
    !$omp target teams distribute parallel do simd collapse(2) map(to:total_solar_irradiance, def_tsi) map(tofrom:toa_flux)
    do igpt = 1, ngpt
      do icol = 1, block_size
        toa_flux(icol,igpt) = toa_flux(icol,igpt) * total_solar_irradiance(icol,b)/def_tsi(icol)
      end do
    end do
    !
    ! Expand the spectrally-constant surface albedo to a per-band albedo for each column
    !
    !$acc parallel loop collapse(2) copyin(surface_albedo)
    !$omp target teams distribute parallel do simd collapse(2) map(to:surface_albedo)
    do icol = 1, block_size
      do ibnd = 1, nbnd
        sfc_alb_spec(ibnd,icol) = surface_albedo(icol,b)
      end do
    end do
    !
    ! Cosine of the solar zenith angle
    !
    !$acc parallel loop copyin(solar_zenith_angle, usecol)
    !$omp target teams distribute parallel do simd map(to:solar_zenith_angle, usecol)
    do icol = 1, block_size
      mu0(icol) = merge(cos(solar_zenith_angle(icol,b)*deg_to_rad), 1._wp, usecol(icol,b))
    end do

    !
    ! ... and compute the spectrally-resolved fluxes, providing reduced values
    !    via ty_fluxes_broadband
    !
    call stop_on_err(rte_sw(optical_props,   &
                            mu0,             &
                            toa_flux,        &
                            sfc_alb_spec,    &
                            sfc_alb_spec,    &
                            fluxes))
    !
    ! Zero out fluxes for which the original solar zenith angle is > 90 degrees.
    !
    do icol = 1, block_size
      if(.not. usecol(icol,b)) then
        flux_up(icol,:,b)  = 0._wp
        flux_dn(icol,:,b)  = 0._wp
      end if
    end do
  end do
  !
  ! End timers
  !
  !$acc exit data delete(optical_props%tau, optical_props%ssa, optical_props%g, optical_props)
  !$omp target exit data map(release:optical_props%tau, optical_props%ssa, optical_props%g)
  !$acc exit data delete(sfc_alb_spec, mu0)
  !$omp target exit data map(release:sfc_alb_spec, mu0)
  !$acc exit data delete(toa_flux, def_tsi)
  !$omp target exit data map(release:toa_flux, def_tsi)
  ! --------------------------------------------------
  call unblock_and_write(trim(flxup_file), 'rsu', flux_up)
  call unblock_and_write(trim(flxdn_file), 'rsd', flux_dn)
end program rrtmgp_rfmip_sw
