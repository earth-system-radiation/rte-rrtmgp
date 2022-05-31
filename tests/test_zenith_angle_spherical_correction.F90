subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "test_solar_zenith_angle stopping"
    error stop 1
  end if
end subroutine stop_on_err
! ------------------------------------------------------------------------------
program test_solar_zenith_angle
  use mo_rte_kind,           only: wp, wl
  use mo_zenith_angle_spherical_correction, &
                             only: zenith_angle_with_height

  use mo_rcemip_profiles,    only: make_rcemip_profiles
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_optical_props,      only: ty_optical_props_2str
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_heating_rates,      only: compute_heating_rate

  integer            :: i
  integer, parameter :: nlay = 80, ncol = 4
  real(wp), dimension(ncol        ) :: refAlt, refMu
  real(wp), dimension(ncol, 0:nlay) ::    alt,    mu

  type(ty_gas_concs)          :: gas_concs
  type(ty_gas_optics_rrtmgp)  :: gas_optics
  type(ty_optical_props_2str) :: atmos
  type(ty_fluxes_broadband)   :: fluxes
  real(wp)                    :: p_lay(ncol, nlay), t_lay(ncol, nlay), &
                                 z_lay(ncol, nlay), mu0  (ncol, nlay), &
                                 heating_rate(ncol, nlay),             &
                                 p_lev(ncol, nlay+1)
  real(wp), dimension(ncol, nlay+1), target &
                              :: flux_up, flux_dn, flux_dir
  logical(wl)                 :: top_at_1
  character(len=128)          :: k_dist_file = "../rrtmgp/data/rrtmgp-data-sw-g112-210809.nc"
  real(wp), dimension(:,:), allocatable &
                              :: toa_flux, sfc_alb_dir, sfc_alb_dif
  !
  ! Geometric calculation
  !
  print *, "Geometry "
  refAlt(:) = nlay * 1000._wp   ! Reference altitude is TOA
  refMu (:) = .02
  alt =  spread([(i, i = 0, nlay)] * 1000._wp, dim = 1, ncopies=ncol)
  call stop_on_err(zenith_angle_with_height(refAlt, refMu, alt, mu))
  print *, "alt     mu0"
  do i = nlay, 0, -1
    print *, int(alt(1,i)/1000._wp), mu(1,i)
  end do

  !
  ! fluxes
  !
  call stop_on_err(make_rcemip_profiles(p_lay(1,:), p_lev(1,:), t_lay(1,:), z_lay(1,:), gas_concs))
  p_lay = spread(p_lay(1,:), dim=1, ncopies=ncol)
  t_lay = spread(t_lay(1,:), dim=1, ncopies=ncol)
  z_lay = spread(z_lay(1,:), dim=1, ncopies=ncol)
  p_lev = spread(p_lev(1,:), dim=1, ncopies=ncol)
  top_at_1 = p_lay(1, 1) < p_lay(1,nlay)
  print *, "Top is at 1? ", top_at_1

  call load_and_init(gas_optics, k_dist_file, gas_concs)
  print *, "temperature limits (K):", gas_optics%get_temp_min(),  gas_optics%get_temp_max()

  call stop_on_err(atmos%alloc_2str(ncol, nlay, gas_optics))
  if(gas_optics%source_is_internal()) call stop_on_err("k-distribution file is for the longwave")
  nbnd = gas_optics%get_nband()
  ngpt = gas_optics%get_ngpt()

  allocate(toa_flux(ncol, ngpt))
  call stop_on_err(gas_optics%gas_optics(p_lay, p_lev, &
                                         t_lay,        &
                                         gas_concs,    &
                                         atmos,        &
                                         toa_flux))

  allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol))
  sfc_alb_dir = 0._wp
  sfc_alb_dif = 0._wp
  fluxes%flux_up     => flux_up (:,:)
  fluxes%flux_dn     => flux_dn (:,:)
  fluxes%flux_dn_dir => flux_dir(:,:)

  call stop_on_err(zenith_angle_with_height(z_lay(:,1), [0.01_wp, 0.02_wp, 0.03_wp, 0.04_wp], z_lay, mu0))
  call stop_on_err(rte_sw(atmos, top_at_1, &
                          mu0,   toa_flux, &
                          sfc_alb_dir, sfc_alb_dif, &
                          fluxes))
  call stop_on_err(compute_heating_rate(flux_up, flux_dn, p_lev, heating_rate))

  write(*, '("plev (Pa)  play (Pa)    Tlay (K) zlay (m)  flux_dir/dn/up")')
  write(*, '(f7.0, 33x, 4(5x, f5.0), "  dir")')  p_lev(1,1), flux_dir(:,1)
  write(*, '(      40x, 4(5x, f5.0), "  dn" )')              flux_dn (:,1)
  write(*, '(      40x, 4(5x, f5.0), "  up" )')              flux_up (:,1)
  do i = 1, nlay
    write(*, '(6x, 3(4x, f7.0), 4(4x, f5.2))') p_lay(1,i), t_lay(1,i), z_lay(1,i), &
                                               heating_rate(:,i) * 86400._wp
    write(*, '(f7.0, 33x, 4(5x, f5.2), "  dir")')   p_lev(1,i+1), flux_dir(:,i+1)
    write(*, '(      40x, 4(5x, f5.2), "  dn" )')                 flux_dn (:,i+1)
    write(*, '(      40x, 4(5x, f5.2), "  up" )')                 flux_up (:,i+1)
  end do



end program test_solar_zenith_angle
