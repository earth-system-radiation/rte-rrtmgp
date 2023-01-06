! ----------------------------------------------------------------------------------
! Test aerosol optics with optical properties derived from lookup tables.
! ----------------------------------------------------------------------------------
program test_aerosol_optics
    use mo_rte_kind,                    only: wp
    use mo_aerosol_optics,              only: ty_aerosol_optics
    use mo_optical_props,               only: ty_optical_props_arry, &
                                              ty_optical_props_1scl, &
                                              ty_optical_props_2str
    use mo_load_aerosol_coefficients,   only: load_aero_lutcoeff, &
                                              read_aero_state, write_aero_op, is_lw
    use mo_rrtmgp_constants,            only: m_h2o, m_dry, grav

    implicit none

    ! ----------------------------------------------------------------------------------
    integer                            :: ncol, nlay, nlev, nbnd, nbin
    integer                            :: i, ibnd

    integer, dimension(:,:), allocatable :: aero_type
                                            ! MERRA2/GOCART aerosol type
                                            ! 0: no aerosol
                                            ! 1: dust
                                            ! 2: sea salt
                                            ! 3: sulfate
                                            ! 4: black carbon, hydrophobic
                                            ! 5: black carbon, hydrophilic
                                            ! 6: organic carbon, hydrophobic
                                            ! 7: organic carbon, hydrophilic
    real(wp), dimension(:,:), allocatable :: aero_size
                                            ! Aerosol size for dust and sea salt
                                            ! Allowable range: 0 - 10 microns
    real(wp), dimension(:,:), allocatable :: aero_mmr
                                            ! Aerosol mass mixing ratio
    real(wp), dimension(:,:), allocatable :: aero_mass
                                            ! Aerosol mass column (g/m2)
    real(wp), dimension(:,:), allocatable :: relhum
                                            ! Relative humidity (fraction)
    real(wp), dimension(:,:), allocatable :: p_lay   ! layer pressure (Pa)
    real(wp), dimension(:,:), allocatable :: p_lev   ! level pressure (Pa)
    real(wp), dimension(:,:), allocatable :: t_lay   ! layer temperature (K)
    real(wp), dimension(:,:), allocatable :: vmr_h2o ! water volume mixing ratio

    type(ty_aerosol_optics)                   :: aerosol_spec
    class(ty_optical_props_arry), allocatable :: aerosol_optical_props

    ! Atmospheric state input and aerosol optical property output
    character(len=128)            :: aerosol_state_file = 'rrtmgp-inputs-outputs-aerosol.nc'
    ! Aerosol property lookup table inputs
    character(len=128)            :: aerosol_coeff_file = 'aero_coefficients.nc'

    integer                       :: l
    ! ----------------------------------------------------------------------------------
    !
    ! Input aerosol optical property coefficients from lookup tables
    !   and call aerosol optics intitialization
    ! Lookup tables are used if they are available
    !
    print *, "Deriving aerosol optics from lookup tables"
    call load_aero_lutcoeff (aerosol_spec, aerosol_coeff_file)

    !
    ! Input aerosol optical properties from NetCDF file
    !
    call read_aero_state(aerosol_state_file, p_lay, p_lev, t_lay, vmr_h2o)
!    call stop_on_err(aerosol_spec%get_aerosol_size_bin(aero_size, nbin))

    ! Define output array sizes
    ncol = size(p_lay,dim=1)
    nlay = size(p_lay,dim=2)
    nlev = size(p_lev,dim=2)

    !
    ! Generate layer relative humidity for aerosol optics calculations
    !
    call get_relhum(ncol, nlay, p_lay, t_lay, vmr_h2o, relhum)
    ! ----------------------------------------------------------------------------------
    ! Set sample user-requested aerosol specification
    !
    ! aerosol type
    !
    allocate(aero_type(ncol,nlay))
    aero_type(:ncol,:nlay) = 0
    aero_type(:ncol,5:7) = 1
    aero_type(:ncol,1) = 2
    !
    ! aerosol size
    !
    allocate(aero_size(ncol,nlay))
    aero_size(:ncol,:nlay) = 0.0_wp
    aero_size(:ncol,5:7) = 2.5_wp
    aero_size(:ncol,1) = 1.5_wp
    !
    ! aerosol mmr
    !
    allocate(aero_mmr(ncol,nlay))
    aero_mmr(:ncol,:nlay) = 0.0_wp
    aero_mmr(:ncol,5:7) = 1.e-7_wp
    aero_mmr(:ncol,1) = 1.e-7_wp
    !
    ! Convert aerosol mmr to aerosol mass 
    !
    call aero_mmr2mass(ncol, nlay, nlev, p_lev, vmr_h2o, aero_mmr, aero_mass)

    ! ----------------------------------------------------------------------------------
    if (is_lw(trim(aerosol_state_file))) then
       print*, 'Calculating longwave aerosol optical depths'
       allocate(ty_optical_props_1scl::aerosol_optical_props)
    else
       print*, 'Calculating shortwave aerosol optical depths'
       allocate(ty_optical_props_2str::aerosol_optical_props)
    endif

    select type (aerosol_optical_props)
      type is (ty_optical_props_1scl)       ! two-stream calculation: tau only
        call stop_on_err(aerosol_optical_props%alloc_1scl(ncol, nlay, aerosol_spec))
      type is (ty_optical_props_2str)       ! two-stream calculation: tau, ssa, g
        call stop_on_err(aerosol_optical_props%alloc_2str(ncol, nlay, aerosol_spec))
    end select

    call stop_on_err(aerosol_spec%aerosol_optics(aero_type, aero_size, aero_mass, relhum, &
                     aerosol_optical_props) )

    call write_aero_op(trim(aerosol_state_file), aerosol_optical_props%get_nband(), aerosol_optical_props)

    ! all done
    print *, 'aerosol optics test end.'

contains
  ! --------------------------------------------------------------------------------------
  !
  ! Calculate layer relative humidity for aerosol optics calculations
  !
  subroutine get_relhum(ncol, nlay, p_lay, t_lay, vmr_h2o, relhum)

    integer,  intent(in) :: ncol, nlay
    real(wp), intent(in) :: p_lay(ncol,nlay)    ! layer pressure (Pa)
    real(wp), intent(in) :: t_lay(ncol,nlay)    ! layer temperature (K)
    real(wp), intent(in) :: vmr_h2o(ncol,nlay)  ! water volume mixing ratio

    real(wp), dimension(:,:), allocatable, intent(inout) :: relhum ! relative humidity (fraction, 0-1)

    ! Local variables 
    integer :: i, k

    real(wp) :: mmr_h2o             ! water mass mixing ratio
    real(wp) :: q_lay               ! water specific humidity
    real(wp) :: q_lay_min, q_tmp, es_tmp
    real(wp) :: mwd, t_ref, rh

    ! Set constants
    mwd = m_h2o/m_dry            ! ratio of water to dry air molecular weights
    t_ref = 273.16_wp            ! reference temperature (K)
    q_lay_min = 1.e-7_wp         ! minimum water mass mixing ratio
    ! -------------------

    ! Allocate output relative humidity array
    allocate(relhum(ncol,nlay))

    ! Derive layer virtual temperature
    do i = 1, ncol 
       do k = 1, nlay
          ! Convert h2o vmr to mmr
          mmr_h2o = vmr_h2o(i,k) * mwd
          q_lay = mmr_h2o / (1 + mmr_h2o)
          q_tmp = max(q_lay_min, q_lay)
          es_tmp = exp( (17.67_wp * (t_lay(i,k)-t_ref)) / (t_lay(i,k)-29.65_wp) )
          rh = (0.263_wp * p_lay(i,k) * q_tmp) / es_tmp
          ! Convert rh from percent to fraction
          relhum(i,k) = 0.01 * rh
       enddo
    enddo
    
  end subroutine get_relhum
  ! --------------------------------------------------------------------------------------
  !
  ! Convert layer aerosol mmr to mass for aerosol optics calculations
  !
  subroutine aero_mmr2mass(ncol, nlay, nlev, p_lev, vmr_h2o, aero_mmr, aero_mass)

    integer,  intent(in) :: ncol, nlay, nlev
    real(wp), intent(in) :: p_lev(ncol,nlay)    ! level pressure (Pa)
    real(wp), intent(in) :: vmr_h2o(ncol,nlay)  ! water volume mixing ratio
    real(wp), intent(in) :: aero_mmr(ncol,nlay) ! aerosol mass mixing ratio

    real(wp), dimension(:,:), allocatable, intent(inout) :: aero_mass ! aerosol column integrated mass (g/m2)

    ! Local variables 
    integer :: i, k

    real(wp) :: mmr_h2o           ! water mass mixing ratio
    real(wp) :: mwd, rg, mmr2mass
    real(wp) :: q_lay             ! water specific humidity
    real(wp) :: p_del, p_deldry   ! layer pressure thickness; moist, dry (Pa, kg m-1 s-2)

    ! Set constants
    mwd = m_h2o/m_dry             ! ratio of water to dry air molecular weights
    rg = 1._wp/grav               ! inverse of accelration of gravity (s2 m-1)
    ! -------------------

    ! Allocate output aerosol mass array
    allocate(aero_mass(ncol,nlay))

    ! Convert aerosol mmr to mass over layer (g/m2)
    do i = 1, ncol 
       do k = 1, nlay
          ! Convert h2o vmr to mmr
          mmr_h2o = vmr_h2o(i,k) * mwd
          q_lay = mmr_h2o / (1._wp + mmr_h2o)
          p_del = p_lev(i,k) - p_lev(i,k+1)
          p_deldry = p_del * (1._wp - q_lay)
          mmr2mass = rg * p_deldry
          ! Derive mass, convert from kg/m2 to g/m2
          aero_mass(i,k) = aero_mmr(i,k) * mmr2mass * 1.e3_wp
       enddo
    enddo
    
  end subroutine aero_mmr2mass
! -----------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop
    !
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then
      write (error_unit,*) trim(msg)
      write (error_unit,*) "test_aerosol_optics stopping"
      stop
    end if
  end subroutine

end program test_aerosol_optics
