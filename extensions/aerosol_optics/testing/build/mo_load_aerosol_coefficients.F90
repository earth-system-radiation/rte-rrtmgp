module mo_load_aerosol_coefficients
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_aerosol_optics,  only: ty_aerosol_optics
  use mo_simple_netcdf, only: read_field, read_string, var_exists, get_dim_size, &
                              write_field, create_dim, create_var
  use netcdf

  implicit none
  private
  public :: load_aero_lutcoeff, read_aero_state, write_aero_op
  public :: read_aero_op_nbnd, write_aero_op_ngpt, is_lw, is_sw
  ! ----------------------------------------------------------------------------------

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read aerosol optical property LUT coefficients from NetCDF file
  !
  subroutine load_aero_lutcoeff(aerosol_spec, aero_coeff_file)
    class(ty_aerosol_optics),   intent(inout) :: aerosol_spec
    character(len=*),           intent(in   ) :: aero_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrh, nbin

    real(wp), dimension(:,:),   allocatable :: band_lims_wvn
    ! LUT interpolation arrays
    real(wp), dimension(:),     allocatable :: aero_bin_lo      ! aerosol size bin lower bounds
    real(wp), dimension(:),     allocatable :: aero_bin_hi      ! aerosol size bin upper bounds
    real(wp), dimension(:),     allocatable :: aero_bin_eff     ! aerosol size bin effective radius
    real(wp), dimension(:),     allocatable :: aero_rh          ! aerosol relative humidity values
    ! LUT coefficients 
    real(wp), dimension(:,:),   allocatable :: aero_dust_ext    ! extinction: dust
    real(wp), dimension(:,:),   allocatable :: aero_dust_ssa    ! single scattering albedo: dust
    real(wp), dimension(:,:),   allocatable :: aero_dust_asy    ! asymmetry parameter: dust
    real(wp), dimension(:,:,:), allocatable :: aero_salt_ext    ! extinction: sea salt
    real(wp), dimension(:,:,:), allocatable :: aero_salt_ssa    ! single scattering albedo: sea salt
    real(wp), dimension(:,:,:), allocatable :: aero_salt_asy    ! asymmetry parameter: sea salt
    real(wp), dimension(:,:),   allocatable :: aero_sulf_ext    ! extinction: sulfate
    real(wp), dimension(:,:),   allocatable :: aero_sulf_ssa    ! single scattering albedo: sulfate
    real(wp), dimension(:,:),   allocatable :: aero_sulf_asy    ! asymmetry parameter: sulfate
    real(wp), dimension(:),     allocatable :: aero_bcar_ext    ! extinction: black carbon, hydrophobic
    real(wp), dimension(:),     allocatable :: aero_bcar_ssa    ! single scattering albedo: black carbon, hydrophobic
    real(wp), dimension(:),     allocatable :: aero_bcar_asy    ! asymmetry parameter: black carbon, hydrophobic
    real(wp), dimension(:,:),   allocatable :: aero_bcar_rh_ext ! extinction: black carbon, hydrophilic
    real(wp), dimension(:,:),   allocatable :: aero_bcar_rh_ssa ! single scattering albedo: black carbon, hydrophilic
    real(wp), dimension(:,:),   allocatable :: aero_bcar_rh_asy ! asymmetry parameter: black carbon, hydrophilic
    real(wp), dimension(:),     allocatable :: aero_ocar_ext    ! extinction: organic carbon, hydrophobic
    real(wp), dimension(:),     allocatable :: aero_ocar_ssa    ! single scattering albedo: organic carbon, hydrophobic
    real(wp), dimension(:),     allocatable :: aero_ocar_asy    ! asymmetry parameter: organic carbon, hydrophobic
    real(wp), dimension(:,:),   allocatable :: aero_ocar_rh_ext ! extinction: organic carbon, hydrophilic
    real(wp), dimension(:,:),   allocatable :: aero_ocar_rh_ssa ! single scattering albedo: organic carbon, hydrophilic
    real(wp), dimension(:,:),   allocatable :: aero_ocar_rh_asy ! asymmetry parameter: organic carbon, hydrophilic
    ! -----------------
    ! Open aerosol optical property coefficient file
    if(nf90_open(trim(aero_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_aero_lutcoeff(): can't open file " // trim(aero_coeff_file))

    ! Read LUT coefficient dimensions
    nband = get_dim_size(ncid,'nband')
    nrh   = get_dim_size(ncid,'nrh')
    nbin  = get_dim_size(ncid,'nbin')

    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)

    ! Allocate aerosol property lookup table interpolation arrays
    allocate(aero_bin_lo(nbin), &
             aero_bin_hi(nbin), &
             aero_bin_eff(nbin), &
             aero_rh(nrh))

    ! Allocate aerosol property lookup table input arrays
    allocate(aero_dust_ext(nbin, nband), &
             aero_dust_ssa(nbin, nband), &
             aero_dust_asy(nbin, nband), &
             aero_salt_ext(nbin, nrh, nband), &
             aero_salt_ssa(nbin, nrh, nband), &
             aero_salt_asy(nbin, nrh, nband), &
             aero_sulf_ext(nrh, nband), &
             aero_sulf_ssa(nrh, nband), &
             aero_sulf_asy(nrh, nband), &
             aero_bcar_ext(nband), &
             aero_bcar_ssa(nband), &
             aero_bcar_asy(nband), &
             aero_bcar_rh_ext(nrh, nband), &
             aero_bcar_rh_ssa(nrh, nband), &
             aero_bcar_rh_asy(nrh, nband), &
             aero_ocar_ext(nband), &
             aero_ocar_ssa(nband), &
             aero_ocar_asy(nband), &
             aero_ocar_rh_ext(nrh, nband), &
             aero_ocar_rh_ssa(nrh, nband), &
             aero_ocar_rh_asy(nrh, nband))

    ! Read LUT interpolation arryas
    aero_bin_lo = read_field(ncid, 'aero_bin_lo',  nbin)
    aero_bin_hi = read_field(ncid, 'aero_bin_hi',  nbin)
    aero_bin_eff= read_field(ncid, 'aero_bin_eff',  nbin)
    aero_rh = read_field(ncid, 'aero_rh',  nrh)

    ! Read LUT coefficients
    aero_dust_ext = read_field(ncid, 'aero_dust_ext',  nbin, nband)
    aero_dust_ssa = read_field(ncid, 'aero_dust_ssa',  nbin, nband)
    aero_dust_asy = read_field(ncid, 'aero_dust_asy',  nbin, nband)
    aero_salt_ext = read_field(ncid, 'aero_salt_ext',  nbin, nrh, nband)
    aero_salt_ssa = read_field(ncid, 'aero_salt_ssa',  nbin, nrh, nband)
    aero_salt_asy = read_field(ncid, 'aero_salt_asy',  nbin, nrh, nband)
    aero_sulf_ext = read_field(ncid, 'aero_sulf_ext',  nrh, nband)
    aero_sulf_ssa = read_field(ncid, 'aero_sulf_ssa',  nrh, nband)
    aero_sulf_asy = read_field(ncid, 'aero_sulf_asy',  nrh, nband)
    aero_bcar_ext = read_field(ncid, 'aero_bcar_ext',  nband)
    aero_bcar_ssa = read_field(ncid, 'aero_bcar_ssa',  nband)
    aero_bcar_asy = read_field(ncid, 'aero_bcar_asy',  nband)
    aero_bcar_rh_ext = read_field(ncid, 'aero_bcar_rh_ext',  nrh, nband)
    aero_bcar_rh_ssa = read_field(ncid, 'aero_bcar_rh_ssa',  nrh, nband)
    aero_bcar_rh_asy = read_field(ncid, 'aero_bcar_rh_asy',  nrh, nband)
    aero_ocar_ext = read_field(ncid, 'aero_ocar_ext',  nband)
    aero_ocar_ssa = read_field(ncid, 'aero_ocar_ssa',  nband)
    aero_ocar_asy = read_field(ncid, 'aero_ocar_asy',  nband)
    aero_ocar_rh_ext = read_field(ncid, 'aero_ocar_rh_ext',  nrh, nband)
    aero_ocar_rh_ssa = read_field(ncid, 'aero_ocar_rh_ssa',  nrh, nband)
    aero_ocar_rh_asy = read_field(ncid, 'aero_ocar_rh_asy',  nrh, nband)

    ncid = nf90_close(ncid)

    call stop_on_err(aerosol_spec%load(band_lims_wvn, &
                                  aero_bin_lo, aero_bin_hi, aero_bin_eff, aero_rh, &
                                  aero_dust_ext, aero_dust_ssa, aero_dust_asy, &
                                  aero_salt_ext, aero_salt_ssa, aero_salt_asy, &
                                  aero_sulf_ext, aero_sulf_ssa, aero_sulf_asy, &
                                  aero_bcar_ext, aero_bcar_ssa, aero_bcar_asy, &
                                  aero_bcar_rh_ext, aero_bcar_rh_ssa, aero_bcar_rh_asy, &
                                  aero_ocar_ext, aero_ocar_ssa, aero_ocar_asy, &
                                  aero_ocar_rh_ext, aero_ocar_rh_ssa, aero_ocar_rh_asy))

  end subroutine load_aero_lutcoeff
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_aero_state(filename, p_lay, p_lev, t_lay, vmr_h2o)
    character(len=*),                      intent(in   ) :: filename
    real(wp), dimension(:,:), allocatable, intent(inout) :: p_lay      ! layer pressure 
    real(wp), dimension(:,:), allocatable, intent(inout) :: p_lev      ! level pressure 
    real(wp), dimension(:,:), allocatable, intent(inout) :: t_lay      ! layer temperature
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_h2o    ! water volume mixing ratio

    ! -----------------
    integer :: ncid
    integer :: ncol, nlay, nlev
    ! -----------------
    ! Open aerosol optical property coefficient file
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("read_aero_state(): can't open file " // trim(filename))

    ncol = get_dim_size(ncid,'col')
    nlay = get_dim_size(ncid,'lay')
    nlev = get_dim_size(ncid,'lev')
    ! Allocate aerosol physical property arrays
    allocate(p_lay(ncol,nlay), &
             p_lev(ncol,nlev), &
             t_lay(ncol,nlay), &
             vmr_h2o(ncol,nlay))

!    icergh     = read_field(ncid, 'icergh')
    p_lay    = read_field(ncid, 'p_lay', ncol, nlay)
    p_lev    = read_field(ncid, 'p_lev', ncol, nlev)
    t_lay    = read_field(ncid, 't_lay', ncol, nlay)
    vmr_h2o  = read_field(ncid, 'vmr_h2o', ncol, nlay)

    ncid = nf90_close(ncid)

  end subroutine read_aero_state
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write aerosol optical properties (nbnd)
  !
  subroutine write_aero_op(filename, nbnd, aerosol_optical_props)
    character(len=*),             intent(in) :: filename
    integer,                      intent(in) :: nbnd
    class(ty_optical_props_arry), intent(in) :: aerosol_optical_props
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nmom
    ! -------------------
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_aero_op: can't open file " // trim(filename))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol   = get_dim_size(ncid, 'col')
    nlay   = get_dim_size(ncid, 'lay')
    nmom   = get_dim_size(ncid, 'mom')

    select type(aerosol_optical_props)
      type is (ty_optical_props_1scl)
        call create_var(ncid,    "tauaer",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call stop_on_err(write_field(ncid, "tauaer",  aerosol_optical_props%tau))
      type is (ty_optical_props_2str)
        call create_var(ncid,    "tauaer",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "ssaaer",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "asyaer",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call stop_on_err(write_field(ncid, "tauaer",  aerosol_optical_props%tau))
        call stop_on_err(write_field(ncid, "ssaaer",  aerosol_optical_props%ssa))
        call stop_on_err(write_field(ncid, "asyaer",  aerosol_optical_props%g  ))
      type is (ty_optical_props_nstr)
        call create_var(ncid,    "tauaer",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "ssaaer",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "paer",    ["mom ", "col ", "lay ", "band"], [nmom, ncol, nlay, nbnd])
        call stop_on_err(write_field(ncid, "tauaer",  aerosol_optical_props%tau))
        call stop_on_err(write_field(ncid, "ssaaer",  aerosol_optical_props%ssa))
        call stop_on_err(write_field(ncid, "paer",    aerosol_optical_props%p  ))
    end select

    ncid = nf90_close(ncid)
  end subroutine write_aero_op

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read aerosol optical properties (nbnd)
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_aero_op_nbnd(aerosol_optical_props, opt_props, filename, play, plev, tlay, vmr_h2o)

    class(ty_optical_props_arry),          intent(inout) :: aerosol_optical_props
    class(ty_optical_props),               intent(inout) :: opt_props
!    logical,                               intent(in   ) :: is_lw
    character(len=*),                      intent(in   ) :: filename
    real(wp), dimension(:,:), allocatable, intent(  out) :: play, plev, tlay, vmr_h2o

    ! -----------------
    integer :: ncid
    integer :: ncol, nlay, nlev, nmom, nbnd, ngpt, npai
    logical :: is_initialized

    integer, dimension(:,:), allocatable                 :: band_lims_gpt
    real(wp), dimension(:,:), allocatable                :: band_lims_wvn

    ! -----------------
    ! Open aerosol optical property coefficient file
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("read_aero_op_nbnd(): can't open file " // trim(filename))

    ncol = get_dim_size(ncid,'col')
    nlay = get_dim_size(ncid,'lay')
    nlev = get_dim_size(ncid,'lev')
    nmom = get_dim_size(ncid,'mom')
    npai = get_dim_size(ncid,'pair')

    allocate(play(ncol,nlay), plev(ncol,nlev), tlay(ncol,nlay))
    play = read_field(ncid, 'p_lay', ncol, nlay)
    plev = read_field(ncid, 'p_lev', ncol, nlev)
    tlay = read_field(ncid, 't_lay', ncol, nlay)

    if (is_lw(trim(filename))) then
       nbnd = get_dim_size(ncid, 'band')
       ngpt = get_dim_size(ncid, 'gpt')
       select type(aerosol_optical_props)
         type is (ty_optical_props_1scl)
            aerosol_optical_props%tau    = read_field(ncid, 'tauaer', ncol, nlay, nbnd)
         type is (ty_optical_props_2str)
            aerosol_optical_props%tau    = read_field(ncid, 'tauaer', ncol, nlay, nbnd)
            aerosol_optical_props%ssa    = read_field(ncid, 'ssaaer', ncol, nlay, nbnd)
            aerosol_optical_props%g      = read_field(ncid, 'asyaer', ncol, nlay, nbnd)
         type is (ty_optical_props_nstr)
            aerosol_optical_props%tau    = read_field(ncid, 'tauaer', ncol, nlay, nbnd)
            aerosol_optical_props%ssa    = read_field(ncid, 'ssaaer', ncol, nlay, nbnd)
            aerosol_optical_props%p      = read_field(ncid, 'paer', nmom, ncol, nlay, nbnd)
       end select
    end if

    if (.not. is_lw(trim(filename))) then
       nbnd = get_dim_size(ncid, 'band')
       ngpt = get_dim_size(ncid, 'gpt')
       select type(aerosol_optical_props)
         type is (ty_optical_props_1scl)
            aerosol_optical_props%tau    = read_field(ncid, 'tauaer', ncol, nlay, nbnd)
         type is (ty_optical_props_2str)
            aerosol_optical_props%tau    = read_field(ncid, 'tauaer', ncol, nlay, nbnd)
            aerosol_optical_props%ssa    = read_field(ncid, 'ssaaer', ncol, nlay, nbnd)
            aerosol_optical_props%g      = read_field(ncid, 'asyaer', ncol, nlay, nbnd)
         type is (ty_optical_props_nstr)
            aerosol_optical_props%tau    = read_field(ncid, 'tauaer', ncol, nlay, nbnd)
            aerosol_optical_props%ssa    = read_field(ncid, 'ssaaer', ncol, nlay, nbnd)
            aerosol_optical_props%p      = read_field(ncid, 'paer', nmom, ncol, nlay, nbnd)
       end select
    end if

    ncid = nf90_close(ncid)

  end subroutine read_aero_op_nbnd

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write aerosol optical properties (ngpt)
  !
  subroutine write_aero_op_ngpt(filename, aerosol_optical_props)
    character(len=*),                   intent(in) :: filename
!    logical,                            intent(in) :: is_lw
    class(ty_optical_props_arry), intent(in) :: aerosol_optical_props
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngptlw, ngptsw, nmom
    ! -------------------
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_aero_op_ngpt: can't open file " // trim(filename))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol   = get_dim_size(ncid, 'col')
    nlay   = get_dim_size(ncid, 'lay')
    nmom   = get_dim_size(ncid, 'mom')

    if (is_lw(trim(filename))) then
       ngptlw = get_dim_size(ncid, 'gpt')
       select type(aerosol_optical_props)
         type is (ty_optical_props_1scl)
           call create_var(ncid,    "tauaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call stop_on_err(write_field(ncid, "tauaerosol",  aerosol_optical_props%tau ))
         type is (ty_optical_props_2str)
           call create_var(ncid,    "tauaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "ssaaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "asyaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call stop_on_err(write_field(ncid, "tauaerosol",  aerosol_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssaaerosol",  aerosol_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "asyaerosol",  aerosol_optical_props%g ))
         type is (ty_optical_props_nstr)
           call create_var(ncid,    "tauaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "ssaaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "paerosol",    ["mom", "col", "lay", "gpt"], [nmom, ncol, nlay, ngptlw])
           call stop_on_err(write_field(ncid, "tauaerosol",  aerosol_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssaaerosol",  aerosol_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "paerosol",    aerosol_optical_props%p ))
       end select
   endif

    if (.not. is_lw(trim(filename))) then
       ngptsw = get_dim_size(ncid, 'gpt')
       select type(aerosol_optical_props)
         type is (ty_optical_props_1scl)
           call create_var(ncid,    "tauaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call stop_on_err(write_field(ncid, "tauaerosol",  aerosol_optical_props%tau ))
         type is (ty_optical_props_2str)
           call create_var(ncid,    "tauaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "ssaaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "asyaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call stop_on_err(write_field(ncid, "tauaerosol",  aerosol_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssaaerosol",  aerosol_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "asyaerosol",  aerosol_optical_props%g ))
         type is (ty_optical_props_nstr)
           call create_var(ncid,    "tauaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "ssaaerosol",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "paerosol",    ["mom", "col", "lay", "gpt"], [nmom, ncol, nlay, ngptsw])
           call stop_on_err(write_field(ncid, "tauaerosol",  aerosol_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssaaerosol",  aerosol_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "paerosol",    aerosol_optical_props%p ))
       end select
    endif

    ncid = nf90_close(ncid)
  end subroutine write_aero_op_ngpt

  !------------------------------------------------------------------------------------------------------
  !
  ! Does this file contain variables needed to do SW calculations ?
  !
  function is_sw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_sw

    integer :: ncid

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("is_sw: can't find file " // trim(fileName))

    is_sw = var_exists(ncid, 'solar_zenith_angle')
    ncid = nf90_close(ncid)
  end function is_sw

  !------------------------------------------------------------------------------------------------------
  function is_lw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw

    is_lw = .not. is_sw(fileName)
  end function is_lw

  ! -----------------------------------------------------------------------------------
    subroutine stop_on_err(msg)
      !
      ! Print error message and stop
      !
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      if(len_trim(msg) > 0) then
        write (error_unit,*) trim(msg)
        stop
      end if
    end subroutine

end module
