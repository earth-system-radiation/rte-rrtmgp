! This code is part of Radiative Transfer for Energetics (RTE)
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
! Provides aerosol optical properties as a function of effective radius for the RRTMGP bands
!   Based on climatoligical aerosol optical properties used in MERRA2 as derived from
!     the GOCART model for 15 aerosol types, including five dust types, five sea salt types,
!     one sulfate type, and both hydrophobic and hydrophilic black carbon and organic carbon.
!   Input aerosol optical data are stored in look-up tables.
!
!   References for gocart interactive aerosols:                        
!     Chin et al., 2000 - jgr, v105, 24671-24687                       
!     Colarco et al., 2010 - jgr, v115, D14207                         
!                                                                      
!   References for merra2 aerosol reanalysis:                          
!     Randles et al., 2017 - jclim, v30, 6823-6850                     
!     Buchard et al., 2017 - jclim, v30, 6851-6871                     
!                                                                      
! The class can be used as-is but is also intended as an example of how to extend the RTE framework
! -------------------------------------------------------------------------------------------------

module mo_aerosol_optics
  use mo_rte_kind,      only: wp, wl
  use mo_rte_config,    only: check_values, check_extents
  use mo_rte_util_array,only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  implicit none
  private
  ! -----------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_aerosol_optics
    private
    !
    ! Aerosol type identifiers (MERRA2/GOCART aerosol climotology)
    !
    ! dust
    integer :: aero_dust = 1
    ! sea-salt
    integer :: aero_salt = 2
    ! sulfate
    integer :: aero_sulf = 3
    ! black carbon (hydrophilic)
    integer :: aero_bcar_rh = 4
    ! black carbon (hydrophobic)
    integer :: aero_bcar = 5
    ! organic carbon (hydrophilic)
    integer :: aero_ocar_rh = 6
    ! organic carbon (hydrophobic)
    integer :: aero_ocar = 7
    !
    ! Lookup table information
    !
    ! Table upper and lower size bin limits and effective (center) values
    real(wp),dimension(:), allocatable :: aero_bin_hi(:)
    real(wp),dimension(:), allocatable :: aero_bin_lo(:)
    real(wp),dimension(:), allocatable :: aero_bin_eff(:)
    ! Table relative humidity values
    real(wp),dimension(:), allocatable :: aero_rh(:)
    integer, dimension(:), allocatable :: rh_nsteps
    real(wp),dimension(:), allocatable :: rh_step_sizes, rh_offsets
    !
    ! The tables themselves.
    !
    real(wp), dimension(:,:  ), allocatable :: aero_dust_ext, aero_dust_ssa, aero_dust_asy          ! (nbin, nbnd)
    real(wp), dimension(:,:,:), allocatable :: aero_salt_ext, aero_salt_ssa, aero_salt_asy          ! (nbin, nrh, nbnd)
    real(wp), dimension(:,:  ), allocatable :: aero_sulf_ext, aero_sulf_ssa, aero_sulf_asy          ! (nrh, nbnd)
    real(wp), dimension(:    ), allocatable :: aero_bcar_ext, aero_bcar_ssa, aero_bcar_asy          ! (nbnd)
    real(wp), dimension(:,:  ), allocatable :: aero_bcar_rh_ext, aero_bcar_rh_ssa, aero_bcar_rh_asy ! (nrh, nbnd)
    real(wp), dimension(:    ), allocatable :: aero_ocar_ext, aero_ocar_ssa, aero_ocar_asy          ! (nbnd)
    real(wp), dimension(:,:  ), allocatable :: aero_ocar_rh_ext, aero_ocar_rh_ssa, aero_ocar_rh_asy ! (nrh, nbnd)
    !
    ! -----
  contains
    generic,   public :: load  => load_lut
    procedure, public :: finalize
    procedure, public :: aerosol_optics

    ! Internal procedures
    procedure, private :: load_lut
  end type ty_aerosol_optics

contains
  ! ------------------------------------------------------------------------------
  !
  ! Routines to load data needed for aerosol optics calculations from lookup-tables.
  !
  ! ------------------------------------------------------------------------------
  function load_lut(this, band_lims_wvn, &
                    aero_bin_lo, aero_bin_hi, aero_bin_eff, aero_rh, &
                    aero_dust_ext, aero_dust_ssa, aero_dust_asy, &
                    aero_salt_ext, aero_salt_ssa, aero_salt_asy, &
                    aero_sulf_ext, aero_sulf_ssa, aero_sulf_asy, &
                    aero_bcar_ext, aero_bcar_ssa, aero_bcar_asy, &
                    aero_bcar_rh_ext, aero_bcar_rh_ssa, aero_bcar_rh_asy, &
                    aero_ocar_ext, aero_ocar_ssa, aero_ocar_asy, &
                    aero_ocar_rh_ext, aero_ocar_rh_ssa, aero_ocar_rh_asy) &
                    result(error_msg)
    class(ty_aerosol_optics),   intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn ! spectral discretization
    ! Lookup table interpolation constants
    real(wp), dimension(:),     intent(in   ) :: aero_bin_lo   ! aerosol lut size bin lower boundary
    real(wp), dimension(:),     intent(in   ) :: aero_bin_hi   ! aerosol lut size bin upper boundary
    real(wp), dimension(:),     intent(in   ) :: aero_bin_eff  ! aerosol lut size bin effective size
    real(wp), dimension(:),     intent(in   ) :: aero_rh       ! relative humidity LUT dimension values
    ! LUT coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter for aerosol types
    real(wp), dimension(:,:),   intent(in)    :: aero_dust_ext, aero_dust_ssa, aero_dust_asy
    real(wp), dimension(:,:,:), intent(in)    :: aero_salt_ext, aero_salt_ssa, aero_salt_asy
    real(wp), dimension(:,:),   intent(in)    :: aero_sulf_ext, aero_sulf_ssa, aero_sulf_asy
    real(wp), dimension(:),     intent(in)    :: aero_bcar_ext, aero_bcar_ssa, aero_bcar_asy
    real(wp), dimension(:,:),   intent(in)    :: aero_bcar_rh_ext, aero_bcar_rh_ssa, aero_bcar_rh_asy
    real(wp), dimension(:),     intent(in)    :: aero_ocar_ext, aero_ocar_ssa, aero_ocar_asy
    real(wp), dimension(:,:),   intent(in)    :: aero_ocar_rh_ext, aero_ocar_rh_ssa, aero_ocar_rh_asy
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer               :: nbin, nrh, nband

    error_msg = this%init(band_lims_wvn, name="RRTMGP aerosol optics")
    !
    ! LUT coefficient dimensions
    !
    nbin  = size(aero_salt_ext,dim=1)
    nrh   = size(aero_salt_ext,dim=2)
    nband = size(aero_salt_ext,dim=3)
    !
    ! Error checking
    !
!    if(nband /= this%get_nband()) &
!      error_msg = "aerosol_optics%init(): number of bands inconsistent between lookup tables, spectral discretization"
    error_msg = ''
    if(size(aero_dust_ext, 1) /= nbin) &
      error_msg = "aerosol_optics%init(): array aero_dust_ext has the wrong number of size bins"
    if(size(aero_dust_ext, 2) /= nband) &
      error_msg = "aerosol_optics%init(): array aero_dust_ext has the wrong number of bands"
    if(size(aero_salt_ext, 3) /= nband) &
      error_msg = "aerosol_optics%init(): array aero_salt_ext has the wrong number of bands"
    if(size(aero_sulf_ext, 2) /= nband) &
      error_msg = "aerosol_optics%init(): array aero_sulf_ext has the wrong number of bands"
    if(size(aero_bcar_ext, 1) /= nband) &
      error_msg = "aerosol_optics%init(): array aero_bcar_ext has the wrong number of bands"
    if(size(aero_ocar_ext, 1) /= nband) &
      error_msg = "aerosol_optics%init(): array aero_ocar_ext has the wrong number of bands"
    if(.not. extents_are(aero_salt_ext, nbin, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_salt_ext isn't consistently sized"
    if(.not. extents_are(aero_salt_ssa, nbin, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_salt_ssa isn't consistently sized"
    if(.not. extents_are(aero_salt_asy, nbin, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_salt_asy isn't consistently sized"
    if(.not. extents_are(aero_sulf_ext, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_sulf_ext isn't consistently sized"
    if(.not. extents_are(aero_sulf_ssa, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_sulf_ssa isn't consistently sized"
    if(.not. extents_are(aero_sulf_asy, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_sulf_asy isn't consistently sized"
    if(.not. extents_are(aero_bcar_rh_ext, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_bcar_rh_ext isn't consistently sized"
    if(.not. extents_are(aero_bcar_rh_ssa, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_bcar_rh_ssa isn't consistently sized"
    if(.not. extents_are(aero_bcar_rh_asy, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_bcar_rh_asy isn't consistently sized"
    if(.not. extents_are(aero_ocar_rh_ext, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_ocar_rh_ext isn't consistently sized"
    if(.not. extents_are(aero_ocar_rh_ssa, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_ocar_rh_ssa isn't consistently sized"
    if(.not. extents_are(aero_ocar_rh_asy, nrh, nband)) &
      error_msg = "aerosol_optics%init(): array aero_ocar_rh_asy isn't consistently sized"
    if(error_msg /= "") return

    ! Allocate LUT parameters
    allocate(this%aero_bin_lo(nbin), &
             this%aero_bin_hi(nbin), &
             this%aero_bin_eff(nbin))
    this%aero_bin_lo = aero_bin_lo
    this%aero_bin_hi = aero_bin_hi
    this%aero_bin_eff = aero_bin_eff
    allocate(this%aero_rh(nrh), &
             this%rh_nsteps(2), &
             this%rh_step_sizes(2), &
             this%rh_offsets(2))
    this%aero_rh = aero_rh
    this%rh_nsteps(1) = 17
    this%rh_nsteps(2) = 20
    this%rh_step_sizes(1) = 0.05_wp
    this%rh_step_sizes(2) = 0.01_wp
    this%rh_offsets(1) = 0.0_wp
    this%rh_offsets(2) = 0.8_wp
    ! Allocate LUT coefficients
    allocate(this%aero_dust_ext(nbin, nband), &
             this%aero_dust_ssa(nbin, nband), &
             this%aero_dust_asy(nbin, nband), &
             this%aero_salt_ext(nbin, nrh, nband), &
             this%aero_salt_ssa(nbin, nrh, nband), &
             this%aero_salt_asy(nbin, nrh, nband), &
             this%aero_sulf_ext(nrh, nband), &
             this%aero_sulf_ssa(nrh, nband), &
             this%aero_sulf_asy(nrh, nband), &
             this%aero_bcar_ext(nband), &
             this%aero_bcar_ssa(nband), &
             this%aero_bcar_asy(nband), &
             this%aero_bcar_rh_ext(nrh, nband), &
             this%aero_bcar_rh_ssa(nrh, nband), &
             this%aero_bcar_rh_asy(nrh, nband), &
             this%aero_ocar_ext(nband), &
             this%aero_ocar_ssa(nband), &
             this%aero_ocar_asy(nband), &
             this%aero_ocar_rh_ext(nrh, nband), &
             this%aero_ocar_rh_ssa(nrh, nband), &
             this%aero_ocar_rh_asy(nrh, nband))

    ! Load LUT coefficients
    !$acc kernels
    !$omp target
    this%aero_dust_ext = aero_dust_ext
    this%aero_dust_ssa = aero_dust_ssa
    this%aero_dust_asy = aero_dust_asy
    this%aero_salt_ext = aero_salt_ext
    this%aero_salt_ssa = aero_salt_ssa
    this%aero_salt_asy = aero_salt_asy
    this%aero_sulf_ext = aero_sulf_ext
    this%aero_sulf_ssa = aero_sulf_ssa
    this%aero_sulf_asy = aero_sulf_asy
    this%aero_bcar_ext = aero_bcar_ext
    this%aero_bcar_ssa = aero_bcar_ssa
    this%aero_bcar_asy = aero_bcar_asy
    this%aero_bcar_rh_ext = aero_bcar_rh_ext
    this%aero_bcar_rh_ssa = aero_bcar_rh_ssa
    this%aero_bcar_rh_asy = aero_bcar_rh_asy
    this%aero_ocar_ext = aero_ocar_ext
    this%aero_ocar_ssa = aero_ocar_ssa
    this%aero_ocar_asy = aero_ocar_asy
    this%aero_ocar_rh_ext = aero_ocar_rh_ext
    this%aero_ocar_rh_ssa = aero_ocar_rh_ssa
    this%aero_ocar_rh_asy = aero_ocar_rh_asy
    !$acc end kernels
    !$omp end target
    !
  end function load_lut
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Finalize
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine finalize(this)
    class(ty_aerosol_optics), intent(inout) :: this

    ! Lookup table aerosol optics interpolation arrays
    if(allocated(this%aero_bin_lo)) then

      deallocate(this%aero_bin_lo, this%aero_bin_hi, this%aero_bin_eff, &
                 this%aero_rh)
    end if

    ! Lookup table aerosol optics coefficients
    if(allocated(this%aero_dust_ext)) then

      deallocate(this%aero_dust_ext, this%aero_dust_ssa, this%aero_dust_asy, &
                 this%aero_salt_ext, this%aero_salt_ssa, this%aero_salt_asy, &
                 this%aero_sulf_ext, this%aero_sulf_ssa, this%aero_sulf_asy, &
                 this%aero_bcar_ext, this%aero_bcar_ssa, this%aero_bcar_asy, &
                 this%aero_bcar_rh_ext, this%aero_bcar_rh_ssa, this%aero_bcar_rh_asy, &
                 this%aero_ocar_ext, this%aero_ocar_ssa, this%aero_ocar_asy, &
                 this%aero_ocar_rh_ext, this%aero_ocar_rh_ssa, this%aero_ocar_rh_asy)
    end if

  end subroutine finalize
  ! ------------------------------------------------------------------------------
  !
  ! Derive aerosol optical properties from provided aerosol input properties
  !
  ! ------------------------------------------------------------------------------
  !
  ! Compute single-scattering properties
  !
  function aerosol_optics(this, aero_type, aero_size, aero_mass, relhum, &
                          optical_props) result(error_msg)

    class(ty_aerosol_optics), &
              intent(in  ) :: this
    integer,  intent(in  ) :: aero_type(:,:)   ! MERRA2/GOCART aerosol type 
                                               ! Dimensions: (ncol,nlay)
                                               ! 1 = aero_dust    (dust)
                                               ! 2 = aero_salt    (salt)
                                               ! 3 = aero_sulf    (sulfate)
                                               ! 4 = aero_bcar_rh (black carbon, hydrophilic)
                                               ! 5 = aero_bcar    (black carbon, hydrophobic)
                                               ! 6 = aero_ocar_rh (organic carbon, hydrophilic)
                                               ! 7 = aero_ocar    (organic carbon, hydrophobic)
    real(wp), intent(in  ) :: aero_size(:,:)   ! aerosol size for dust and sea-salt
                                               ! Dimensions: (ncol,nlay)
    real(wp), intent(in  ) :: aero_mass(:,:)   ! aerosol mass column (g/m2)
                                               ! Dimensions: (ncol,nlay)
    real(wp), intent(in  ) :: relhum(:,:)      ! relative humidity (fraction, 0-1)
                                               ! Dimensions: (ncol,nlay)

    class(ty_optical_props_arry), &
              intent(inout) :: optical_props
                                               ! Dimensions: (ncol,nlay,nbnd)

    character(len=128)      :: error_msg
    ! ------- Local -------
    logical(wl), dimension(size(aero_type,1), size(aero_type,2)) :: aeromsk
    real(wp),    dimension(size(aero_type,1), size(aero_type,2), size(this%aero_dust_ext,2)) :: &
                 atau, ataussa, ataussag
    integer  :: ncol, nlay, nbnd, nbin, nrh
    integer  :: icol, ilay, ibnd, ibin
    ! scalars for total tau, tau*ssa
    real(wp) :: tau, taussa

    ! ----------------------------------------
    !
    ! Error checking
    !
    ! ----------------------------------------

    error_msg = ''
    if(.not.(allocated(this%aero_dust_ext))) then
      error_msg = 'aerosol optics: no data has been initialized'
      return
    end if

    ncol = size(aero_type,1)
    nlay = size(aero_type,2)
    nbin = size(this%aero_bin_lo,1)
    nrh = size(this%aero_rh,1)
    nbnd = size(this%aero_dust_ext,2)

    !
    ! Array sizes
    !
    if (check_extents) then
      if(size(aero_type, 1) /= ncol .or. size(aero_type, 2) /= nlay) &
        error_msg = "aerosol optics: aero_type has wrong extents"
      if(size(aero_size, 1) /= ncol .or. size(aero_size, 2) /= nlay) &
        error_msg = "aerosol optics: aero_size has wrong extents"
      if(size(aero_mass, 1) /= ncol .or. size(aero_mass, 2) /= nlay) &
        error_msg = "aerosol optics: aero_mass has wrong extents"
      if(size(relhum,  1) /= ncol .or. size(relhum,  2) /= nlay) &
        error_msg = "aerosol optics: relhumn has wrong extents"
      if(optical_props%get_ncol() /= ncol .or. optical_props%get_nlay() /= nlay) &
        error_msg = "aerosol optics: optical_props have wrong extents"
      if(error_msg /= "") return
    end if

    !
    ! Spectral consistency
    !
    if(check_values) then
      if(.not. this%bands_are_equal(optical_props)) &
        error_msg = "aerosol optics: optical properties don't have the same band structure"
      if(optical_props%get_nband() /= optical_props%get_ngpt() ) &
        error_msg = "aerosol optics: optical properties must be requested by band not g-points"
      if(error_msg /= "") return
    end if

    !$acc data copyin(aero_type, aero_size, aero_mass, relhum)                         &
    !$acc      create(atau, ataussa, ataussag) &
    !$acc      create(aeromsk)
    !$omp target data map(to:aero_type, aero_size, aero_mass, relhum) &
    !$omp map(alloc:atau, ataussa, ataussag) &
    !$omp map(alloc:aeromsk)
    !
    ! Aerosol mask; don't need aerosol optics if there's no aerosol
    !
    !$acc parallel loop gang vector default(none) collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do ilay = 1, nlay
      do icol = 1, ncol
        aeromsk(icol,ilay) = aero_type(icol,ilay) > 0
      end do
    end do

    !
    ! Aerosol size, relative humidity
    !
    if(check_values) then
      if(any_vals_outside(aero_size, aeromsk, this%aero_bin_lo(1), this%aero_bin_hi(nbin))) &
        error_msg = 'aerosol optics: requested aerosol size is out of bounds'
      if(any_vals_outside(relhum, 0._wp, 1._wp)) &
        error_msg = 'aerosol optics: relative humidity fraction is out of bounds'
    end if
    if(error_msg == "") then
      !
      !
      ! ----------------------------------------
      !
      ! The lookup tables determining extinction coefficient, single-scattering albedo,
      !   and asymmetry parameter g as a function of aerosol type, aerosol size and
      !   relative humidity.
      ! We compute the optical depth tau (= exintinction coeff * aerosol mass ) and the
      !    products tau*ssa and tau*ssa*g separately for each aerosol type requested.
      ! These are used to determine the aerosol optical properties.
      !
      if (allocated(this%aero_dust_ext)) then
        !
        ! Aerosol
        !
        call compute_all_from_table(ncol, nlay, nbnd, nbin, nrh,                                &
                                    aeromsk, aero_type, aero_size, aero_mass, relhum,           &
                                    this%aero_bin_lo, this%aero_bin_hi, this%aero_rh,           &
                                    this%rh_nsteps, this%rh_step_sizes, this%rh_offsets,        &
                                    this%aero_dust, this%aero_salt, this%aero_sulf,             &
                                    this%aero_bcar, this%aero_bcar_rh, this%aero_ocar, this%aero_ocar_rh, &
                                    this%aero_dust_ext, this%aero_dust_ssa, this%aero_dust_asy, &
                                    this%aero_salt_ext, this%aero_salt_ssa, this%aero_salt_asy, &
                                    this%aero_sulf_ext, this%aero_sulf_ssa, this%aero_sulf_asy, &
                                    this%aero_bcar_rh_ext, this%aero_bcar_rh_ssa, this%aero_bcar_rh_asy, &
                                    this%aero_bcar_ext, this%aero_bcar_ssa, this%aero_bcar_asy, &
                                    this%aero_ocar_rh_ext, this%aero_ocar_rh_ssa, this%aero_ocar_rh_asy, &
                                    this%aero_ocar_ext, this%aero_ocar_ssa, this%aero_ocar_asy, &
                                    atau, ataussa, ataussag)

      endif

      !
      ! Combine liquid and ice contributions into total cloud optical properties
      !   See also the increment routines in mo_optical_props_kernels
      !
      select type(optical_props)
      type is (ty_optical_props_1scl)
        !$acc parallel loop gang vector default(none) collapse(3) &
        !$acc               copyin(optical_props) copyout(optical_props%tau)
        !$omp target teams distribute parallel do simd collapse(3) &
        !$omp map(from:optical_props%tau)

        do ibnd = 1, nbnd
          do ilay = 1, nlay
            do icol = 1, ncol
              ! Absorption optical depth  = (1-ssa) * tau = tau - taussa
              optical_props%tau(icol,ilay,ibnd) = (atau(icol,ilay,ibnd) - ataussa(icol,ilay,ibnd))
            end do
          end do
        end do
      type is (ty_optical_props_2str)
        !$acc parallel loop gang vector default(none) collapse(3) &
        !$acc               copyin(optical_props) copyout(optical_props%tau, optical_props%ssa, optical_props%g)
        !$omp target teams distribute parallel do simd collapse(3) &
        !$omp map(from:optical_props%tau, optical_props%ssa, optical_props%g)
        do ibnd = 1, nbnd
          do ilay = 1, nlay
            do icol = 1, ncol
              tau    = atau   (icol,ilay,ibnd)
              taussa = ataussa(icol,ilay,ibnd)
              optical_props%g  (icol,ilay,ibnd) = (ataussag(icol,ilay,ibnd)) / &
                                                        max(epsilon(tau), taussa)
              optical_props%ssa(icol,ilay,ibnd) = taussa/max(epsilon(tau), tau)
              optical_props%tau(icol,ilay,ibnd) = tau
            end do
          end do
        end do
      type is (ty_optical_props_nstr)
        error_msg = "aerosol optics: n-stream calculations not yet supported"
      end select
    end if 
    !$acc end data
    !$omp end target data
  end function aerosol_optics
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function get_aerosol_size_bin(nbin, aero_bin_lo, aero_bin_hi, aero_size) result(ibin)
  !--------------------------------------------------------------------------------------------------------------------
  ! Purpose: For an input aerosol particle size, return the corresponding bin number
  !          from the MERRA aerosol optical property particle size dimension as defined
  !          by aero_bin_lo and aero_bin_hi. 
  !--------------------------------------------------------------------------------------------------------------------
    integer,                   intent(in   ) :: nbin
    real(wp), dimension(nbin), intent(in   ) :: aero_bin_lo, aero_bin_hi
    real(wp),                  intent(in   ) :: aero_size

    integer  :: i, ibin

    character(len=128)      :: error_msg

    ibin = 0

    error_msg = ''
    if(any_vals_outside(aero_size, aero_bin_lo(1), aero_bin_hi(nbin))) &
      error_msg = "get_aerosol_size_bin(): requested aerosol size outside of allowable range"
 
    do i=1,nbin 
       if (aero_size .ge. aero_bin_lo(i) .and. aero_size .le. aero_bin_hi(i)) then
          ibin = i
       endif
    enddo

  end function get_aerosol_size_bin
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! For size dimension, select size bin appropriate for the requested aerosol size.
  ! For rh dimension, linearly interpolate values from a lookup table with "nrh" 
  !   unevenly-spaced elements "aero_rh". The last dimension for all tables is band.
  ! Returns 0 where the mask is false.
  !
  subroutine compute_all_from_table(ncol, nlay, nbnd, nbin, nrh,                    &
                                    mask, type, size,  mass, rh,                    &
                                    aero_bin_lo, aero_bin_hi, aero_rh,              &
                                    rh_nsteps, rh_step_sizes, rh_offsets,           &
                                    aero_dust, aero_salt, aero_sulf,                &
                                    aero_bcar, aero_bcar_rh, aero_ocar, aero_ocar_rh, &
                                    tau_dust_table, ssa_dust_table, asy_dust_table, &
                                    tau_salt_table, ssa_salt_table, asy_salt_table, &
                                    tau_sulf_table, ssa_sulf_table, asy_sulf_table, &
                                    tau_bcar_rh_table, ssa_bcar_rh_table, asy_bcar_rh_table, &
                                    tau_bcar_table, ssa_bcar_table, asy_bcar_table, &
                                    tau_ocar_rh_table, ssa_ocar_rh_table, asy_ocar_rh_table, &
                                    tau_ocar_table, ssa_ocar_table, asy_ocar_table, &
                                    tau, taussa, taussag)

    integer,                               intent(in) :: ncol, nlay, nbnd, nbin, nrh
    logical(wl), dimension(ncol,nlay),     intent(in) :: mask
    integer,     dimension(ncol,nlay),     intent(in) :: type
    real(wp),    dimension(ncol,nlay),     intent(in) :: size, mass, rh

    real(wp),    dimension(nbin),          intent(in) :: aero_bin_lo, aero_bin_hi
    real(wp),    dimension(nrh),           intent(in) :: aero_rh
    integer,     dimension(2),             intent(in) :: rh_nsteps
    real(wp),    dimension(2),             intent(in) :: rh_step_sizes, rh_offsets

    integer,                               intent(in) :: aero_dust, aero_salt, aero_sulf
    integer,                               intent(in) :: aero_bcar, aero_bcar_rh, aero_ocar, aero_ocar_rh
 
    real(wp),    dimension(nbin,    nbnd), intent(in) :: tau_dust_table, ssa_dust_table, asy_dust_table
    real(wp),    dimension(nbin,nrh,nbnd), intent(in) :: tau_salt_table, ssa_salt_table, asy_salt_table
    real(wp),    dimension(     nrh,nbnd), intent(in) :: tau_sulf_table, ssa_sulf_table, asy_sulf_table
    real(wp),    dimension(     nrh,nbnd), intent(in) :: tau_bcar_rh_table, ssa_bcar_rh_table, asy_bcar_rh_table
    real(wp),    dimension(         nbnd), intent(in) :: tau_bcar_table, ssa_bcar_table, asy_bcar_table
    real(wp),    dimension(     nrh,nbnd), intent(in) :: tau_ocar_rh_table, ssa_ocar_rh_table, asy_ocar_rh_table
    real(wp),    dimension(         nbnd), intent(in) :: tau_ocar_table, ssa_ocar_table, asy_ocar_table

    real(wp),    dimension(ncol,nlay,nbnd)            :: tau, taussa, taussag
    ! ---------------------------
    integer  :: icol, ilay, ibnd, ibin
    integer  :: index
    real(wp) :: fint
    real(wp) :: t, ts, tsg  ! tau, tau*ssa, tau*ssa*g
    ! ---------------------------
    !$acc parallel loop gang vector default(present) collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do ibnd = 1, nbnd
      do ilay = 1,nlay
        do icol = 1, ncol
          ibin = get_aerosol_size_bin(nbin, aero_bin_lo, aero_bin_hi, size(icol,ilay))
          if (mask(icol,ilay) .and. rh(icol,ilay) .le. rh_offsets(2)) then
            index = min(floor((rh(icol,ilay) - rh_offsets(1))/rh_step_sizes(1))+1, rh_nsteps(1)-1)
            fint = (rh(icol,ilay) - rh_offsets(1))/rh_step_sizes(1) - (index-1)
          elseif (mask(icol,ilay) .and. rh(icol,ilay) .gt. rh_offsets(2)) then
            index = min(floor((rh(icol,ilay) - rh_offsets(2))/rh_step_sizes(2))+1, rh_nsteps(2)-1)
            fint = (rh(icol,ilay) - rh_offsets(2))/rh_step_sizes(2) - (index-1)
          endif
          ! dust
          if (mask(icol,ilay) .and. type(icol,ilay) == aero_dust) then
            t   = mass(icol,ilay) * tau_dust_table(ibin,ibnd)
            ts  = t               * ssa_dust_table(ibin,ibnd)
            taussag(icol,ilay,ibnd) =  &
                  ts              * asy_dust_table(ibin,ibnd)
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          ! sea-salt  (IN PROGRESS)
          elseif (mask(icol,ilay) .and. type(icol,ilay) == aero_salt) then
            t   = mass(icol,ilay) * &
                  (tau_salt_table(ibin,index,ibnd) + fint * (tau_salt_table(ibin,index+1,ibnd) - tau_salt_table(ibin,index,ibnd)))
            ts  = t              * &
                  (ssa_salt_table(ibin,index,ibnd) + fint * (ssa_salt_table(ibin,index+1,ibnd) - ssa_salt_table(ibin,index,ibnd)))
            taussag(icol,ilay,ibnd) =  &
                  ts             * &
                  (asy_salt_table(ibin,index,ibnd) + fint * (asy_salt_table(ibin,index+1,ibnd) - asy_salt_table(ibin,index,ibnd)))
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          ! sulfate  (IN PROGRESS)
          elseif (mask(icol,ilay) .and. type(icol,ilay) == aero_sulf) then
            t   = mass(icol,ilay) * &
                  (tau_sulf_table(index,  ibnd) + fint * (tau_sulf_table(index+1,ibnd) - tau_sulf_table(index,ibnd)))
            ts  = t              * &
                  (ssa_sulf_table(index,  ibnd) + fint * (ssa_sulf_table(index+1,ibnd) - ssa_sulf_table(index,ibnd)))
            taussag(icol,ilay,ibnd) =  &
                  ts             * &
                  (asy_sulf_table(index,  ibnd) + fint * (asy_sulf_table(index+1,ibnd) - asy_sulf_table(index,ibnd)))
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          ! black carbon - hydrophilic  (IN PROGRESS)
          elseif (mask(icol,ilay) .and. type(icol,ilay) == aero_bcar_rh) then
            t   = mass(icol,ilay) * &
                  (tau_bcar_rh_table(index,  ibnd) + fint * (tau_bcar_rh_table(index+1,ibnd) - tau_bcar_rh_table(index,ibnd)))
            ts  = t              * &
                  (ssa_bcar_rh_table(index,  ibnd) + fint * (ssa_bcar_rh_table(index+1,ibnd) - ssa_bcar_rh_table(index,ibnd)))
            taussag(icol,ilay,ibnd) =  &
                  ts             * &
                  (asy_bcar_rh_table(index,  ibnd) + fint * (asy_bcar_rh_table(index+1,ibnd) - asy_bcar_rh_table(index,ibnd)))
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          ! black carbon - hydrophobic
          elseif (mask(icol,ilay) .and. type(icol,ilay) == aero_bcar) then
            t   = mass(icol,ilay) * tau_bcar_table(ibnd)
            ts  = t               * ssa_bcar_table(ibnd)
            taussag(icol,ilay,ibnd) =  &
                  ts              * asy_bcar_table(ibnd)
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          ! organic carbon - hydrophilic  (IN PROGRESS)
          elseif (mask(icol,ilay) .and. type(icol,ilay) == aero_ocar_rh) then
            t   = mass(icol,ilay) * &
                  (tau_ocar_rh_table(index,  ibnd) + fint * (tau_ocar_rh_table(index+1,ibnd) - tau_ocar_rh_table(index,ibnd)))
            ts  = t              * &
                  (ssa_ocar_rh_table(index,  ibnd) + fint * (ssa_ocar_rh_table(index+1,ibnd) - ssa_ocar_rh_table(index,ibnd)))
            taussag(icol,ilay,ibnd) =  &
                  ts             * &
                  (asy_ocar_rh_table(index,  ibnd) + fint * (asy_ocar_rh_table(index+1,ibnd) - asy_ocar_rh_table(index,ibnd)))
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          ! organic carbon - hydrophobic
          elseif (mask(icol,ilay) .and. type(icol,ilay) == aero_ocar) then
            t   = mass(icol,ilay) * tau_ocar_table(ibnd)
            ts  = t               * ssa_ocar_table(ibnd)
            taussag(icol,ilay,ibnd) =  &
                  ts              * asy_ocar_table(ibnd)
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          else
            tau    (icol,ilay,ibnd) = 0._wp
            taussa (icol,ilay,ibnd) = 0._wp
            taussag(icol,ilay,ibnd) = 0._wp
          end if
        end do
      end do
    end do
  end subroutine compute_all_from_table

end module mo_aerosol_optics
