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
!   References for the gocart interactive aerosols:   
!     Chin et al., jgr, 2000 (https://doi.org/10.1029/2000jd900384)
!     Chin et al., jas, 2002 (https://doi.org/10.1175/1520-0469(2002)059<0461:TAOTFT>2.0.CO;2)
!     Colarco et al., jgr, 2010 (https://doi.org/10.1029/2009jd012820)
!                                                                      
!   References for merra2 aerosol reanalysis:                          
!     Randles et al., j. clim., 2017 (https://doi.org/10.1175/jcli-d-16-0609.1) 
!     Buchard et al., j. clim., 2017 (https://doi.org/10.1175/jcli-d-16-0613.1)
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

  ! MERRA2/GOCART aerosol types
  integer, parameter, public :: merra_ntype = 7          ! Number of MERRA aerosol types
  integer, parameter, public :: merra_aero_dust = 1      ! Dust
  integer, parameter, public :: merra_aero_salt = 2      ! Salt
  integer, parameter, public :: merra_aero_sulf = 3      ! sulfate
  integer, parameter, public :: merra_aero_bcar_rh = 4   ! black carbon, hydrophilic
  integer, parameter, public :: merra_aero_bcar = 5      ! black carbon, hydrophobic
  integer, parameter, public :: merra_aero_ocar_rh = 6   ! organic carbon, hydrophilic
  integer, parameter, public :: merra_aero_ocar = 7      ! organic carbon, hydrophobic

  private
  ! -----------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_aerosol_optics
    private
    !
    ! Lookup table information
    !
    ! Table upper and lower size bin limits and effective (center) values
    real(wp),dimension(:,:), allocatable :: merra_aero_bin_lims     ! Dimensions (pair,nbin)
    ! Table relative humidity values
    real(wp),dimension(:), allocatable :: aero_rh(:)
    !
    ! The aerosol tables themselves.
    ! extinction (m2/kg)
    ! single scattering albedo (unitless)
    ! asymmetry parameter (unitless)
    !
    real(wp), dimension(:,:,:  ), allocatable :: aero_dust_tbl      ! ext, ssa, g (nval, nbin, nbnd)
    real(wp), dimension(:,:,:,:), allocatable :: aero_salt_tbl      ! ext, ssa, g (nval, nrh, nbin, nbnd)
    real(wp), dimension(:,:,:  ), allocatable :: aero_sulf_tbl      ! ext, ssa, g (nval, nrh, nbnd)
    real(wp), dimension(:,:    ), allocatable :: aero_bcar_tbl      ! ext, ssa, g (nval, nbnd)
    real(wp), dimension(:,:,:  ), allocatable :: aero_bcar_rh_tbl   ! ext, ssa, g (nval, nrh, nbnd)
    real(wp), dimension(:,:    ), allocatable :: aero_ocar_tbl      ! ext, ssa, g (nval, nbnd)
    real(wp), dimension(:,:,:  ), allocatable :: aero_ocar_rh_tbl   ! ext, ssa, g (nval, nrh, nbnd)
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
                    merra_aero_bin_lims, aero_rh, &
                    aero_dust_tbl, aero_salt_tbl, aero_sulf_tbl, &
                    aero_bcar_tbl, aero_bcar_rh_tbl, &
                    aero_ocar_tbl, aero_ocar_rh_tbl) &
                    result(error_msg)

    class(ty_aerosol_optics),   intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn ! spectral discretization
    ! Lookup table interpolation constants
    real(wp), dimension(:,:),   intent(in   ) :: merra_aero_bin_lims ! aerosol lut size bin limiits (pair,nbin)
    real(wp), dimension(:),     intent(in   ) :: aero_rh       ! relative humidity LUT dimension values
    ! LUT coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter for aerosol types
    real(wp), dimension(:,:,:),   intent(in)  :: aero_dust_tbl
    real(wp), dimension(:,:,:,:), intent(in)  :: aero_salt_tbl
    real(wp), dimension(:,:,:),   intent(in)  :: aero_sulf_tbl
    real(wp), dimension(:,:),     intent(in)  :: aero_bcar_tbl
    real(wp), dimension(:,:,:),   intent(in)  :: aero_bcar_rh_tbl
    real(wp), dimension(:,:),     intent(in)  :: aero_ocar_tbl
    real(wp), dimension(:,:,:),   intent(in)  :: aero_ocar_rh_tbl
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer               :: npair, nbin, nrh, nband

    error_msg = this%init(band_lims_wvn, name="RRTMGP aerosol optics")
    !
    ! LUT coefficient dimensions
    !
    npair = size(merra_aero_bin_lims,dim=1)
    nval  = size(aero_salt_tbl,dim=1)
    nrh   = size(aero_salt_tbl,dim=2)
    nbin  = size(aero_salt_tbl,dim=3)
    nband = size(aero_salt_tbl,dim=4)
    !
    ! Error checking
    !
    if (check_extents) then
!    if(nband /= this%get_nband()) &
!      error_msg = "aerosol_optics%init(): number of bands inconsistent between lookup tables, spectral discretization"
      error_msg = ''
      if(.not. extents_are(aero_dust_tbl, nval, nbin, nband)) &
        error_msg = "aerosol_optics%init(): array aero_dust_tbl isn't consistently sized"
      if(.not. extents_are(aero_salt_tbl, nval, nrh, nbin, nband)) &
        error_msg = "aerosol_optics%init(): array aero_salt_tbl isn't consistently sized"
      if(.not. extents_are(aero_sulf_tbl, nval, nrh, nband)) &
        error_msg = "aerosol_optics%init(): array aero_sulf_tbl isn't consistently sized"
      if(.not. extents_are(aero_bcar_rh_tbl, nval, nrh, nband)) &
        error_msg = "aerosol_optics%init(): array aero_bcar_rh_tbl isn't consistently sized"
      if(.not. extents_are(aero_bcar_tbl, nval, nband)) &
        error_msg = "aerosol_optics%init(): array aero_bcar_tbl isn't consistently sized"
      if(.not. extents_are(aero_ocar_rh_tbl, nval, nrh, nband)) &
        error_msg = "aerosol_optics%init(): array aero_ocar_rh_tbl isn't consistently sized"
      if(.not. extents_are(aero_ocar_tbl, nval, nband)) &
        error_msg = "aerosol_optics%init(): array aero_ocar_tbl isn't consistently sized"
      if(error_msg /= "") return
    endif

    ! Allocate LUT parameters
    allocate(this%merra_aero_bin_lims(npair,nbin))
    this%merra_aero_bin_lims = merra_aero_bin_lims

    allocate(this%aero_rh(nrh))
    this%aero_rh = aero_rh

    ! Allocate LUT coefficients
    allocate(this%aero_dust_tbl(nval, nbin, nband), &
             this%aero_salt_tbl(nval, nrh, nbin, nband), &
             this%aero_sulf_tbl(nval, nrh, nband), &
             this%aero_bcar_tbl(nval, nband), &
             this%aero_bcar_rh_tbl(nval, nrh, nband), &
             this%aero_ocar_tbl(nval, nband), &
             this%aero_ocar_rh_tbl(nval, nrh, nband))

    ! Load LUT coefficients
    !$acc kernels
    !$omp target
    this%aero_dust_tbl = aero_dust_tbl
    this%aero_salt_tbl = aero_salt_tbl
    this%aero_sulf_tbl = aero_sulf_tbl
    this%aero_bcar_tbl = aero_bcar_tbl
    this%aero_bcar_rh_tbl = aero_bcar_rh_tbl
    this%aero_ocar_tbl = aero_ocar_tbl
    this%aero_ocar_rh_tbl = aero_ocar_rh_tbl
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
    if(allocated(this%merra_aero_bin_lims)) then

      deallocate(this%merra_aero_bin_lims, this%aero_rh)
    end if

    ! Lookup table aerosol optics coefficients
    if(allocated(this%aero_dust_tbl)) then

      deallocate(this%aero_dust_tbl, this%aero_salt_tbl, this%aero_sulf_tbl, &
                 this%aero_bcar_tbl, this%aero_bcar_rh_tbl, &
                 this%aero_ocar_tbl, this%aero_ocar_rh_tbl)
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
                                               ! 1 = merra_aero_dust    (dust)
                                               ! 2 = merra_aero_salt    (salt)
                                               ! 3 = merra_aero_sulf    (sulfate)
                                               ! 4 = merra_aero_bcar_rh (black carbon, hydrophilic)
                                               ! 5 = merra_aero_bcar    (black carbon, hydrophobic)
                                               ! 6 = merra_aero_ocar_rh (organic carbon, hydrophilic)
                                               ! 7 = merra_aero_ocar    (organic carbon, hydrophobic)
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
    real(wp),    dimension(size(aero_type,1), size(aero_type,2), size(this%aero_dust_tbl,3)) :: &
                 atau, ataussa, ataussag
    integer  :: ncol, nlay, nbnd, npair, nbin, nrh
    integer  :: icol, ilay, ibnd, ibin
    ! scalars for total tau, tau*ssa
    real(wp) :: tau, taussa

    ! ----------------------------------------
    !
    ! Error checking
    !
    ! ----------------------------------------

    error_msg = ''
    if(.not.(allocated(this%aero_dust_tbl))) then
      error_msg = 'aerosol optics: no data has been initialized'
      return
    end if

    ncol = size(aero_type,1)
    nlay = size(aero_type,2)
    npair= size(this%merra_aero_bin_lims,1)
    nbin = size(this%merra_aero_bin_lims,2)
    nrh  = size(this%aero_rh,1)
    nval = size(this%aero_dust_tbl,1)
    nbnd = size(this%aero_dust_tbl,3)

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
      if(any_vals_outside(aero_type, 0, merra_ntype)) &
        error_msg = 'aerosol optics: aerosol type is out of bounds'
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
      if(any_vals_outside(aero_size, aeromsk, &
        this%merra_aero_bin_lims(0,1), this%merra_aero_bin_lims(1,nbin))) &
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
      if (allocated(this%aero_dust_tbl)) then
        !
        ! Aerosol
        !
        call compute_all_from_table(ncol, nlay, nbnd, npair, nbin, nrh,                         &
                                    aero_type, aero_size, aero_mass, relhum,                    &
                                    this%merra_aero_bin_lims, this%aero_rh,                     &
                                    this%aero_dust_tbl, this%aero_salt_tbl, this%aero_sulf_tbl, &
                                    this%aero_bcar_rh_tbl, this%aero_bcar_tbl,                  &
                                    this%aero_ocar_rh_tbl, this%aero_ocar_tbl,                  &
                                    atau, ataussa, ataussag)

      endif

      !
      ! Derive total aerosol optical properties
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
  function get_aerosol_size_bin(npair, nbin, merra_aero_bin_lims, aero_size) result(ibin)
  !--------------------------------------------------------------------------------------------------------------------
  ! Purpose: For an input aerosol particle size, return the corresponding bin number
  !          from the MERRA aerosol optical property particle size dimension as defined
  !          by merra_aero_bin_lims. 
  !--------------------------------------------------------------------------------------------------------------------
    integer,                          intent(in   ) :: npair, nbin
    real(wp), dimension(npair, nbin), intent(in   ) :: merra_aero_bin_lims
    real(wp),                         intent(in   ) :: aero_size

    integer  :: i, ibin

    character(len=128)      :: error_msg

    ibin = 0

    error_msg = ''
    if(any_vals_outside(aero_size, merra_aero_bin_lims(0,1), merra_aero_bin_lims(1,nbin))) &
      error_msg = "get_aerosol_size_bin(): requested aerosol size outside of allowable range"
 
    do i=1,nbin 
       if (aero_size .ge. merra_aero_bin_lims(0,i) .and. &
           aero_size .le. merra_aero_bin_lims(1,i)) then
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
  ! Returns zero where no aerosol is present.
  !
  subroutine compute_all_from_table(ncol, nlay, nbnd, npair, nbin, nrh,                &
                                    type, size, mass, rh,                              &
                                    merra_aero_bin_lims, aero_rh,                      &
                                    aero_dust_table, aero_salt_table, aero_sulf_table, &
                                    aero_bcar_rh_table, aero_bcar_table,               &
                                    aero_ocar_rh_table, aero_ocar_table,               &
                                    tau, taussa, taussag)

    integer,                               intent(in) :: ncol, nlay, nbnd, npair, nbin, nrh
    integer,     dimension(ncol,nlay),     intent(in) :: type
    real(wp),    dimension(ncol,nlay),     intent(in) :: size, mass, rh

    real(wp),    dimension(npair,nbin),    intent(in) :: merra_aero_bin_lims
    real(wp),    dimension(nrh),           intent(in) :: aero_rh

    real(wp),    dimension(nval,    nbin,nbnd), intent(in) :: aero_dust_table
    real(wp),    dimension(nval,nrh,nbin,nbnd), intent(in) :: aero_salt_table
    real(wp),    dimension(nval,nrh,     nbnd), intent(in) :: aero_sulf_table
    real(wp),    dimension(nval,nrh,     nbnd), intent(in) :: aero_bcar_rh_table
    real(wp),    dimension(nval,         nbnd), intent(in) :: aero_bcar_table
    real(wp),    dimension(nval,nrh,     nbnd), intent(in) :: aero_ocar_rh_table
    real(wp),    dimension(nval,         nbnd), intent(in) :: aero_ocar_table

    real(wp),    dimension(ncol,nlay,nbnd)            :: tau, taussa, taussag
    ! ---------------------------
    integer  :: icol, ilay, ibnd, ibin
    integer  :: itype, irh1, irh2
    real(wp) :: drh0, drh1, rdrh
    real(wp) :: t, ts, tsg  ! tau, tau*ssa, tau*ssa*g
    ! ---------------------------
    !$acc parallel loop gang vector default(present) collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do ibnd = 1, nbnd
      do ilay = 1,nlay
        do icol = 1, ncol
          ibin = get_aerosol_size_bin(nbin, merra_aero_bin_lims, size(icol,ilay))
          itype = type(icol,ilay)
          ! relative humidity linear interpolation coefficients
          if (itype .ne. 0) then
             irh2 = 1
             do while (rh(icol,ilay) .gt. aero_rh(irh2))
                irh2 = irh2 + 1
                if (irh2 .gt. nrh) exit
             enddo
             irh1 = max(1, irh2-1)
             irh2 = min(nrh, irh2)
             drh0 = aero_rh(irh2) - aero_rh(irh1)
             drh1 = rh(icol,ilay) - aero_rh(irh1)
             if (irh1 == irh2) then
                rdrh = 0._wp
             else
                rdrh = drh1 / drh0
             endif
          endif

          ! Set aerosol optical properties where aerosol present. Use aerosol type array as the mask. 
          select case (itype)

             ! dust
             case(merra_aero_dust)
               t   = mass(icol,ilay) * aero_dust_table(1,ibin,ibnd)
               ts  = t               * aero_dust_table(2,ibin,ibnd)
               taussag(icol,ilay,ibnd) =  &
                     ts              * aero_dust_table(3,ibin,ibnd)
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! sea-salt
             case(merra_aero_salt)
               t   = mass(icol,ilay) * &
                     (aero_salt_table(1,irh1,ibin,ibnd) + rdrh * (aero_salt_table(1,irh2,ibin,ibnd) - aero_salt_table(1,irh1,ibin,ibnd)))
               ts  = t              * &
                     (aero_salt_table(2,irh1,ibin,ibnd) + rdrh * (aero_salt_table(2,irh2,ibin,ibnd) - aero_salt_table(2,irh1,ibin,ibnd)))
               taussag(icol,ilay,ibnd) =  &
                     ts             * &
                     (aero_salt_table(3,irh1,ibin,ibnd) + rdrh * (aero_salt_table(3,irh2,ibin,ibnd) - aero_salt_table(3,irh1,ibin,ibnd)))
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! sulfate
             case(merra_aero_sulf)
               t   = mass(icol,ilay) * &
                     (aero_sulf_table(1,irh1,ibnd) + rdrh * (aero_sulf_table(1,irh2,ibnd) - aero_sulf_table(1,irh1,ibnd)))
               ts  = t              * &
                     (aero_sulf_table(2,irh1,ibnd) + rdrh * (aero_sulf_table(2,irh2,ibnd) - aero_sulf_table(2,irh1,ibnd)))
               taussag(icol,ilay,ibnd) =  &
                     ts             * &
                     (aero_sulf_table(3,irh1,ibnd) + rdrh * (aero_sulf_table(3,irh2,ibnd) - aero_sulf_table(3,irh1,ibnd)))
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! black carbon - hydrophilic
             case(merra_aero_bcar_rh)
               t   = mass(icol,ilay) * &
                     (aero_bcar_rh_table(1,irh1,ibnd) + rdrh * (aero_bcar_rh_table(1,irh2,ibnd) - aero_bcar_rh_table(1,irh1,ibnd)))
               ts  = t              * &
                     (aero_bcar_rh_table(2,irh1,ibnd) + rdrh * (aero_bcar_rh_table(2,irh2,ibnd) - aero_bcar_rh_table(2,irh1,ibnd)))
               taussag(icol,ilay,ibnd) =  &
                     ts             * &
                     (aero_bcar_rh_table(3,irh1,ibnd) + rdrh * (aero_bcar_rh_table(3,irh2,ibnd) - aero_bcar_rh_table(3,irh1,ibnd)))
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! black carbon - hydrophobic
             case(merra_aero_bcar)
               t   = mass(icol,ilay) * aero_bcar_table(1,ibnd)
               ts  = t               * aero_bcar_table(2,ibnd)
               taussag(icol,ilay,ibnd) =  &
                     ts              * aero_bcar_table(3,ibnd)
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! organic carbon - hydrophilic
             case(merra_aero_ocar_rh)
               t   = mass(icol,ilay) * &
                     (aero_ocar_rh_table(1,irh1,ibnd) + rdrh * (aero_ocar_rh_table(1,irh2,ibnd) - aero_ocar_rh_table(1,irh1,ibnd)))
               ts  = t              * &
                     (aero_ocar_rh_table(2,irh1,ibnd) + rdrh * (aero_ocar_rh_table(2,irh2,ibnd) - aero_ocar_rh_table(2,irh1,ibnd)))
               taussag(icol,ilay,ibnd) =  &
                     ts             * &
                     (aero_ocar_rh_table(3,irh1,ibnd) + rdrh * (aero_ocar_rh_table(3,irh2,ibnd) - aero_ocar_rh_table(3,irh1,ibnd)))
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! organic carbon - hydrophobic
             case(merra_aero_ocar)
               t   = mass(icol,ilay) * aero_ocar_table(1,ibnd)
               ts  = t               * aero_ocar_table(2,ibnd)
               taussag(icol,ilay,ibnd) =  &
                     ts              * aero_ocar_table(3,ibnd)
               taussa (icol,ilay,ibnd) = ts
               tau    (icol,ilay,ibnd) = t
             ! no aerosol
             case default
               tau    (icol,ilay,ibnd) = 0._wp
               taussa (icol,ilay,ibnd) = 0._wp
               taussag(icol,ilay,ibnd) = 0._wp

          end select

        end do
      end do
    end do
  end subroutine compute_all_from_table
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Function for interpolating MERRA 2-d aerosol optics tables in the rh dimension
  ! for each spectral band, vertical layer and column. Interpolation is performed
  ! only where aerosol in present using aerosol type as the mask. 
  ! Returns zero where no aerosol is present.
  !
  function interp_aero_table_2d()

  end function interp_aero_table_2d

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Function for interpolating MERRA 3-d aerosol optics tables in the rh dimension
  ! for each spectral band, vertical layer and column. Interpolation is performed
  ! only where aerosol in present using aerosol type as the mask. 
  ! Returns zero where no aerosol is present.
  !
  function interp_aero_table_3d()

  end function interp_aero_table_3d

end module mo_aerosol_optics
