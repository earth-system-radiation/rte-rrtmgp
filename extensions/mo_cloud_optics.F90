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
! Provides cloud optical properties as a function of effective radius for the RRTMGP bands
!   Based on Mie calculations for liquid and Ping Yang lookup-tables for ice
!   Can use either look-up tables or Pade approximates according to the do_lut flag
!   Mike Iacono is the original author
! -------------------------------------------------------------------------------------------------

module mo_cloud_optics
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  implicit none
  private

  ! -----------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_cloud_optics
    ! Cloud physical properties                         ! (ncol,nlay)
    real(wp), dimension(:,:), allocatable :: cldfrac    ! cloud fraction
    real(wp), dimension(:,:), allocatable :: ciwp       ! cloud ice water path
    real(wp), dimension(:,:), allocatable :: clwp       ! cloud liquid water path
    real(wp), dimension(:,:), allocatable :: rei        ! cloud ice particle effective size (microns)
    real(wp), dimension(:,:), allocatable :: rel        ! cloud liquid particle effective radius (microns)

    ! All other data should be privide
    !
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    integer :: icergh                                   ! (1 = none, 2 = medium, 3 = high)

    ! Method for interpolation of cloud optical property coefficients to particle size
    logical :: do_lut                                   ! (.True. = LUT, .False. = Pade)

    ! Particle size boundary limits
    real(wp) :: radliq_lwr                     ! liquid particle size lower bound for interpolation
    real(wp) :: radliq_upr                     ! liquid particle size upper bound for interpolation
    real(wp) :: radice_lwr                     ! ice particle size lower bound for interpolation
    real(wp) :: radice_upr                     ! ice particle size upper bound for interpolation
    ! Lookup table interpolation constants
    real(wp) :: radliq_fac                     ! constant for calculating interpolation indices for liquid
    real(wp) :: radice_fac                     ! constant for calculating interpolation indices for ice
    ! Lookup table cloud optics coefficients
    real(wp), dimension(:,:    ), allocatable :: lut_extliq     ! (nsize_liq, nband)
    real(wp), dimension(:,:    ), allocatable :: lut_ssaliq     ! (nsize_liq, nband)
    real(wp), dimension(:,:    ), allocatable :: lut_asyliq     ! (nsize_liq, nband)
    real(wp), dimension(:,:,:  ), allocatable :: lut_extice     ! (nsize_ice, nband, nrghice)
    real(wp), dimension(:,:,:  ), allocatable :: lut_ssaice     ! (nsize_ice, nband, nrghice)
    real(wp), dimension(:,:,:  ), allocatable :: lut_asyice     ! (nsize_ice, nband, nrghice)

    ! Pade cloud optics coefficients
    real(wp), dimension(:,:,:  ), allocatable :: pade_extliq    ! (nband, nsizereg, ncoeff_ext)
    real(wp), dimension(:,:,:  ), allocatable :: pade_ssaliq    ! (nband, nsizereg, ncoeff_ssa_g)
    real(wp), dimension(:,:,:  ), allocatable :: pade_asyliq    ! (nband, nsizereg, ncoeff_ssa_g)
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice    ! (nband, nsizereg, ncoeff_ext, nrghice)
    real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice    ! (nband, nsizereg, ncoeff_ssa_g, nrghice)
    real(wp), dimension(:,:,:,:), allocatable :: pade_asyice    ! (nband, nsizereg, ncoeff_ssa_g, nrghice)
    ! Particle size regimes for Pade formulations
    integer, dimension(:), allocatable :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq  ! (nbound)
    integer, dimension(:), allocatable :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice  ! (nbound)

! ------------------------------------------------------------------------------------------
  contains
    generic,   public :: load  => load_lut, load_pade
    procedure, public :: cloud_optics
    ! Internal procedures
    procedure, private :: load_lut
    procedure, private :: load_pade
    procedure, private :: pade_ext
    procedure, private :: pade_ssa
    procedure, private :: pade_asy
  end type ty_cloud_optics

contains
  ! ------------------------------------------------------------------------------
  !
  ! Routines to load data needed for cloud optics calculations. Two routines: one to load
  !    lookup-tables and one for coefficients for Pade approximates
  !
  ! ------------------------------------------------------------------------------
  function load_lut(this, band_lims_wvn, &
                    radliq_lwr, radliq_upr, radliq_fac, &
                    radice_lwr, radice_upr, radice_fac, &
                    lut_extliq, lut_ssaliq, lut_asyliq, &
                    lut_extice, lut_ssaice, lut_asyice) result(error_msg)
    class(ty_cloud_optics),     intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn ! Spectral discretization
    ! Lookup table interpolation constants
    ! Lower and upper bounds of the tables; also the constant for calculating interpolation indices for liquid
    real(wp),                   intent(in   ) :: radliq_lwr, radliq_upr, radliq_fac
    real(wp),                   intent(in   ) :: radice_lwr, radice_upr, radice_fac
    ! LUT coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter for liquid and ice respectively
    real(wp), dimension(:,:),   intent(in)    :: lut_extliq, lut_ssaliq, lut_asyliq
    real(wp), dimension(:,:,:), intent(in)    :: lut_extice, lut_ssaice, lut_asyice
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer               :: nband, nrghice, nsize_liq, nsize_ice

    error_msg = this%init(band_lims_wvn, name="RRTMGP cloud optics")
    !
    ! LUT coefficient dimensions
    !
    nsize_liq = size(lut_extliq,dim=1)
    nsize_ice = size(lut_extice,dim=1)
    nband     = size(lut_extliq,dim=2)
    nrghice   = size(lut_extice,dim=3)
    !
    ! Error checking
    !   Can we check for consistency between table bounds and _fac?
    !
    if(nband /= this%get_nband()) &
      error_msg = "cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization"
    if(size(lut_extice, 2) /= nband) &
      error_msg = "cloud_optics%init(): array lut_extice has the wrong number of bands"
    if(any([size(lut_ssaliq, 1), size(lut_ssaliq, 2)] /= [nsize_liq, nband])) &
      error_msg = "cloud_optics%init(): array lut_ssaliq isn't consistently sized"
    if(any([size(lut_asyliq, 1), size(lut_asyliq, 2)] /= [nsize_liq, nband])) &
      error_msg = "cloud_optics%init(): array ssaliq isn't consistently sized"
    if(any([size(lut_ssaice, 1), size(lut_ssaice, 2), size(lut_ssaice, 3)] /= [nsize_ice, nband, nrghice])) &
      error_msg = "cloud_optics%init(): array ssaice isn't consistently sized"
    if(any([size(lut_asyice, 1), size(lut_asyice, 2), size(lut_asyice, 3)] /= [nsize_ice, nband, nrghice])) &
      error_msg = "cloud_optics%init(): array ssaice isn't consistently sized"
    if(error_msg /= "") return

    ! Allocate LUT coefficients
    allocate(this%lut_extliq(nsize_liq, nband), &
             this%lut_ssaliq(nsize_liq, nband), &
             this%lut_asyliq(nsize_liq, nband), &
             this%lut_extice(nsize_ice, nband, nrghice), &
             this%lut_ssaice(nsize_ice, nband, nrghice), &
             this%lut_asyice(nsize_ice, nband, nrghice))

    ! Load LUT constants
    this%radliq_lwr = radliq_lwr
    this%radliq_upr = radliq_upr
    this%radliq_fac = radliq_fac
    this%radice_lwr = radice_lwr
    this%radice_upr = radice_upr
    this%radice_fac = radice_fac

    ! Load LUT coefficients
    this%lut_extliq = lut_extliq
    this%lut_ssaliq = lut_ssaliq
    this%lut_asyliq = lut_asyliq
    this%lut_extice = lut_extice
    this%lut_ssaice = lut_ssaice
    this%lut_asyice = lut_asyice
  end function load_lut
  ! ------------------------------------------------------------------------------
  !
  ! Cloud optics initialization function - Pade
  !
  ! ------------------------------------------------------------------------------
  function load_pade(this, band_lims_wvn, &
                     radliq_lwr, radliq_upr, radice_lwr, radice_upr, &
                     pade_extliq, pade_ssaliq, pade_asyliq, &
                     pade_extice, pade_ssaice, pade_asyice, &
                     pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
                     pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice) &
                     result(error_msg)
    class(ty_cloud_optics),       intent(inout) :: this          ! cloud specification data
    real(wp), dimension(:,:),     intent(in   ) :: band_lims_wvn ! Spectral discretization
    !
    ! Particle size boundary limits
    !
    real(wp),                     intent(in   ) :: radliq_lwr, radliq_upr
    real(wp),                     intent(in   ) :: radice_lwr, radice_upr
    !
    ! Pade coefficients: extinction, single-scattering albedo, and asymmetry factor for liquid and ice
    !
    real(wp), dimension(:,:,:),   intent(in)    :: pade_extliq, pade_ssaliq, pade_asyliq
    real(wp), dimension(:,:,:,:), intent(in)    :: pade_extice, pade_ssaice, pade_asyice
    !
    ! Boundaries of size regimes. Liquid and ice are separate;
    !   extinction is fit to different numbers of size bins than single-scattering albedo and asymmetry factor
    !
    integer,  dimension(:),       intent(in)    :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq
    integer,  dimension(:),       intent(in)    :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice
    character(len=128)    :: error_msg

! ------- Local -------

    integer               :: nband, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

! ------- Definitions -------

    ! Pade coefficient dimensions
    nband        = size(pade_extliq,dim=1)
    nsizereg     = size(pade_extliq,dim=2)
    ncoeff_ext   = size(pade_extliq,dim=3)
    ncoeff_ssa_g = size(pade_ssaliq,dim=3)
    nrghice      = size(pade_extice,dim=4)
    nbound       = size(pade_sizreg_extliq,dim=1)
    error_msg = this%init(band_lims_wvn, name="RRTMGP cloud optics")
    !
    ! Error checking
    !
    if(nband /= this%get_nband()) &
      error_msg = "cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization"
    if(any([size(pade_ssaliq, 1), size(pade_ssaliq, 2), size(pade_ssaliq, 3)] /= [nband, nsizereg, ncoeff_ssa_g])) &
      error_msg = "cloud_optics%init(): array ssaliq isn't consistently sized"
    if(any([size(pade_asyliq, 1), size(pade_asyliq, 2), size(pade_asyliq, 3)] /= [nband, nsizereg, ncoeff_ssa_g])) &
      error_msg = "cloud_optics%init(): array pade_asyliq isn't consistently sized"
    if(any([size(pade_extice, 1), size(pade_extice, 2), size(pade_extice, 3), size(pade_extice, 4)] /= &
           [nband,                nsizereg,             ncoeff_ext,           nrghice]))               &
      error_msg = "cloud_optics%init(): array pade_extice isn't consistently sized"
    if(any([size(pade_ssaice, 1), size(pade_ssaice, 2), size(pade_ssaice, 3), size(pade_ssaice, 4)] /= &
           [nband,                nsizereg,             ncoeff_ssa_g,         nrghice]))               &
      error_msg = "cloud_optics%init(): array pade_ssaice isn't consistently sized"
    if(any([size(pade_asyice, 1), size(pade_asyice, 2), size(pade_asyice, 3), size(pade_asyice, 4)] /= &
           [nband,                nsizereg,             ncoeff_ssa_g,         nrghice]))               &
      error_msg = "cloud_optics%init(): array pade_asyice isn't consistently sized"
    if(any([                          size(pade_sizreg_ssaliq), size(pade_sizreg_asyliq),               &
            size(pade_sizreg_extice), size(pade_sizreg_ssaice), size(pade_sizreg_asyice)] /= nbound))   &
      error_msg = "cloud_optics%init(): one or more Pade size regime arrays are inconsistently sized"
    if(error_msg /= "") return
    !
    ! Consistency between size regimes and lower bounds
    !
    if(.false) then  ! Commented out for now because the files we're reading are inconsistent in this way
    if(any([pade_sizreg_extliq(1), pade_sizreg_ssaliq(1), pade_sizreg_asyliq(1)] < radliq_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_lwr"
    if(any([pade_sizreg_extice(1), pade_sizreg_ssaice(1), pade_sizreg_asyice(1)] < radice_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_lwr"
    if(any([pade_sizreg_extliq(nbound), pade_sizreg_ssaliq(nbound), pade_sizreg_asyliq(nbound)] > radliq_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_upr"
    if(any([pade_sizreg_extice(nbound), pade_sizreg_ssaice(nbound), pade_sizreg_asyice(nbound)] > radice_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_upr"
    if(error_msg /= "") return
    end if

    ! Load particle size boundaries
    this%radliq_lwr = radliq_lwr
    this%radliq_upr = radliq_upr
    this%radice_lwr = radice_lwr
    this%radice_upr = radice_upr

    ! Allocate Pade coefficients
    allocate(this%pade_extliq(nband, nsizereg, ncoeff_ext), &
             this%pade_ssaliq(nband, nsizereg, ncoeff_ssa_g), &
             this%pade_asyliq(nband, nsizereg, ncoeff_ssa_g), &
             this%pade_extice(nband, nsizereg, ncoeff_ext,   nrghice), &
             this%pade_ssaice(nband, nsizereg, ncoeff_ssa_g, nrghice), &
             this%pade_asyice(nband, nsizereg, ncoeff_ssa_g, nrghice))

    ! Load Pade coefficients
    this%pade_extliq = pade_extliq
    this%pade_ssaliq = pade_ssaliq
    this%pade_asyliq = pade_asyliq
    this%pade_extice = pade_extice
    this%pade_ssaice = pade_ssaice
    this%pade_asyice = pade_asyice

    ! Allocate Pade coefficient particle size regime boundaries
    allocate(this%pade_sizreg_extliq(nbound), &
             this%pade_sizreg_ssaliq(nbound), &
             this%pade_sizreg_asyliq(nbound), &
             this%pade_sizreg_extice(nbound), &
             this%pade_sizreg_ssaice(nbound), &
             this%pade_sizreg_asyice(nbound))

    ! Load Pade coefficient particle size regime boundaries
    this%pade_sizreg_extliq = pade_sizreg_extliq
    this%pade_sizreg_ssaliq = pade_sizreg_ssaliq
    this%pade_sizreg_asyliq = pade_sizreg_asyliq
    this%pade_sizreg_extice = pade_sizreg_extice
    this%pade_sizreg_ssaice = pade_sizreg_ssaice
    this%pade_sizreg_asyice = pade_sizreg_asyice
  end function load_pade

  ! ------------------------------------------------------------------------------
  !
  ! Derive cloud optical properties from provided cloud physical properties
  !
  ! ------------------------------------------------------------------------------
  function cloud_optics(cloud_spec, &
                        ncol, nlayers, nbnd, nrghice, &
                        cldfrac, clwp, ciwp, rel, rei, optical_props) result(error_msg)
  ! ------------------------------------------------------------------------------

  ! Purpose:  Compute the cloud optical properties for each cloudy layer.

  ! ------- Input -------

    class(ty_cloud_optics),             intent(inout) :: cloud_spec
                                               ! cloud specification data
    integer, intent(in) :: ncol                ! total number of columns
    integer, intent(in) :: nlayers             ! total number of layers
    integer, intent(in) :: nbnd                ! number of bands
    integer, intent(in) :: nrghice             ! number of ice roughness categories

    real(wp), intent(in) :: cldfrac(:,:)       ! cloud fraction
                                               !    Dimensions: (ncol,nlayers)
    real(wp), intent(in) :: ciwp(:,:)          ! cloud ice water path
                                               !    Dimensions: (ncol,nlayers)
    real(wp), intent(in) :: clwp(:,:)          ! cloud liquid water path
                                               !    Dimensions: (ncol,nlayers)
    real(wp), intent(in) :: rei(:,:)           ! cloud ice particle effective size (microns)
                                               !    Dimensions: (ncol,nlayers)
    real(wp), intent(in) :: rel(:,:)           ! cloud liquid particle effective radius (microns)
                                               !    Dimensions: (ncol,nlayers)

  ! ------- Output -------
    class(ty_optical_props_arry), intent(inout) :: optical_props
                                               ! Dimensions: (ncol,nlayers,nbnd)

    character(len=128)    :: error_msg
    ! ------- Local -------
    real(wp) :: radliq                         ! cloud liquid droplet radius (microns)
    real(wp) :: radice                         ! cloud ice effective size (microns)
    real(wp) :: factor, fint

    logical :: cldmsk(ncol,nlayers)            ! cloud mask
    logical :: liqmsk(ncol,nlayers)            ! liquid cloud mask
    logical :: icemsk(ncol,nlayers)            ! ice cloud mask

    integer :: index                           !
    integer :: icol, ilyr, ibnd                !
    integer :: irad, irade, irads, iradg       !
    integer :: icergh                          ! ice surface roughness
                                               ! (1 = none, 2 = medium, 3 = high)

    real(wp) :: extliq(ncol,nlayers,nbnd)      ! liquid extinction coefficient
    real(wp) :: ssaliq(ncol,nlayers,nbnd)      ! liquid single scattering albedo
    real(wp) :: asyliq(ncol,nlayers,nbnd)      ! liquid asymmetry parameter

    real(wp) :: extice(ncol,nlayers,nbnd)      ! ice extinction coefficients
    real(wp) :: ssaice(ncol,nlayers,nbnd)      ! ice single scattering albedo
    real(wp) :: asyice(ncol,nlayers,nbnd)      ! ice asymmetry parameter

    real(wp) :: tauliq                         ! liquid cloud extinction optical depth
    real(wp) :: tauice                         ! ice cloud extinction optical depth

    real(wp) :: scatice                        ! Ice scattering term
    real(wp) :: scatliq                        ! Liquid scattering term

    real(wp) :: g                              ! asymmetry parameter - local

! ------- Definitions -------

! ------- Cloud masks -------
   cldmsk = (cldfrac >= 1.e-20_wp      .and. clwp + ciwp >= 1.e-20_wp)
   liqmsk = (cldmsk .and. rel > 0.0_wp .and. clwp >= 1.e-20_wp)
   icemsk = (cldmsk .and. rei > 0.0_wp .and. ciwp >= 1.e-20_wp)

! ------- Error checking -------
    error_msg = ''
    icergh = cloud_spec%icergh
    if (icergh < 1 .or. icergh > nrghice) then
       error_msg = 'cloud optics: cloud ice surface roughness flag is out of bounds'
       return
    endif
! For liquid OP, particle size is limited to radliq_lwr to radliq_upr microns
   if( any(liqmsk .and. rel < cloud_spec%radliq_lwr) .or. &
       any(liqmsk .and. rel > cloud_spec%radliq_upr) ) then
      error_msg = 'cloud optics: liquid effective radius is out of bounds'
      return
   endif
! For Yang (2013) ice OP, particle size is limited to radice_lwr to radice_upr microns
   if( any(icemsk .and. rei < cloud_spec%radice_lwr) .or. &
       any(icemsk .and. rei > cloud_spec%radice_upr) ) then
      error_msg = 'cloud optics: ice effective radius is out of bounds'
      return
   endif

    ! Initialize
    extliq(:,:,:) = 0.0_wp
    ssaliq(:,:,:) = 0.0_wp
    asyliq(:,:,:) = 0.0_wp
    extice(:,:,:) = 0.0_wp
    ssaice(:,:,:) = 0.0_wp
    asyice(:,:,:) = 0.0_wp
    !
    ! Cloud optical properties from LUT method
    !
    if (cloud_spec%do_lut) then
      do icol = 1, ncol
        do ilyr = 1, nlayers
          !
          ! Liquid optical properties
          !
          if (liqmsk(icol,ilyr)) then
            radliq = rel(icol,ilyr)
            factor = radliq - cloud_spec%radliq_fac
            index = int(factor)
            if (index .eq. 0) index = 1
            fint = factor - real(index)

            do ibnd = 1, nbnd
              extliq(icol,ilyr,ibnd) = cloud_spec%lut_extliq(index,ibnd) + &
                               fint * (cloud_spec%lut_extliq(index+1,ibnd) - &
                                       cloud_spec%lut_extliq(index,ibnd))
              ssaliq(icol,ilyr,ibnd) = cloud_spec%lut_ssaliq(index,ibnd) + &
                               fint * (cloud_spec%lut_ssaliq(index+1,ibnd) - &
                                       cloud_spec%lut_ssaliq(index,ibnd))
              asyliq(icol,ilyr,ibnd) = cloud_spec%lut_asyliq(index,ibnd) + &
                               fint * (cloud_spec%lut_asyliq(index+1,ibnd) - &
                                       cloud_spec%lut_asyliq(index,ibnd))
            enddo
          endif
          !
          ! Ice optical propertiesP for requested ice roughness (icergh)
          !
          if (icemsk(icol,ilyr)) then
            radice = rei(icol,ilyr)
            factor = radice * cloud_spec%radice_fac
            index = int(factor)
            fint = factor - real(index)

            do ibnd = 1, nbnd
              extice(icol,ilyr,ibnd) = cloud_spec%lut_extice(index,ibnd,icergh) + &
                               fint * (cloud_spec%lut_extice(index+1,ibnd,icergh) - &
                                       cloud_spec%lut_extice(index,ibnd,icergh))
              ssaice(icol,ilyr,ibnd) = cloud_spec%lut_ssaice(index,ibnd,icergh) + &
                               fint * (cloud_spec%lut_ssaice(index+1,ibnd,icergh) - &
                                       cloud_spec%lut_ssaice(index,ibnd,icergh))
              asyice(icol,ilyr,ibnd) = cloud_spec%lut_asyice(index,ibnd,icergh) + &
                               fint * (cloud_spec%lut_asyice(index+1,ibnd,icergh) - &
                                       cloud_spec%lut_asyice(index,ibnd,icergh))
            enddo
          endif
        enddo
      enddo
    else
      !
      ! Cloud optical properties from Pade coefficient method
      !
      do icol = 1, ncol
        do ilyr = 1, nlayers
          !
          ! Liquid optical properties
          !
          if (liqmsk(icol,ilyr)) then
            radliq = rel(icol,ilyr)
            !
            ! Define coefficient particle size regime for current size: extinction, ssa
            !
            error_msg = get_irad(cloud_spec, radliq, .True., 'ext', irade)
            if(error_msg /= "") return
            error_msg = get_irad(cloud_spec, radliq, .True., 'ssa', irads)
            if(error_msg /= "") return
            error_msg = get_irad(cloud_spec, radliq, .True., 'asy', iradg)
            if(error_msg /= "") return

            do ibnd = 1, nbnd
              extliq(icol,ilyr,ibnd) = cloud_spec%pade_ext(radliq, .True., ibnd, irade)
              ssaliq(icol,ilyr,ibnd) = cloud_spec%pade_ssa(radliq, .True., ibnd, irads)
              asyliq(icol,ilyr,ibnd) = cloud_spec%pade_asy(radliq, .True., ibnd, iradg)
            enddo
          endif
          !
          ! Ice optical properties
          !
          if (icemsk(icol,ilyr)) then
            radice = rei(icol,ilyr)
            ! Define coefficient particle size regime for current size: extinction, ssa
            error_msg = get_irad(cloud_spec, radice, .False., 'ext', irade)
            if(error_msg /= "") return
            error_msg = get_irad(cloud_spec, radice, .False., 'ssa', irads)
            if(error_msg /= "") return
            error_msg = get_irad(cloud_spec, radice, .False., 'asy', iradg)
            if(error_msg /= "") return

            do ibnd = 1, nbnd
             extice(icol,ilyr,ibnd) = cloud_spec%pade_ext(radice, .False., ibnd, irade, icergh)
             ssaice(icol,ilyr,ibnd) = cloud_spec%pade_ssa(radice, .False., ibnd, irads, icergh)
             asyice(icol,ilyr,ibnd) = cloud_spec%pade_asy(radice, .False., ibnd, iradg, icergh)
            enddo
          endif
        enddo
      enddo
    endif
    !
    ! Combine liquid and ice contributions into total cloud optical properties
    !
    select type(optical_props)
      type is (ty_optical_props_1scl)
        optical_props%tau(:,:,:) = 0.0_wp
      type is (ty_optical_props_2str)
        optical_props%tau(:,:,:) = 0.0_wp
        optical_props%ssa(:,:,:) = 0.0_wp
        optical_props%g  (:,:,:) = 0.0_wp
      type is (ty_optical_props_nstr)
        optical_props%tau(:,:,:) = 0.0_wp
        optical_props%ssa(:,:,:) = 0.0_wp
        optical_props%p  (:,:,:,:) = 0.0_wp
    end select
     do icol = 1, ncol
        do ilyr = 1, nlayers
           if (cldmsk(icol,ilyr)) then
              do ibnd = 1, nbnd
                 tauice = ciwp(icol,ilyr) * extice(icol,ilyr,ibnd)
                 tauliq = clwp(icol,ilyr) * extliq(icol,ilyr,ibnd)

                 select type(optical_props)
                 type is (ty_optical_props_1scl)
                   optical_props%tau(icol,ilyr,ibnd) = tauice + tauliq
                 type is (ty_optical_props_2str)
                   optical_props%tau(icol,ilyr,ibnd) = tauice + tauliq
                   scatice = ssaice(icol,ilyr,ibnd) * tauice
                   scatliq = ssaliq(icol,ilyr,ibnd) * tauliq
                   optical_props%ssa(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                         optical_props%tau(icol,ilyr,ibnd)
                   optical_props%g  (icol,ilyr,ibnd) = &
                        (scatice * asyice(icol,ilyr,ibnd) + scatliq * asyliq(icol,ilyr,ibnd)) / &
                        (scatice + scatliq)
                 type is (ty_optical_props_nstr)
                   optical_props%tau(icol,ilyr,ibnd) = tauice + tauliq
                   scatice = ssaice(icol,ilyr,ibnd) * tauice
                   scatliq = ssaliq(icol,ilyr,ibnd) * tauliq
                   optical_props%ssa(icol,ilyr,ibnd) = (scatice + scatliq) / &
                                         optical_props%tau(icol,ilyr,ibnd)
                   g   = &
                        (scatice * asyice(icol,ilyr,ibnd) + scatliq * asyliq(icol,ilyr,ibnd)) / &
                        (scatice + scatliq)
                   optical_props%p  (1,icol,ilyr,ibnd) = 1.0_wp
                   optical_props%p  (2,icol,ilyr,ibnd) = g
                   optical_props%p  (3,icol,ilyr,ibnd) = 0.1_wp
                 end select
              enddo
           endif
        enddo
     enddo

  end function cloud_optics

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function get_irad(cloud_spec,rad,is_liquid,param,irad_out) result(error_msg)

    class(ty_cloud_optics), intent(inout)      :: cloud_spec

    real(wp), intent(in)                       :: rad        ! particle radius
    logical, intent(in)                        :: is_liquid  ! T = liquid; F = ice
    character(len=3), intent(in)               :: param      ! ext/ssa/asy
    integer, intent(out)                       :: irad_out
    character(len=128)                         :: error_msg

    ! Local variables
    integer                                    :: irad, nrad
    integer, dimension(4)                      :: sizreg

    irad_out = 0
    error_msg = ""

    if (is_liquid) then
       nrad = 2
       if (param .eq. 'ext') sizreg = cloud_spec%pade_sizreg_extliq(:)
       if (param .eq. 'ssa') sizreg = cloud_spec%pade_sizreg_ssaliq(:)
       if (param .eq. 'asy') sizreg = cloud_spec%pade_sizreg_asyliq(:)
    else
       nrad = 3
       if (param .eq. 'ext') sizreg = cloud_spec%pade_sizreg_extice(:)
       if (param .eq. 'ssa') sizreg = cloud_spec%pade_sizreg_ssaice(:)
       if (param .eq. 'asy') sizreg = cloud_spec%pade_sizreg_asyice(:)
    endif

    do irad = 1, nrad
       if (rad .ge. sizreg(irad) .and. rad .le. sizreg(irad+1)) then
          irad_out = irad
          return
       endif
    enddo

! No size regime found
    error_msg = "cloud_optics, get_irad: out-of-range particle size value"

  end function get_irad

  !---------------------------------------------------------------------------
  function pade_ext(cloud_spec,reff,is_liquid,ibnd,irad,icergh)

    class(ty_cloud_optics), intent(inout) :: cloud_spec

    real(wp), intent(in) :: reff            ! particle radius (microns)
    logical, intent(in)  :: is_liquid       ! T = liquid; F = ice
    integer, intent(in) :: ibnd             ! band number
    integer, intent(in) :: irad             ! particle size regime
    integer, intent(in), optional :: icergh ! ice roughness index

    real(wp) :: pade_ext
    real(wp) :: reff2, reff3

    real(wp) :: p(6)                        ! extinction Pade coefficients

    if (is_liquid) then
       p(:) = cloud_spec%pade_extliq(ibnd,irad,:)
    else
       p(:) = cloud_spec%pade_extice(ibnd,irad,:,icergh)
    endif
    reff2 = reff * reff
    reff3 = reff2 * reff
! Pade formulation: Extinction Coefficient (Hogan and Bozzo, ECMWF, TM787, 2016)
    pade_ext = (p(1) + p(2)*reff + p(3)*reff2) / &
               (1.0_wp + p(4)*reff + p(5)*reff2 + p(6)*reff3)

  end function pade_ext

  !---------------------------------------------------------------------------
  function pade_ssa(cloud_spec,reff,is_liquid,ibnd,irad,icergh)

    class(ty_cloud_optics), intent(inout) :: cloud_spec

    real(wp), intent(in) :: reff            ! particle radius (microns)
    logical, intent(in)  :: is_liquid       ! T = liquid; F = ice
    integer, intent(in) :: ibnd             ! band number
    integer, intent(in) :: irad             ! particle size regime
    integer, intent(in), optional :: icergh ! ice roughness index

    real(wp) :: pade_ssa
    real(wp) :: reff2

    real(wp) :: p(5)                        ! ssa Pade coefficients

    if (is_liquid) then
       p(:) = cloud_spec%pade_ssaliq(ibnd,irad,:)
    else
       p(:) = cloud_spec%pade_ssaice(ibnd,irad,:,icergh)
    endif
    reff2 = reff * reff
! Pade formulation: Single Scattering Albedo (Hogan and Bozzo, ECMWF, TM787, 2016)
    pade_ssa = 1.0_wp - (p(1)+p(2)*reff+p(3)*reff2) / &
                        (1.0_wp+p(4)*reff+p(5)*reff2)
! Some values of ssa going slightly over 1.0 for Pade method. Replace small departures with 1.0 for now.
    if (pade_ssa >= 1.0_wp .and. pade_ssa <= 1.0005_wp) pade_ssa = 1.0_wp

  end function pade_ssa

  !---------------------------------------------------------------------------
  function pade_asy(cloud_spec,reff,is_liquid,ibnd,irad,icergh)

    class(ty_cloud_optics), intent(inout) :: cloud_spec

    real(wp), intent(in) :: reff            ! particle radius (microns)
    logical, intent(in)  :: is_liquid       ! T = liquid; F = ice
    integer, intent(in) :: ibnd             ! band number
    integer, intent(in) :: irad             ! particle size regime
    integer, intent(in), optional :: icergh ! ice roughness index

    real(wp) :: pade_asy
    real(wp) :: reff2

    real(wp) :: p(5)                        ! g Pade coefficients

    if (is_liquid) then
       p(:) = cloud_spec%pade_asyliq(ibnd,irad,:)
    else
       p(:) = cloud_spec%pade_asyice(ibnd,irad,:,icergh)
    endif
    reff2 = reff * reff
! Pade formulation: Asymmetry Parameter (Hogan and Bozzo, ECMWF, TM787, 2016)
    pade_asy = ((p(1)+p(2)*reff+p(3)*reff2) / &
               (1.0_wp+p(4)*reff+p(5)*reff2))

  end function pade_asy

end module mo_cloud_optics
