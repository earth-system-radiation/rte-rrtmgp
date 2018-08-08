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
    private
    !
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    integer, public :: icergh = 0                                ! (1 = none, 2 = medium, 3 = high)

    ! Method for interpolation of cloud optical property coefficients to particle size
    logical, public :: do_lut                                   ! (.True. = LUT, .False. = Pade)

    ! Particle size boundary limits
    real(wp) :: radliq_lwr = 0._wp              ! liquid particle size lower bound for interpolation
    real(wp) :: radliq_upr = 0._wp              ! liquid particle size upper bound for interpolation
    real(wp) :: radice_lwr = 0._wp              ! ice particle size lower bound for interpolation
    real(wp) :: radice_upr = 0._wp              ! ice particle size upper bound for interpolation
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
    real(wp), dimension(:), allocatable :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq  ! (nbound)
    real(wp), dimension(:), allocatable :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice  ! (nbound)

! ------------------------------------------------------------------------------------------
  contains
    generic,   public :: load  => load_lut, load_pade
    procedure, public :: finalize
    procedure, public :: cloud_optics
    procedure, public :: get_min_radius_liq
    procedure, public :: get_min_radius_ice
    procedure, public :: get_max_radius_liq
    procedure, public :: get_max_radius_ice
    procedure, public :: get_num_ice_roughness_types
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
                     pade_extliq, pade_ssaliq, pade_asyliq, &
                     pade_extice, pade_ssaice, pade_asyice, &
                     pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
                     pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice) &
                     result(error_msg)
    class(ty_cloud_optics),       intent(inout) :: this          ! cloud specification data
    real(wp), dimension(:,:),     intent(in   ) :: band_lims_wvn ! Spectral discretization
    !
    ! Pade coefficients: extinction, single-scattering albedo, and asymmetry factor for liquid and ice
    !
    real(wp), dimension(:,:,:),   intent(in)    :: pade_extliq, pade_ssaliq, pade_asyliq
    real(wp), dimension(:,:,:,:), intent(in)    :: pade_extice, pade_ssaice, pade_asyice
    !
    ! Boundaries of size regimes. Liquid and ice are separate;
    !   extinction is fit to different numbers of size bins than single-scattering albedo and asymmetry factor
    !
    real(wp),  dimension(:),       intent(in)    :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq
    real(wp),  dimension(:),       intent(in)    :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice
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
    nbound       = size(pade_sizreg_extliq)
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
      if(nbound /= 4) &
        error_msg = "cloud_optics%init(): Expecting precisely three size regimes for Pade approximants"
    if(error_msg /= "") return
    !
    ! Consistency among size regimes
    !
    this%radliq_lwr = pade_sizreg_extliq(1)
    this%radliq_upr = pade_sizreg_extliq(nbound)
    this%radice_lwr = pade_sizreg_extice(1)
    this%radice_upr = pade_sizreg_extice(nbound)

    if(any([pade_sizreg_ssaliq(1), pade_sizreg_asyliq(1)] < this%radliq_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have inconsistent lowest values"
    if(any([pade_sizreg_ssaice(1), pade_sizreg_asyice(1)] < this%radice_lwr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have inconsistent lower values"

    if(any([pade_sizreg_ssaliq(nbound), pade_sizreg_asyliq(nbound)] > this%radliq_upr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_upr"
    if(any([pade_sizreg_ssaice(nbound), pade_sizreg_asyice(nbound)] > this%radice_upr)) &
      error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_upr"
    if(error_msg /= "") return

    ! Allocate Pade coefficients
    allocate(this%pade_extliq(nband, nsizereg, ncoeff_ext),   &
             this%pade_ssaliq(nband, nsizereg, ncoeff_ssa_g), &
             this%pade_asyliq(nband, nsizereg, ncoeff_ssa_g), &
             this%pade_extice(nband, nsizereg, ncoeff_ext,   nrghice), &
             this%pade_ssaice(nband, nsizereg, ncoeff_ssa_g, nrghice), &
             this%pade_asyice(nband, nsizereg, ncoeff_ssa_g, nrghice))

    ! Allocate Pade coefficient particle size regime boundaries
    allocate(this%pade_sizreg_extliq(nbound), &
             this%pade_sizreg_ssaliq(nbound), &
             this%pade_sizreg_asyliq(nbound), &
             this%pade_sizreg_extice(nbound), &
             this%pade_sizreg_ssaice(nbound), &
             this%pade_sizreg_asyice(nbound))

    ! Load Pade coefficients
    this%pade_extliq = pade_extliq
    this%pade_ssaliq = pade_ssaliq
    this%pade_asyliq = pade_asyliq
    this%pade_extice = pade_extice
    this%pade_ssaice = pade_ssaice
    this%pade_asyice = pade_asyice

    ! Load Pade coefficient particle size regime boundaries
    this%pade_sizreg_extliq = pade_sizreg_extliq
    this%pade_sizreg_ssaliq = pade_sizreg_ssaliq
    this%pade_sizreg_asyliq = pade_sizreg_asyliq
    this%pade_sizreg_extice = pade_sizreg_extice
    this%pade_sizreg_ssaice = pade_sizreg_ssaice
    this%pade_sizreg_asyice = pade_sizreg_asyice
  end function load_pade
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Finalize
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine finalize(this)
    class(ty_cloud_optics), intent(inout) :: this

    this%radliq_lwr = 0._wp
    this%radliq_upr = 0._wp
    this%radice_lwr = 0._wp
    this%radice_upr = 0._wp
    this%radliq_fac = 0._wp
    this%radice_fac = 0._wp

    ! Lookup table cloud optics coefficients
    if(allocated(this%lut_extliq)) &
      deallocate(this%lut_extliq, this%lut_ssaliq, this%lut_asyliq, &
                 this%lut_extice, this%lut_ssaice, this%lut_asyice)

  ! Pade cloud optics coefficients
    if(allocated(this%pade_extliq)) then
      deallocate(this%pade_extliq, this%pade_ssaliq, this%pade_asyliq, &
                 this%pade_extice, this%pade_ssaice, this%pade_asyice, &
                 this%pade_sizreg_extliq, this%pade_sizreg_ssaliq, this%pade_sizreg_asyliq, &
                 this%pade_sizreg_extice, this%pade_sizreg_ssaice, this%pade_sizreg_asyice)
    end if
  end subroutine finalize
  ! ------------------------------------------------------------------------------
  !
  ! Derive cloud optical properties from provided cloud physical properties
  !
  ! ------------------------------------------------------------------------------
  !
  ! Compute single-scattering properties
  !
  function cloud_optics(this, &
                        ncol, nlay, nbnd, nrghice, &
                        liqmsk, icemsk,   &
                        clwp, ciwp, rel, rei, optical_props) result(error_msg)
    class(ty_cloud_optics), &
              intent(inout) :: this
    integer,  intent(in   ) :: ncol, nlay, nbnd
    integer,  intent(in   ) :: nrghice              ! number of ice roughness categories
    logical,  intent(in   ) :: liqmsk(ncol,nlay), & ! Cloud mask for liquid and ice clouds respectively
                               icemsk(ncol,nlay)
    real(wp), intent(in   ) :: ciwp(ncol,nlay), &     ! cloud ice water path
                               clwp(ncol,nlay), &     ! cloud liquid water path
                               rei(ncol,nlay), &      ! cloud ice particle effective size (microns)
                               rel(ncol,nlay)      ! cloud liquid particle effective radius (microns)
    class(ty_optical_props_arry), &
              intent(inout) :: optical_props
                                               ! Dimensions: (ncol,nlay,nbnd)

    character(len=128)    :: error_msg
    ! ------- Local -------
    integer, parameter                  :: max_re_moments = 3
    real(wp), dimension(max_re_moments) :: re_moments ! re, re**2, re**3 etc.
    real(wp) :: radliq                         ! cloud liquid droplet radius (microns)
    real(wp) :: radice                         ! cloud ice effective size (microns)
    real(wp) :: factor, fint

    integer :: index                           !
    integer :: icol, ilyr, ibnd, i             !
    integer :: irad, irade, irads, iradg       !
    integer :: icergh                          ! ice surface roughness
                                               ! (1 = none, 2 = medium, 3 = high)

    real(wp) :: extliq(ncol,nlay,nbnd)      ! liquid extinction coefficient
    real(wp) :: ssaliq(ncol,nlay,nbnd)      ! liquid single scattering albedo
    real(wp) :: asyliq(ncol,nlay,nbnd)      ! liquid asymmetry parameter

    real(wp) :: extice(ncol,nlay,nbnd)      ! ice extinction coefficients
    real(wp) :: ssaice(ncol,nlay,nbnd)      ! ice single scattering albedo
    real(wp) :: asyice(ncol,nlay,nbnd)      ! ice asymmetry parameter

    real(wp) :: tauliq                         ! liquid cloud extinction optical depth
    real(wp) :: tauice                         ! ice cloud extinction optical depth

    real(wp) :: scatice                        ! Ice scattering term
    real(wp) :: scatliq                        ! Liquid scattering term

    real(wp) :: g                              ! asymmetry parameter - local
    integer  :: nsizereg
! ------- Definitions -------

! ------- Error checking -------
    error_msg = ''
    icergh = this%icergh
    if (icergh < 1 .or. icergh > nrghice) &
       error_msg = 'cloud optics: cloud ice surface roughness flag is out of bounds'

    if(any(liqmsk .and. (rel < this%radliq_lwr .or. rel > this%radliq_upr))) &
      error_msg = 'cloud optics: liquid effective radius is out of bounds'

    if(any(icemsk .and. (rei < this%radice_lwr .or. rei > this%radice_upr))) &
      error_msg = 'cloud optics: ice effective radius is out of bounds'

    if(any((liqmsk .and.  clwp < 0._wp) .or. (icemsk .and.  ciwp < 0._wp))) &
      error_msg = 'cloud optics: negative clwp or ciwp where clouds are supposed to be'
    if(error_msg /= "") return

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
    if (this%do_lut) then
      do icol = 1, ncol
        do ilyr = 1, nlay
          !
          ! Liquid optical properties
          !
          if (liqmsk(icol,ilyr)) then
            radliq = rel(icol,ilyr)
            factor = radliq - this%radliq_fac
            index = int(factor)
            if (index .eq. 0) index = 1
            fint = factor - real(index)

            do ibnd = 1, nbnd
              extliq(icol,ilyr,ibnd) = table_interp(this%lut_extliq(:,ibnd), index, fint)
              ssaliq(icol,ilyr,ibnd) = table_interp(this%lut_ssaliq(:,ibnd), index, fint)
              asyliq(icol,ilyr,ibnd) = table_interp(this%lut_asyliq(:,ibnd), index, fint)
           end do
          endif
          !
          ! Ice optical properties for requested ice roughness (icergh)
          !
          if (icemsk(icol,ilyr)) then
            radice = rei(icol,ilyr)
            factor = radice * this%radice_fac
            index = int(factor)
            fint = factor - real(index)

            do ibnd = 1, nbnd
              extice(icol,ilyr,ibnd) = table_interp(this%lut_extice(:,ibnd,icergh), index, fint)
              ssaice(icol,ilyr,ibnd) = table_interp(this%lut_ssaice(:,ibnd,icergh), index, fint)
              asyice(icol,ilyr,ibnd) = table_interp(this%lut_asyice(:,ibnd,icergh), index, fint)
           end do
          endif
        enddo
      enddo
    else
      !
      ! Cloud optical properties from Pade coefficient method
      !
      ! This assumes that all the Pade treaments have the same number of size regimes
      nsizereg = size(this%pade_sizreg_extliq)-1
      do icol = 1, ncol
        do ilyr = 1, nlay
          !
          ! Liquid optical properties
          !
          if (liqmsk(icol,ilyr)) then
            radliq = rel(icol,ilyr)
            re_moments(:max_re_moments) = [radliq, radliq*radliq, radliq*radliq*radliq]

            !
            ! Define coefficient particle size regime for current size: extinction, ssa
            !
            irade = get_irad(radliq, this%pade_sizreg_extliq)
            irads = get_irad(radliq, this%pade_sizreg_ssaliq)
            iradg = get_irad(radliq, this%pade_sizreg_asyliq)

            do ibnd = 1, nbnd
              extliq(icol,ilyr,ibnd) = this%pade_ext(radliq, .True., ibnd, irade)
              ssaliq(icol,ilyr,ibnd) = this%pade_ssa(radliq, .True., ibnd, irads)
              asyliq(icol,ilyr,ibnd) = this%pade_asy(radliq, .True., ibnd, iradg)
            enddo
!            extliq(icol,ilyr,1:nbnd) = pade_eval(nbnd,nsizereg,2,3,irade,re_moments,this%pade_extliq)
!            ssaliq(icol,ilyr,1:nbnd) = 1._wp - max(0._wp,
!                                       pade_eval(nbnd,nsizereg,2,3,irads,re_moments,this%pade_ssaliq)
!                                                  )
!            asyliq(icol,ilyr,1:nbnd) = pade_eval(nbnd,nsizereg,2,2,iradg,re_moments,this%pade_asyliq)
          endif
          !
          ! Ice optical properties
          !
          if (icemsk(icol,ilyr)) then
            radice = rei(icol,ilyr)
            re_moments(:max_re_moments) = [radice, radice*radice, radice*radice*radice]
            irade = get_irad(radice, this%pade_sizreg_extice)
            irads = get_irad(radice, this%pade_sizreg_ssaice)
            iradg = get_irad(radice, this%pade_sizreg_asyice)

            do ibnd = 1, nbnd
             extice(icol,ilyr,ibnd) = this%pade_ext(radice, .False., ibnd, irade, icergh)
             ssaice(icol,ilyr,ibnd) = this%pade_ssa(radice, .False., ibnd, irads, icergh)
             asyice(icol,ilyr,ibnd) = this%pade_asy(radice, .False., ibnd, iradg, icergh)
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
        do ilyr = 1, nlay
           if (liqmsk(icol,ilyr) .or. icemsk(icol,ilyr)) then
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
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function get_min_radius_liq(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radliq_lwr
  end function get_min_radius_liq
  !-----------------------------------------------
  function get_max_radius_liq(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radliq_upr
  end function get_max_radius_liq
  !-----------------------------------------------
  function get_min_radius_ice(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radice_lwr
  end function get_min_radius_ice
  !-----------------------------------------------
  function get_max_radius_ice(this) result(r)
    class(ty_cloud_optics), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radice_upr
  end function get_max_radius_ice
  !-----------------------------------------------
  function get_num_ice_roughness_types(this) result(i)
    class(ty_cloud_optics), intent(in   ) :: this
    integer                               :: i

    i = 0
    if(allocated(this%pade_extice)) i = size(this%pade_extice, dim=4)
    if(allocated(this%lut_extice )) i = size(this%lut_extice,  dim=3)
  end function get_num_ice_roughness_types
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Ancillary functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function get_irad(rad, sizereg) result(irad)
    real(wp),               intent(in) :: rad
    real(wp), dimension(:), intent(in) :: sizereg
    integer                            :: irad
    !
    ! Finds index into size regime table
    ! This works only if there are precisely three size regimes (four bounds) and it's
    !   previously guaranteed that sizereg(1) <= rad <= sizereg(4)
    !
    irad = min(floor((rad - sizereg(2))/sizereg(3)) + 2, 3)
  end function get_irad
  !---------------------------------------------------------------------------
  !
  ! Evaluate Pade approximant of order [m/n] = [2/2] or [2/3]
  !   It might be better to write this as a loop for general order
  !
  function pade_eval(nband, nrads, m, n, irad, re_moments, coeff_table)
    integer,                intent(in) :: nband, nrads, m, n, irad
    real(wp), dimension(nband, nrads, m+n+1), &
                            intent(in) :: coeff_table
    real(wp), dimension(           max(m,n)), &
                            intent(in) :: re_moments ! [re, re**2, re**3] etc.
    real(wp), dimension(nband)         :: pade_eval

    integer :: iband

    if(n+m == 5) then
      do iband = 1, nband
        pade_eval(iband) = (coeff_table(iband,irad,1)               +  &
                            coeff_table(iband,irad,2)*re_moments(1) +  &
                            coeff_table(iband,irad,3)*re_moments(2)) / &
                           (1.0_wp                                  +  &
                            coeff_table(iband,irad,4)*re_moments(1) +  &
                            coeff_table(iband,irad,5)*re_moments(2))
      end do
    else if(n+m == 6) then
      pade_eval(iband) = (coeff_table(iband,irad,1)               +  &
                          coeff_table(iband,irad,2)*re_moments(1) +  &
                          coeff_table(iband,irad,3)*re_moments(2)) / &
                         (1.0_wp                                  +  &
                          coeff_table(iband,irad,4)*re_moments(1) +  &
                          coeff_table(iband,irad,5)*re_moments(2) +  &
                          coeff_table(iband,irad,6)*re_moments(3))
    end if
!    if(any(pade_eval < 0._wp)) then
!      print *, re_moments(1), irad
!      print *, pack([(iband, iband = 1, nband)], pade_eval < 0._wp)
!    end if
  end function pade_eval
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
  !---------------------------------------------------------------------------
  function table_interp(table, index, fint) result(x)
    real(wp), dimension(:), intent(in) :: table
    integer,                  intent(in) :: index
    real(wp),                 intent(in) :: fint
    real(wp)                             :: x

    x = table(index) + fint * (table(index+1) - table(index))
  end function table_interp
end module mo_cloud_optics
