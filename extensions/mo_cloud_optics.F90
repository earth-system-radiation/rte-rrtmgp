! Module: mo_cloud_optics

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:
! This is the interface for routines that receive cloud physical properties
! and return cloud optical properties by band using LUT input data.
!

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
!    private

    ! User input
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    integer :: icergh                                   ! (1 = none, 2 = medium, 3 = high)

    ! Method for interpolation of cloud optical property coefficients to particle size
    logical :: do_lut                                   ! (.True. = LUT, .False. = Pade)

    ! Cloud physical properties                         ! (ncol,nlay)
    real(wp), dimension(:,:), allocatable :: cldfrac    ! cloud fraction
    real(wp), dimension(:,:), allocatable :: ciwp       ! cloud ice water path
    real(wp), dimension(:,:), allocatable :: clwp       ! cloud liquid water path
    real(wp), dimension(:,:), allocatable :: rei        ! cloud ice particle effective size (microns)
    real(wp), dimension(:,:), allocatable :: rel        ! cloud liquid particle effective radius (microns)

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
    generic,   public :: init_cldopt  => init_lut, init_pade
    procedure, public :: cloud_optics
    ! Internal procedures
    procedure, private :: init_lut
    procedure, private :: init_pade
    procedure, private :: pade_ext
    procedure, private :: pade_ssa
    procedure, private :: pade_asy
  end type ty_cloud_optics

contains
  ! ------------------------------------------------------------------------------
  !
  ! Cloud optics initialization function - LUT
  !
  ! ------------------------------------------------------------------------------
  function init_lut(cloud_spec, &
                    radliq_lwr, radliq_upr, radliq_fac, &
                    radice_lwr, radice_upr, radice_fac, &
                    lut_extliq, lut_ssaliq, lut_asyliq, &
                    lut_extice, lut_ssaice, lut_asyice) &
                    result(error_msg)
  ! ------------------------------------------------------------------------------
  ! Purpose:  Load lookup table cloud optics coefficients into cloud optics object

  ! ------- Input -------

    class(ty_cloud_optics),     intent(inout) :: cloud_spec   ! cloud specification data

    ! Lookup table interpolation constants
    real(wp),                   intent(in   ) :: radliq_lwr   ! liquid particle size lower bound for interpolation
    real(wp),                   intent(in   ) :: radliq_upr   ! liquid particle size upper bound for interpolation
    real(wp),                   intent(in   ) :: radliq_fac   ! constant for calculating interpolation indices for liquid
    real(wp),                   intent(in   ) :: radice_lwr   ! ice particle size lower bound for interpolation
    real(wp),                   intent(in   ) :: radice_upr   ! ice particle size upper bound for interpolation
    real(wp),                   intent(in   ) :: radice_fac   ! constant for calculating interpolation indices for ice
    ! LUT coefficients
    real(wp), dimension(:,:),   intent(in)    :: lut_extliq   ! extinction: liquid
    real(wp), dimension(:,:),   intent(in)    :: lut_ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:),   intent(in)    :: lut_asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:), intent(in)    :: lut_extice   ! extinction: ice
    real(wp), dimension(:,:,:), intent(in)    :: lut_ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:), intent(in)    :: lut_asyice   ! asymmetry parameter: ice

! ------- Local -------

    integer               :: nband, nrghice, nsize_liq, nsize_ice
    character(len=128)    :: error_msg

! ------- Definitions -------

! ------- Error checking -------
    error_msg = ''

    ! LUT coefficient dimensions
    nband     = size(lut_extliq,dim=2)
    nsize_liq = size(lut_extliq,dim=1)
    nsize_ice = size(lut_extice,dim=1)
    nrghice   = size(lut_extice,dim=3)

    ! Load LUT constants
    cloud_spec%radliq_lwr = radliq_lwr
    cloud_spec%radliq_upr = radliq_upr
    cloud_spec%radliq_fac = radliq_fac
    cloud_spec%radice_lwr = radice_lwr
    cloud_spec%radice_upr = radice_upr
    cloud_spec%radice_fac = radice_fac

    ! Allocate LUT coefficients
    allocate(cloud_spec%lut_extliq(nsize_liq, nband))
    allocate(cloud_spec%lut_ssaliq(nsize_liq, nband))
    allocate(cloud_spec%lut_asyliq(nsize_liq, nband))
    allocate(cloud_spec%lut_extice(nsize_ice, nband, nrghice))
    allocate(cloud_spec%lut_ssaice(nsize_ice, nband, nrghice))
    allocate(cloud_spec%lut_asyice(nsize_ice, nband, nrghice))

    ! Load LUT coefficients
    cloud_spec%lut_extliq = lut_extliq
    cloud_spec%lut_ssaliq = lut_ssaliq
    cloud_spec%lut_asyliq = lut_asyliq
    cloud_spec%lut_extice = lut_extice
    cloud_spec%lut_ssaice = lut_ssaice
    cloud_spec%lut_asyice = lut_asyice

  end function init_lut

  ! ------------------------------------------------------------------------------
  !
  ! Cloud optics initialization function - Pade
  !
  ! ------------------------------------------------------------------------------
  function init_pade(cloud_spec, &
                     radliq_lwr, radliq_upr, radice_lwr, radice_upr, &
                     pade_extliq, pade_ssaliq, pade_asyliq, &
                     pade_extice, pade_ssaice, pade_asyice, &
                     pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
                     pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice) &
                     result(error_msg)
  ! ------------------------------------------------------------------------------
  ! Purpose:  Load Pade cloud optics coefficients into cloud optics object

  ! ------- Input -------

    class(ty_cloud_optics),       intent(inout) :: cloud_spec    ! cloud specification data

    ! Particle size boundary limits
    real(wp),                   intent(in   ) :: radliq_lwr   ! liquid particle size lower bound for interpolation
    real(wp),                   intent(in   ) :: radliq_upr   ! liquid particle size upper bound for interpolation
    real(wp),                   intent(in   ) :: radice_lwr   ! ice particle size lower bound for interpolation
    real(wp),                   intent(in   ) :: radice_upr   ! ice particle size upper bound for interpolation

    real(wp), dimension(:,:,:),   intent(in)    :: pade_extliq   ! extinction: liquid
    real(wp), dimension(:,:,:),   intent(in)    :: pade_ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:,:),   intent(in)    :: pade_asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:,:), intent(in)    :: pade_extice   ! extinction: ice
    real(wp), dimension(:,:,:,:), intent(in)    :: pade_ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:,:), intent(in)    :: pade_asyice   ! asymmetry parameter: ice

    integer,  dimension(:),       intent(in)    :: pade_sizreg_extliq ! particle size boundaries, extinction: liquid
    integer,  dimension(:),       intent(in)    :: pade_sizreg_ssaliq ! particle size boundaries, single scattering albedo: liquid
    integer,  dimension(:),       intent(in)    :: pade_sizreg_asyliq ! particle size boundaries, asymmetry parameter: liquid
    integer,  dimension(:),       intent(in)    :: pade_sizreg_extice ! particle size boundaries, extinction: ice
    integer,  dimension(:),       intent(in)    :: pade_sizreg_ssaice ! particle size boundaries, single scattering albedo: ice
    integer,  dimension(:),       intent(in)    :: pade_sizreg_asyice ! particle size boundaries, asymmetry parameter: ice

! ------- Local -------

    integer               :: nband, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound
    character(len=128)    :: error_msg

! ------- Definitions -------

! ------- Error checking -------
    error_msg = ''

    ! Pade coefficient dimensions
    nband        = size(pade_extliq,dim=1)
    nsizereg     = size(pade_extliq,dim=2)
    ncoeff_ext   = size(pade_extliq,dim=3)
    ncoeff_ssa_g = size(pade_ssaliq,dim=3)
    nrghice      = size(pade_extice,dim=4)
    nbound       = size(pade_sizreg_extliq,dim=1)

    ! Load particle size boundaries
    cloud_spec%radliq_lwr = radliq_lwr
    cloud_spec%radliq_upr = radliq_upr
    cloud_spec%radice_lwr = radice_lwr
    cloud_spec%radice_upr = radice_upr

    ! Allocate Pade coefficients
    allocate(cloud_spec%pade_extliq(nband, nsizereg, ncoeff_ext))
    allocate(cloud_spec%pade_ssaliq(nband, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_asyliq(nband, nsizereg, ncoeff_ssa_g))
    allocate(cloud_spec%pade_extice(nband, nsizereg, ncoeff_ext, nrghice))
    allocate(cloud_spec%pade_ssaice(nband, nsizereg, ncoeff_ssa_g, nrghice))
    allocate(cloud_spec%pade_asyice(nband, nsizereg, ncoeff_ssa_g, nrghice))

    ! Load Pade coefficients
    cloud_spec%pade_extliq = pade_extliq
    cloud_spec%pade_ssaliq = pade_ssaliq
    cloud_spec%pade_asyliq = pade_asyliq
    cloud_spec%pade_extice = pade_extice
    cloud_spec%pade_ssaice = pade_ssaice
    cloud_spec%pade_asyice = pade_asyice

    ! Allocate Pade coefficient particle size regime boundaries
    allocate(cloud_spec%pade_sizreg_extliq(nbound))
    allocate(cloud_spec%pade_sizreg_ssaliq(nbound))
    allocate(cloud_spec%pade_sizreg_asyliq(nbound))
    allocate(cloud_spec%pade_sizreg_extice(nbound))
    allocate(cloud_spec%pade_sizreg_ssaice(nbound))
    allocate(cloud_spec%pade_sizreg_asyice(nbound))

    ! Load Pade coefficient particle size regime boundaries
    cloud_spec%pade_sizreg_extliq = pade_sizreg_extliq
    cloud_spec%pade_sizreg_ssaliq = pade_sizreg_ssaliq
    cloud_spec%pade_sizreg_asyliq = pade_sizreg_asyliq
    cloud_spec%pade_sizreg_extice = pade_sizreg_extice
    cloud_spec%pade_sizreg_ssaice = pade_sizreg_ssaice
    cloud_spec%pade_sizreg_asyice = pade_sizreg_asyice
  end function init_pade

  ! ------------------------------------------------------------------------------
  !
  ! Derive cloud optical properties from provided cloud physical properties
  !
  ! ------------------------------------------------------------------------------
  function cloud_optics( &
  ! Input
                   cloud_spec, &
                   ncol, nlayers, nbnd, nrghice, &
                   cldfrac, clwp, ciwp, rel, rei, &
  ! Output
                   optical_props) &
                   result(error_msg)
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


    real(wp) :: cwp(ncol,nlayers)              ! cloud water path (liquid + ice)
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

    real(wp) :: g                           ! asymmetry parameter - local

! ------- Definitions -------

! ------- Cloud masks -------
   cwp(:,:) = ciwp(:,:) + clwp(:,:)
   cldmsk = .False.
   liqmsk = .False.
   icemsk = .False.
   where (cldfrac >= 1.e-20_wp .and. cwp >= 1.e-20_wp) cldmsk = .True.
   where (cldmsk .and. rel > 0.0_wp .and. clwp >= 1.e-20_wp) liqmsk = .True.
   where (cldmsk .and. rei > 0.0_wp .and. ciwp >= 1.e-20_wp) icemsk = .True.

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

!
! Cloud optical properties from LUT method
!
    if (cloud_spec%do_lut) then

       ! Main column loop
       do icol = 1, ncol

          ! Main layer loop
          do ilyr = 1, nlayers
             if (cldmsk(icol,ilyr)) then
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

             endif

          ! End layer loop
          enddo

       ! End column loop
       enddo
!
! Check for non-positive cloud optical depth
!
       if (any(spread(cldmsk(:,:), dim=3, ncopies=nbnd) .and. &
               extliq * spread(clwp(:,:), dim=3, ncopies=nbnd) <= 0.0_wp .and. &
               extice * spread(ciwp(:,:), dim=3, ncopies=nbnd) <= 0.0_wp)) then
          error_msg = 'cloud optics: cloud optical depth is not positive'
          return
       endif
!
! Combine liquid and ice contributions into total cloud optical properties
!
       ! Main column loop
       do icol = 1, ncol

          ! Main layer loop
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

          ! End layer loop
          enddo

       ! End column loop
       enddo
!
! Cloud optical properties from Pade coefficient method
!
    else

       ! Main column loop
       do icol = 1, ncol

          ! Main layer loop
          do ilyr = 1, nlayers
             if (cldmsk(icol,ilyr)) then
!
! Liquid optical properties
!
                if (liqmsk(icol,ilyr)) then
                   radliq = rel(icol,ilyr)
! Define coefficient particle size regime for current size: extinction, ssa
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
! Derive optical properties for selected ice roughness
                      select case (icergh)

                      ! No ice roughness
                      case(1)
                         extice(icol,ilyr,ibnd) = cloud_spec%pade_ext(radice, .False., ibnd, irade, icergh)
                         ssaice(icol,ilyr,ibnd) = cloud_spec%pade_ssa(radice, .False., ibnd, irads, icergh)
                         asyice(icol,ilyr,ibnd) = cloud_spec%pade_asy(radice, .False., ibnd, iradg, icergh)

                      ! Medium ice roughness
                      case(2)
                         extice(icol,ilyr,ibnd) = cloud_spec%pade_ext(radice, .False., ibnd, irade, icergh)
                         ssaice(icol,ilyr,ibnd) = cloud_spec%pade_ssa(radice, .False., ibnd, irads, icergh)
                         asyice(icol,ilyr,ibnd) = cloud_spec%pade_asy(radice, .False., ibnd, iradg, icergh)

                      ! High ice roughness
                      case(3)
                         extice(icol,ilyr,ibnd) = cloud_spec%pade_ext(radice, .False., ibnd, irade, icergh)
                         ssaice(icol,ilyr,ibnd) = cloud_spec%pade_ssa(radice, .False., ibnd, irads, icergh)
                         asyice(icol,ilyr,ibnd) = cloud_spec%pade_asy(radice, .False., ibnd, iradg, icergh)

                      end select
                   enddo
                endif

             endif

          ! End layer loop
          enddo

       ! End column loop
       enddo
!
! Check for non-positive cloud optical depth
!
       if (any(spread(cldmsk(:,:), dim=3, ncopies=nbnd) .and. &
               extliq * spread(clwp(:,:), dim=3, ncopies=nbnd) <= 0.0_wp .and. &
               extice * spread(ciwp(:,:), dim=3, ncopies=nbnd) <= 0.0_wp)) then
          error_msg = 'cloud optics: cloud optical depth is not positive'
          return
       endif
!
! Combine liquid and ice contributions into total cloud optical properties
!
       ! Main column loop
       do icol = 1, ncol

          ! Main layer loop
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

          ! End layer loop
          enddo

       ! End column loop
       enddo

    endif

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
