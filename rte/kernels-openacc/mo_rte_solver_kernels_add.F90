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
!
! Numeric calculations for radiative transfer solvers.
!   Emission/absorption (no-scattering) calculations
!     solver for multi-angle Gaussian quadrature
!     solver for a single angle, calling
!       source function computation (linear-in-tau)
!       transport
!   Extinction-only calculation (direct solar beam)
!   Two-stream calculations
!     solvers for LW and SW with different boundary conditions and source functions
!       source function calculation for LW, SW
!       two-stream calculations for LW, SW (using different assumtions about phase function)
!       transport (adding)
!   Application of boundary conditions
!
! -------------------------------------------------------------------------------------------------
module mo_rte_solver_kernels_add
  use,  intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp, wl
  use mo_rte_solver_kernels, only : apply_BC
  implicit none
  private

  public :: lw_solver_1rescl_GaussQuad,  lw_solver_1rescl

  real(wp), parameter :: pi = acos(-1._wp)
contains
   
! ++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  TANG + IP
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! -------------------------------------------------------------------------------------------------
  !
  ! Top-level longwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW fluxes, no scattering, mu (cosine of integration angle) specified by column
  !   Does radiation calculation at user-supplied angles; converts radiances to flux
  !   using user-supplied weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_1rescl(ncol, nlay, ngpt, top_at_1, D, weight,                             &
                              tau, scaling, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                              radn_up, radn_dn) bind(C, name="lw_solver_1rescl")
  use mo_rte_solver_kernels, only : lw_source_noscat, lw_transport_noscat

    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: D            ! secant of propagation angle  []
    real(wp),                              intent(in   ) :: weight       ! quadrature weight
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: scaling          ! single scattering albedo []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
    ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
    ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
    real(wp), dimension(ncol,nlay,  ngpt), target, &
                                           intent(in   ) :: lev_source_inc, lev_source_dec
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis     ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src      ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn      ! Top level must contain incident flux boundary condition

    ! Local variables, WITH g-point dependency
    real(wp), dimension(ncol,nlay,ngpt) :: tau_loc, &  ! path length (tau/mu)
                                             trans       ! transmissivity  = exp(-tau)
    real(wp), dimension(ncol,nlay,ngpt) :: source_dn, source_up
    real(wp), dimension(ncol,     ngpt) :: source_sfc, sfc_albedo

    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn ! Mapping increasing/decreasing indicies to up/down

    real(wp), parameter :: pi = acos(-1._wp)
    integer             :: ilev, igpt, top_level
    ! ------------------------------------
    real(wp)            :: fact
    real(wp), parameter :: tau_thresh = sqrt(epsilon(tau))
    integer             :: icol
    real(wp), dimension(ncol     ) :: sfcSource
    ! ------------------------------------

    ! Which way is up?
    ! Level Planck sources for upward and downward radiation
    ! When top_at_1, lev_source_up => lev_source_dec
    !                lev_source_dn => lev_source_inc, and vice-versa
    if(top_at_1) then
      top_level = 1
      lev_source_up => lev_source_dec
      lev_source_dn => lev_source_inc
    else
      top_level = nlay+1
      lev_source_up => lev_source_inc
      lev_source_dn => lev_source_dec
    end if

    !$acc enter data copyin(d,tau,sfc_src,sfc_emis,lev_source_dec,lev_source_inc,lay_source,radn_dn)
    !$acc enter data create(tau_loc,trans,source_dn,source_up,source_sfc,sfc_albedo,radn_up)
    !$acc enter data attach(lev_source_up,lev_source_dn)

    ! NOTE: This kernel produces small differences between GPU and CPU
    ! implementations on Ascent with PGI, we assume due to floating point
    ! differences in the exp() function. These differences are small in the
    ! RFMIP test case (10^-6).
    !$acc parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilev = 1, nlay
        do icol = 1, ncol
          !
          ! Optical path and transmission, used in source function and transport calculations
          !
          tau_loc(icol,ilev,igpt) = tau(icol,ilev,igpt)*D(icol,igpt)
          trans  (icol,ilev,igpt) = exp(-tau_loc(icol,ilev,igpt))
        end do
      end do
    end do

    !$acc parallel loop collapse(2)
    do igpt = 1, ngpt
      do icol = 1, ncol
      !
      ! Transport is for intensity
      !   convert flux at top of domain to intensity assuming azimuthal isotropy
      !
        radn_dn(icol,top_level,igpt) = radn_dn(icol,top_level,igpt)/(2._wp * pi * weight)
      !
        ! Surface albedo, surface source function
      !
        sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt)
        source_sfc(icol,igpt) = sfc_emis(icol,igpt) * sfc_src(icol,igpt)
      end do
    end do


    !
    ! Source function for diffuse radiation
    !
    call lw_source_noscat(ncol, nlay, ngpt, &
                          lay_source, lev_source_up, lev_source_dn, &
                          tau_loc, trans, source_dn, source_up)

    !
    ! Transport
    !
    call lw_transport_1rescl(ncol, nlay, ngpt, top_at_1,  &
                             tau_loc, scaling, trans, &
                             sfc_albedo, source_dn, source_up, source_sfc, &
                             radn_up, radn_dn)

      ! Convert intensity to flux assuming azimuthal isotropy and quadrature weight
      !
      ! radn_dn(:,:,igpt) = 2._wp * pi * weight * radn_dn(:,:,igpt)
      ! radn_up(:,:,igpt) = 2._wp * pi * weight * radn_up(:,:,igpt)
  end subroutine lw_solver_1rescl
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW transport, no scattering, multi-angle quadrature
  !   Users provide a set of weights and quadrature angles
  !   Routine sums over single-angle solutions for each sets of angles/weights
  !
  ! ---------------------------------------------------------------
 
  subroutine lw_solver_1rescl_GaussQuad(ncol, nlay, ngpt, top_at_1, nmus, Ds, weights, &
                                   tau, scaling, lay_source, lev_source_inc, lev_source_dec, &
                                   sfc_emis, sfc_src,&
                                  flux_up, flux_dn) &
                                   bind(C, name="lw_solver_1rescl_GaussQuad")
    use mo_rte_solver_kernels, only : lw_solver_noscat
    use mo_rte_solver_kernels, only : apply_BC
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    integer,                               intent(in   ) :: nmus         ! number of quadrature angles
    real(wp), dimension(nmus),             intent(in   ) :: Ds, weights  ! quadrature secants, weights
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: scaling          ! single scattering albedo []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis     ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src      ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn      ! Top level must contain incident flux boundary condition
    ! Local variables
    real(wp), dimension(ncol,nlay+1,ngpt) :: radn_dn, radn_up ! Fluxes per quad angle
    real(wp), dimension(ncol,       ngpt) :: Ds_ncol

    integer :: imu, top_level
    real    :: weight
    ! ------------------------------------
    !
    ! For the first angle output arrays store total flux
    !
    top_level = MERGE(1, nlay+1, top_at_1)
    Ds_ncol(:,:) = Ds(1)
    weight = 2._wp*pi*weights(1)
    radn_dn(1:ncol, top_level, 1:ngpt)  = flux_dn(1:ncol, top_level, 1:ngpt) / weight

    call lw_solver_1rescl(ncol, nlay, ngpt, &
                          top_at_1, Ds_ncol, weights(1), tau, scaling, &
                          lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                          flux_up, flux_dn)

    flux_up = flux_up * weight
    flux_dn = flux_dn * weight

    do imu = 2, nmus
      Ds_ncol(:,:) = Ds(imu)
      weight = 2._wp*pi*weights(imu)
      radn_dn(1:ncol, top_level, 1:ngpt)  = flux_dn(1:ncol, top_level, 1:ngpt) / weight
      call lw_solver_1rescl(ncol, nlay, ngpt, &
                            top_at_1, Ds_ncol, weights(imu), tau, scaling, &
                            lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                            radn_up, radn_dn)

      flux_up(:,:,:) = flux_up(:,:,:) + weight*radn_up(:,:,:)
      flux_dn(:,:,:) = flux_dn(:,:,:) + weight*radn_dn(:,:,:)
    end do
  end subroutine lw_solver_1rescl_GaussQuad

    ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave no-scattering transport
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_transport_1rescl(ncol, nlay, ngpt, top_at_1, &
                                 tau, scaling, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                                 radn_up, radn_dn) bind(C, name="lw_transport_1rescl")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1   !
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                       trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: scaling        ! single scattering
    real(wp), dimension(ncol       ,ngpt), intent(in   ) :: sfc_albedo ! Surface albedo
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: source_dn, &
                                                       source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ,ngpt), intent(in   ) :: source_sfc ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_up    ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn    !Top level must contain incident flux boundary condition
    ! Local variables
    integer :: ilev, icol, igpt
    ! ---------------------------------------------------
    real(wp) :: xx
    if(top_at_1) then
      !
      ! Top of domain is index 1
      !
      ! Downward propagation
      !$acc  parallel loop collapse(2)
      do igpt = 1, ngpt
        do icol = 1, ncol

          do ilev = 1, nlay
            radn_dn(icol,ilev+1,igpt) = trans(icol,ilev,igpt)  *radn_dn(icol,ilev,igpt)  + source_dn(icol,ilev,igpt)
          end do

          ! Surface reflection and emission
          radn_up(icol,nlay+1,igpt) = radn_dn(icol,nlay+1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt)

          ! 1st Upward propagation
          do ilev = nlay, 1, -1
            radn_up(icol,ilev,igpt) = trans(icol,ilev,igpt)*radn_up(icol,ilev+1,igpt) + source_up(icol,ilev,igpt)

            if ( scaling(icol,ilev,igpt) > 1e-6 )  then
          !  
          ! here scaling is used to store parameter wb/[(]1-w(1-b)] of Eq.21 of the Tang's paper
          ! explanation of factor 0.4 note A of Table
          !
              xx = 0.4_wp*scaling(icol,ilev,igpt)*&
                     ( radn_dn(icol,ilev,igpt)*(1.-trans(icol,ilev,igpt)**2 ) - &
                       source_dn(icol,ilev,igpt)  *trans(icol,ilev,igpt ) - &
                       source_up(icol,ilev,igpt))
            endif  
          enddo  
          ! 2nd Downward propagation
          do ilev = 1, nlay
            radn_dn(icol,ilev+1,igpt) = trans(icol,ilev,igpt)*radn_dn(icol,ilev,igpt) + source_dn(icol,ilev,igpt)
            if ( scaling(icol,ilev,igpt) > 1e-6 )  then
          !  
          ! here scaling is used to store parameter wb/[(]1-w(1-b)] of Eq.21 of the Tang's paper
          ! explanation of factor 0.4 note A of Table
          !
                xx = 0.4_wp*scaling(icol,ilev,igpt)*( &
                    radn_up(icol,ilev,igpt)*(1. -trans(icol,ilev,igpt)**2)  - &
                    source_up(icol,ilev,igpt)*trans(icol,ilev,igpt) - &
                    source_dn(icol,ilev,igpt) )
                  radn_dn(icol,ilev+1,igpt) = radn_dn(icol,ilev+1,igpt) + xx
            endif  
          enddo  
        enddo
      enddo
    else
      !
      ! Top of domain is index nlay+1
      !
      ! Downward propagation
      !
      ! --------- N+1
      !                   layer N
      ! ----------N
      !
      !
      !$acc  parallel loop collapse(2)
      do igpt = 1, ngpt
        do icol = 1, ncol
          do ilev = nlay, 1, -1
            radn_dn(icol,ilev,igpt) = trans(icol,ilev,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt)
          end do

          ! Surface reflection and emission
          radn_up(icol,1,igpt) = radn_dn(icol,1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt)
          ! Upward propagation
          do ilev = 1, nlay
            radn_up(icol,ilev+1,igpt) =  trans(icol,ilev,igpt) * radn_up(icol,ilev,igpt) +  source_up(icol,ilev,igpt)
            if ( scaling(icol,ilev,igpt) > 1e-6 )  then
          !  
          ! here scaling is used to store parameter wb/[(]1-w(1-b)] of Eq.21 of the Tang's paper
          ! explanation of factor 0.4 note A of Table
          !
               xx = 0.4_wp*scaling(icol,ilev,igpt)*&
                      ( radn_dn(icol,ilev+1,igpt)*(1.-trans(icol,ilev,igpt)**2 ) - &
                        source_dn(icol,ilev,igpt) *trans(icol,ilev ,igpt) - &
                        source_up(icol,ilev,igpt))
               radn_up(icol,ilev+1,igpt) = radn_up(icol,ilev+1,igpt) + xx
           endif  
          end do

          ! 2st Downward propagation
          do ilev = nlay, 1, -1
            radn_dn(icol,ilev,igpt) = trans(icol,ilev,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt)
             if ( scaling(icol,ilev,igpt) > 1e-6 )  then
          !  
          ! here scaling is used to store parameter wb/[(]1-w(1-b)] of Eq.21 of the Tang's paper
          ! explanation of factor 0.4 note A of Table
          !
                       xx = 0.4_wp*scaling(icol,ilev,igpt)*( &
                        radn_up(icol,ilev,igpt)*(1.-trans(icol,ilev,igpt)**2)  - &
                        source_up(icol,ilev,igpt)*trans(icol,ilev ,igpt ) - &
                        source_dn(icol,ilev,igpt) )
                radn_dn(icol,ilev,igpt) = radn_dn(icol,ilev,igpt) + xx
            endif  
          end do
        enddo
      enddo
    end if
  end subroutine lw_transport_1rescl
end module mo_rte_solver_kernels_add
