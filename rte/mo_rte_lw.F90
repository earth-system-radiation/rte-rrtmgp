! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2021,  Atmospheric and Environmental Research,
! Regents of the University of Colorado, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!  Contains a single routine to compute direct and diffuse fluxes of solar radiation given
!    atmospheric optical properties, spectrally-resolved
!    information about vertical ordering
!    internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
!    boundary conditions: surface emissivity defined per band
!    optionally, a boundary condition for incident diffuse radiation
!    optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected
!
! If optical properties are supplied via class ty_optical_props_1scl (absorption optical thickenss only)
!    then an emission/absorption solver is called
!    If optical properties are supplied via class ty_optical_props_2str fluxes are computed via
!    two-stream calculations and adding.
!
! It is the user's responsibility to ensure that emissivity is on the same
!   spectral grid as the optical properties.
!
! Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
!   whatever summary the user needs.
!
! The routine does error checking and choses which lower-level kernel to invoke based on
!   what kinds of optical properties are supplied
!
! -------------------------------------------------------------------------------------------------
module mo_rte_lw
  use mo_rte_kind,      only: wp, wl
  use mo_rte_config,    only: check_extents, check_values
  use mo_rte_util_array,only: any_vals_less_than, any_vals_outside, extents_are, zero_array
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_source_functions,   &
                        only: ty_source_func_lw
  use mo_fluxes,        only: ty_fluxes, ty_fluxes_broadband
  use mo_rte_solver_kernels, &
                        only: lw_solver_noscat, lw_solver_noscat_GaussQuad, lw_solver_2stream
  implicit none
  private

  public :: rte_lw
contains
  ! --------------------------------------------------
  !
  ! Interface using only optical properties and source functions as inputs; fluxes as outputs.
  !
  ! --------------------------------------------------
  function rte_lw(optical_props, top_at_1, &
                  sources, sfc_emis,       &
                  fluxes,                  &
                  inc_flux, n_gauss_angles, use_2stream, &
                  lw_Ds, flux_up_Jac) result(error_msg)
    class(ty_optical_props_arry), intent(in   ) :: optical_props     ! Set of optical properties as one or more arrays
    logical,                      intent(in   ) :: top_at_1          ! Is the top of the domain at index 1?
                                                                     ! (if not, ordering is bottom-to-top)
    type(ty_source_func_lw),      intent(in   ) :: sources        ! Derived type with Planck source functions
    real(wp), dimension(:,:),     intent(in   ) :: sfc_emis       ! emissivity at surface [] (nband, ncol)
    class(ty_fluxes),             intent(inout) :: fluxes         ! Dervied type for computing spectral integrals from g-point fluxes.
                                                                  ! Default computes broadband fluxes at all levels
                                                                  ! if output arrays are defined. Can be extended per user desires.
    real(wp), dimension(:,:),   &
                target, optional, intent(in   ) :: inc_flux       ! incident flux at domain top [W/m2] (ncol, ngpts)
    integer,            optional, intent(in   ) :: n_gauss_angles ! Number of angles used in Gaussian quadrature (max 3)
                                                                  ! (no-scattering solution)
    logical,            optional, intent(in   ) :: use_2stream    ! When 2-stream parameters (tau/ssa/g) are provided, use 2-stream methods
                                                                  ! Default is to use re-scaled longwave transport
    real(wp), dimension(:,:),   &
                      optional,   intent(in   ) :: lw_Ds          ! User-specifed 1/cos of transport angle per col, g-point
    real(wp), dimension(:,:), target,  &
                      optional,   intent(inout) :: flux_up_Jac    ! surface temperature flux  Jacobian [W/m2/K] (ncol, nlay+1)
    character(len=128)                          :: error_msg      ! If empty, calculation was successful
    ! --------------------------------
    !
    ! Local variables
    !
    integer      :: ncol, nlay, ngpt, nband
    integer      :: n_quad_angs
    logical(wl)  :: using_2stream, do_Jacobians, has_dif_bc, do_broadband
    real(wp), dimension(:,:),   allocatable :: sfc_emis_gpt
    real(wp), dimension(:,:), pointer       :: jacobian
    real(wp), dimension(optical_props%get_ncol(),           &
                        optical_props%get_nlay()+1), target &
                                            :: decoy2D ! Used for optional outputs - needs to be full size.
    real(wp), dimension(1,1,1)              :: empty3D ! Used for optional inputs - size is irrelevant

    real(wp), dimension(:,:,:), pointer     :: gpt_flux_up, gpt_flux_dn, decoy3D
    real(wp), dimension(:,:),   pointer     :: flux_dn_loc, flux_up_loc, flux_net_loc
    real(wp), dimension(:,:),   pointer     :: inc_flux_diffuse
    ! --------------------------------------------------
    !
    ! Weights and angle secants for first order (k=1) Gaussian quadrature.
    !   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    !   after Abramowitz & Stegun 1972, page 921
    !
    integer,  parameter :: max_gauss_pts = 4
    real(wp), parameter,                         &
      dimension(max_gauss_pts, max_gauss_pts) :: &
        gauss_Ds  = RESHAPE([1.66_wp,               0._wp,         0._wp,         0._wp, &  ! Diffusivity angle, not Gaussian angle
                             1.18350343_wp, 2.81649655_wp,         0._wp,         0._wp, &
                             1.09719858_wp, 1.69338507_wp, 4.70941630_wp,         0._wp, &
                             1.06056257_wp, 1.38282560_wp, 2.40148179_wp, 7.15513024_wp], &
                            [max_gauss_pts, max_gauss_pts]),              &
        gauss_wts = RESHAPE([0.5_wp,          0._wp,           0._wp,           0._wp, &
                             0.3180413817_wp, 0.1819586183_wp, 0._wp,           0._wp, &
                             0.2009319137_wp, 0.2292411064_wp, 0.0698269799_wp, 0._wp, &
                             0.1355069134_wp, 0.2034645680_wp, 0.1298475476_wp, 0.0311809710_wp], &
                             [max_gauss_pts, max_gauss_pts])
    ! ------------------------------------------------------------------------------------
    !
    ! Error checking
    !   if inc_flux is present it has the right dimensions, is positive definite
    !
    ! --------------------------------
    ncol  = optical_props%get_ncol()
    nlay  = optical_props%get_nlay()
    ngpt  = optical_props%get_ngpt()
    nband = optical_props%get_nband()
    do_Jacobians = present(flux_up_Jac)
    if(do_Jacobians) then
      jacobian => flux_up_Jac
    else
      jacobian => decoy2D
    end if
    error_msg = ""
    ! ------------------------------------------------------------------------------------
    !
    ! Error checking -- input consistency of sizes and validity of values
    !
    ! --------------------------------

    if(.not. fluxes%are_desired()) &
      error_msg = "rte_lw: no space allocated for fluxes"

    if (do_Jacobians .and. check_extents) then
      if( .not. extents_are(flux_up_Jac, ncol, nlay+1)) &
        error_msg = "rte_lw: flux Jacobian inconsistently sized"
    endif


    if (check_extents) then
      !
      ! Source functions
      !
      if(any([sources%get_ncol(), sources%get_nlay(), sources%get_ngpt()]  /= [ncol, nlay, ngpt])) &
        error_msg = "rte_lw: sources and optical properties inconsistently sized"
      !
      ! Surface emissivity
      !
      if(.not. extents_are(sfc_emis, nband, ncol)) &
        error_msg = "rte_lw: sfc_emis inconsistently sized"
      !
      ! Incident flux, if present
      !
      if(present(inc_flux)) then
        if(.not. extents_are(inc_flux, ncol, ngpt)) &
          error_msg = "rte_lw: inc_flux inconsistently sized"
      end if
    end if

    if(check_values) then
      if(any_vals_outside(sfc_emis, 0._wp, 1._wp)) &
        error_msg = "rte_lw: sfc_emis has values < 0 or > 1"
      if(present(inc_flux)) then
        if(any_vals_less_than(inc_flux, 0._wp)) &
          error_msg = "rte_lw: inc_flux has values < 0"
      end if

      if(present(n_gauss_angles)) then
        if(n_gauss_angles > max_gauss_pts) &
          error_msg = "rte_lw: asking for too many quadrature points for no-scattering calculation"
        if(n_gauss_angles < 1) &
          error_msg = "rte_lw: have to ask for at least one quadrature point for no-scattering calculation"
      end if
    end if
    if(len_trim(error_msg) > 0) return

    !
    ! Number of quadrature points for no-scattering calculation
    !
    n_quad_angs = 1
    if(present(n_gauss_angles)) n_quad_angs = n_gauss_angles
    !
    ! Optionally - use 2-stream methods when low-order scattering properties are provided?
    !
    using_2stream = .false.
    if(present(use_2stream)) using_2stream = use_2stream

    !
    ! Checking that optional arguements are consistent with one another and with optical properties
    !
    select type (optical_props)
      class is (ty_optical_props_1scl)
        if (using_2stream) &
          error_msg = "rte_lw: can't use two-stream methods with only absorption optical depth"
        if (present(lw_Ds)) then
          if(.not. extents_are(lw_Ds, ncol, ngpt)) &
            error_msg = "rte_lw: lw_Ds inconsistently sized"
          if(any_vals_less_than(lw_Ds, 1._wp)) &
            error_msg = "rte_lw: one or more values of lw_Ds < 1."
          if(n_quad_angs /= 1) &
            error_msg = "rte_lw: providing lw_Ds incompatible with specifying n_gauss_angles"
        end if
      class is (ty_optical_props_2str)
        if (present(lw_Ds)) &
          error_msg = "rte_lw: lw_Ds not valid when providing scattering optical properties"
        if (using_2stream .and. n_quad_angs /= 1) &
          error_msg = "rte_lw: using_2stream=true incompatible with specifying n_gauss_angles"
        if (using_2stream .and. do_Jacobians) &
          error_msg = "rte_lw: can't provide Jacobian of fluxes w.r.t surface temperature with 2-stream"
      class default
        error_msg =  "rte_lw: lw_solver(...ty_optical_props_nstr...) not yet implemented"
    end select
    if(len_trim(error_msg) > 0) return

    !
    ! Ensure values of tau, ssa, and g are reasonable if using scattering
    !
    !$acc enter data copyin(optical_props)
    if(check_values)  then
      error_msg =  optical_props%validate()
      if(len_trim(error_msg) > 0) then
        if(len_trim(optical_props%get_name()) > 0) &
          error_msg = trim(optical_props%get_name()) // ': ' // trim(error_msg)
        return
      end if
    end if
    ! ------------------------------------------------------------------------------------
    !
    !    Lower boundary condition -- expand surface emissivity by band to gpoints
    !
    allocate(sfc_emis_gpt(ncol,         ngpt))
    !$acc        enter data create(   sfc_emis_gpt)
    !$omp target enter data map(alloc:sfc_emis_gpt)
    call expand_and_transpose(optical_props, sfc_emis, sfc_emis_gpt)

    !$acc        enter data create(   flux_up_Jac) if(do_Jacobians)
    !$omp target enter data map(alloc:flux_up_Jac) if(do_Jacobians)

    select type(fluxes)
      type is (ty_fluxes_broadband)
        do_broadband = .true._wl
        !
        ! Broadband fluxes class has three possible outputs; allocate memory for local use
        !   if one or more haven't been requested
        !
        if(associated(fluxes%flux_up)) then
          flux_up_loc => fluxes%flux_up
        else
          allocate(flux_up_loc(ncol, nlay+1))
          !$acc        enter data create(   flux_up_loc)
          !$omp target enter data map(alloc:flux_up_loc)
        end if
        if(associated(fluxes%flux_dn)) then
          flux_dn_loc => fluxes%flux_dn
        else
          allocate(flux_dn_loc(ncol, nlay+1))
          !$acc        enter data create(   flux_dn_loc)
          !$omp target enter data map(alloc:flux_dn_loc)
        end if
      class default
        !
        ! If broadband integrals aren't being computed, allocate working space
        !   and decoy addresses for spectrally-integrated fields
        !
        do_broadband = .false._wl
        !$acc        enter data create(   decoy2D)
        !$omp target enter data map(alloc:decoy2D)
        flux_up_loc  => decoy2D
        flux_dn_loc  => decoy2D
    end select
    ! Move the below to the default case later
    allocate(gpt_flux_up (ncol,nlay+1,ngpt), &
             gpt_flux_dn (ncol,nlay+1,ngpt))
    !$acc        enter data create(   gpt_flux_up, gpt_flux_dn)
    !$omp target enter data map(alloc:gpt_flux_up, gpt_flux_dn)

    !
    !   Upper boundary condition -  use values in optional arg or be set to 0
    !
    if (present(inc_flux)) then
      !  inc_flux_dif might have been copied-in before, when checking extents/values,
      !  in which case no copying will occur here
      !$acc        enter data copyin(inc_flux_dif)
      !$omp target enter data map(to:inc_flux_dif)
      inc_flux_diffuse => inc_flux
    else
      allocate(inc_flux_diffuse(ncol, ngpt))
      !$acc        enter data create(   inc_flux_diffuse)
      !$omp target enter data map(alloc:inc_flux_diffuse)
      call zero_array(ncol, ngpt, inc_flux_diffuse)
    end if


    !
    ! Compute the radiative transfer...
    !
    select type (optical_props)
      class is (ty_optical_props_1scl)
        !
        ! No scattering two-stream calculation
        !
        !$acc        enter data copyin(optical_props%tau)
        !$omp target enter data map(to:optical_props%tau)
        error_msg =  optical_props%validate()
        if(len_trim(error_msg) > 0) return

        if (present(lw_Ds)) then
          call lw_solver_noscat(ncol, nlay, ngpt, &
                                logical(top_at_1, wl), &
                                lw_Ds, gauss_wts(1,1), &
                                optical_props%tau, &
                                sources%lay_source, sources%lev_source_inc, sources%lev_source_dec, &
                                sfc_emis_gpt, sources%sfc_source,  &
                                inc_flux_diffuse,                  &
                                gpt_flux_up, gpt_flux_dn,          &
                                do_broadband, flux_up_loc, flux_dn_loc,     &
                                logical(do_Jacobians, wl), sources%sfc_source_Jac, jacobian, &
                                logical(.false., wl),  empty3D, empty3D)
        else
          call lw_solver_noscat_GaussQuad(ncol, nlay, ngpt,                 &
                                logical(top_at_1, wl), n_quad_angs,         &
                                gauss_Ds(1:n_quad_angs,n_quad_angs),        &
                                gauss_wts(1:n_quad_angs,n_quad_angs),       &
                                optical_props%tau,                 &
                                sources%lay_source, sources%lev_source_inc, &
                                sources%lev_source_dec,            &
                                sfc_emis_gpt, sources%sfc_source,  &
                                inc_flux_diffuse,                  &
                                gpt_flux_up, gpt_flux_dn,          &
                                do_broadband, flux_up_loc, flux_dn_loc,     &
                                logical(do_Jacobians, wl), sources%sfc_source_Jac, jacobian, &
                                logical(.false., wl),  empty3D, empty3D)
        end if
        !$acc        exit data delete(     optical_props%tau)
        !$omp target exit data map(release:optical_props%tau)
      class is (ty_optical_props_2str)
        if (using_2stream) then
          !
          ! two-stream calculation with scattering
          !
          !$acc        enter data copyin(optical_props%tau, optical_props%ssa, optical_props%g)
          !$omp target enter data map(to:optical_props%tau, optical_props%ssa, optical_props%g)
          error_msg =  optical_props%validate()
          if(len_trim(error_msg) > 0) return
          call lw_solver_2stream(ncol, nlay, ngpt, logical(top_at_1, wl), &
                                 optical_props%tau, optical_props%ssa, optical_props%g,              &
                                 sources%lay_source, sources%lev_source_inc, sources%lev_source_dec, &
                                 sfc_emis_gpt, sources%sfc_source,       &
                                 inc_flux_diffuse,                       &
                                 gpt_flux_up, gpt_flux_dn)
          !$acc        exit data delete(     optical_props%tau, optical_props%ssa, optical_props%g)
          !$omp target exit data map(release:optical_props%tau, optical_props%ssa, optical_props%g)
        else
          !
          ! Re-scaled solution to account for scattering
          !
          !$acc        enter data copyin(optical_props%tau, optical_props%ssa, optical_props%g)
          !$omp target enter data map(to:optical_props%tau, optical_props%ssa, optical_props%g)
          call lw_solver_noscat_GaussQuad(ncol, nlay, ngpt,                 &
                                logical(top_at_1, wl), n_quad_angs,         &
                                gauss_Ds(1:n_quad_angs,n_quad_angs),        &
                                gauss_wts(1:n_quad_angs,n_quad_angs),       &
                                optical_props%tau,                 &
                                sources%lay_source, sources%lev_source_inc, &
                                sources%lev_source_dec,            &
                                sfc_emis_gpt, sources%sfc_source,  &
                                inc_flux_diffuse,                  &
                                gpt_flux_up, gpt_flux_dn,          &
                                do_broadband, flux_up_loc, flux_dn_loc,     &
                                logical(do_Jacobians, wl), sources%sfc_source_Jac, jacobian, &
                                logical(.true., wl),  optical_props%ssa, optical_props%g)
          !$acc        exit data delete(     optical_props%tau, optical_props%ssa, optical_props%g)
          !$omp target exit data map(release:optical_props%tau, optical_props%ssa, optical_props%g)
        endif
      class is (ty_optical_props_nstr)
        !
        ! n-stream calculation
        !
        error_msg = 'lw_solver(...ty_optical_props_nstr...) not yet implemented'
    end select

    if (error_msg /= '') return

    select type(fluxes)
      !
      ! Tidy up memory for broadband fluxes on GPUs
      !
      type is (ty_fluxes_broadband)
        !$acc        exit data delete(     decoy3D)
        !$omp target exit data map(release:decoy3D)
        if(associated(fluxes%flux_net)) then
          !
          ! Make this OpenACC/MP friendly
          !
          fluxes%flux_net(:,:) = flux_dn_loc(:,:) - flux_up_loc(:,:)
        end if
        if(associated(fluxes%flux_up)) then
          !$acc        exit data copyout( flux_up_loc)
          !$omp target exit data map(from:flux_up_loc)
        else
          !$acc        exit data delete(     flux_up_loc)
          !$omp target exit data map(release:flux_up_loc)
        end if
        if(associated(fluxes%flux_dn)) then
          !$acc        exit data copyout( flux_dn_loc)
          !$omp target exit data map(from:flux_dn_loc)
        else
          !$acc        exit data delete(     flux_dn_loc)
          !$omp target exit data map(release:flux_dn_loc)
        end if
      class default
        !
        ! ...or reduce spectral fluxes to desired output quantities
        !
        error_msg = fluxes%reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)
        !$acc        exit data delete(     gpt_flux_up, gpt_flux_dn, gpt_flux_dir, decoy2D)
        !$omp target exit data map(release:gpt_flux_up, gpt_flux_dn, gpt_flux_dir, decoy2D)
    end select
    !$acc        exit data delete(     sfc_emis_gpt)
    !$omp target exit data map(release:sfc_emis_gpt)
    !$acc        exit data delete(optical_props)
    !$acc        exit data copyout( flux_up_Jac) if(do_Jacobians)
    !$omp target exit data map(from:flux_up_Jac) if(do_Jacobians)

  end function rte_lw
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
  !
  subroutine expand_and_transpose(ops,arr_in,arr_out)
    class(ty_optical_props),  intent(in ) :: ops
    real(wp), dimension(:,:), intent(in ) :: arr_in  ! (nband, ncol)
    real(wp), dimension(:,:), intent(out) :: arr_out ! (ncol, igpt)
    ! -------------
    integer :: ncol, nband, ngpt
    integer :: icol, iband, igpt
    integer, dimension(2,ops%get_nband()) :: limits

    ncol  = size(arr_in, 2)
    nband = ops%get_nband()
    ngpt  = ops%get_ngpt()
    limits = ops%get_band_lims_gpoint()
    !$acc parallel loop collapse(2) copyin(arr_in, limits)
    !$omp target teams distribute parallel do simd collapse(2) map(to:arr_in, limits)
    do iband = 1, nband
      do icol = 1, ncol
        do igpt = limits(1, iband), limits(2, iband)
          arr_out(icol, igpt) = arr_in(iband,icol)
        end do
      end do
    end do

  end subroutine expand_and_transpose
  !--------------------------------------------------------------------------------------------------------------------
end module mo_rte_lw
