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
! Encapsulate optical properties defined on a spectral grid of N bands.
!   The bands are described by their limiting wavenumbers. They need not be contiguous or complete.
!   A band may contain more than one spectral sub-point (g-point) in which case a mapping must be supplied.
!   A name may be provided and will be prepended to error messages.
!   The base class (ty_optical_props) encapsulates only this spectral discretization and must be initialized
!      with the spectral information before use.
!
!   Optical properties may be represented as arrays with dimensions ncol, nlay, ngpt
!   (abstract class ty_optical_props_arry).
!   The type holds arrays depending on how much information is needed
!   There are a possibilites
!      ty_optical_props_1rescl  holds absorption optical depth tau, used in calculations accounting for extinction and emission, as well as scaling coefficient
!
! Optical properties can increment or "add themselves to" a set of properties represented with arrays
!   as long as both sets have the same underlying band structure. Properties defined by band
!   may be added to properties defined by g-point; the same value is assumed for all g-points with each band.
!
! Subsets of optical properties held as arrays may be extracted along the column dimension.
!
! -------------------------------------------------------------------------------------------------
module mo_optical_props_add
  use mo_rte_kind,              only: wp
  use mo_optical_props,         only: ty_optical_props_1scl
 
  ! Implemented based on the paper
  ! Tang G, P Yang, GW Kattawar, X Huang, EJ Mlawer, BA Baum, MD King, 2018: Improvement of 
  ! the Simulation of Cloud Longwave Scattering in Broadband Radiative Transfer Models, 
  ! Journal of the Atmospheric Sciences 75 (7), 2217-2233
  ! https://doi.org/10.1175/JAS-D-18-0014.1
  !
  type, extends(ty_optical_props_1scl) :: ty_optical_props_1rescl
    real(wp), dimension(:,:,:),   allocatable :: scaling ! adjuxtment factor in Tang approximation  ssa b/[1 - ssa(1-b)] (ncol, nlay, ngpt) where b = (1 - g) / 2)
  contains
      procedure, public  :: create_1rescl => create_1rescl_from_2str
      procedure, public  :: validate => validate_1rescl
      procedure, public  :: get_subset => subset_1rescl_range

      procedure, private :: alloc_only_1rescl
      procedure, private :: init_and_alloc_1rescl
      procedure, private :: copy_and_alloc_1rescl
      generic,   public  :: alloc_1rescl => alloc_only_1rescl, init_and_alloc_1rescl, copy_and_alloc_1rescl
  end type  ty_optical_props_1rescl


contains  

  ! ------------------------------------------------------------------------------------------
  function create_1rescl_from_2str(this, dat_2str) result(err_message)
    use mo_optical_props,         only: ty_optical_props_2str
    class(ty_optical_props_1rescl), intent(inout) :: this
    class(ty_optical_props_2str),   intent(in   ) :: dat_2str
    character(128)                                :: err_message

    real(wp) :: wf
    integer  :: icol, ilay, igpt
    integer :: ncol, nlay, ngpt
    ! --------------------------------
    ncol = dat_2str%get_ncol()
    nlay = dat_2str%get_nlay()
    ngpt = dat_2str%get_ngpt()
    err_message = ""

    err_message = this%alloc_1rescl(ncol, nlay, dat_2str, dat_2str%get_name())
    if ( err_message /= '' ) return

    call scalingTang(ncol, nlay, ngpt, this%tau, this%scaling, dat_2str%tau, dat_2str%ssa, dat_2str%g)

  end function create_1rescl_from_2str  
  
  ! can be moved to proper kernel
  pure subroutine scalingTang(ncol, nlay, ngpt, tauLoc, scaling, tau, ssa, g)
    integer ,                              intent(in)    :: ncol
    integer ,                              intent(in)    :: nlay
    integer ,                              intent(in)    :: ngpt
    real(wp), dimension(ncol, nlay, ngpt), intent(in)    :: tau
    real(wp), dimension(ncol, nlay, ngpt), intent(in)    :: ssa
    real(wp), dimension(ncol, nlay, ngpt), intent(in)    :: g

    real(wp), dimension(ncol, nlay, ngpt), intent(inout) :: tauLoc
    real(wp), dimension(ncol, nlay, ngpt), intent(inout) :: scaling

    integer  :: icol, ilay, igpt
    real(wp) :: xx, ssal, yy
    !$acc enter data copyin(tau, ssa, g)
    !$acc enter data create(tauLoc, scaling)

    !$acc parallel loop collapse(3)
    do igpt=1,ngpt
      do ilay=1,nlay
        do icol=1,ncol
          ssal = ssa(icol, ilay, igpt)
          xx = ssal*(1._wp - g(icol, ilay, igpt)) / 2._wp
          yy = (1._wp - ssal + xx )
          ! Eq.15 of the paper
          tauLoc(icol, ilay, igpt) = yy * tau(icol, ilay, igpt)
          ! 
          ! here ssa is used to store parameter wb/[1-w(1-b)] of Eq.21 of the Tang's paper
          ! actually it is in line of parameter rescaling defined in Eq.7
          if (yy > epsilon(1._wp)) then
            scaling(icol, ilay, igpt) = xx / yy
          else
            scaling(icol, ilay, igpt) = 1._wp
          endif
        enddo
      enddo
    enddo
    !$acc exit data copyout(tauLoc, scaling)
    !$acc exit data delete(tau, ssa, g)
  end subroutine scalingTang
  ! ------------------------------------------------------------------------------------------
  !
  ! --- Validation
  !
  ! ------------------------------------------------------------------------------------------
  function validate_1rescl(this) result(err_message)
    use mo_rte_util_array, only: any_vals_less_than, any_vals_outside
    class(ty_optical_props_1rescl), intent(in) :: this
    character(len=128)                         :: err_message

    err_message = ''
    if(.not. allocated(this%tau)) then
      err_message = " tau not allocated"
    end if
    if(.not. allocated(this%scaling)) then
      err_message = trim(err_message) // " scaling not allocated"
    end if
    if(any_vals_less_than(this%tau, 0._wp)) &
      err_message = trim(err_message) // " tau values out of range"

    if(any_vals_outside(this%scaling, 0._wp, 1._wp)) &
      err_message = trim(err_message) // " scaling values out of range"

    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
      err_message = 'validate_1rescl : ' // trim(err_message)

  end function validate_1rescl

! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: initialization, allocation, and finalization
  !    Initialization and allocation can be combined by supplying either
  !
  ! ------------------------------------------------------------------------------------------
  !
  ! Straight allocation routines
  !
  function alloc_only_1rescl(this, ncol, nlay) result(err_message)
    class(ty_optical_props_1rescl) :: this
    integer,          intent(in) :: ncol, nlay
    character(len=128)           :: err_message

    err_message = ""
    if(.not. this%is_initialized()) then
      err_message = "alloc_only_1rescl%optical_props%alloc: spectral discretization hasn't been provided"
      return
    end if
    if(any([ncol, nlay] <= 0)) then
      err_message = "alloc_only_1rescl%optical_props%alloc: must provide positive extents for ncol, nlay"
    else
      if (allocated(this%tau    )) deallocate(this%tau)
      if (allocated(this%scaling)) deallocate(this%scaling)
      allocate(this%tau    (ncol, nlay, this%get_ngpt()))
      allocate(this%scaling(ncol, nlay, this%get_ngpt()))
    end if
  end function alloc_only_1rescl

! ------------------------------------------------------------------------------------------
  !
  ! Combined allocation/initialization routines
  !
  ! ------------------------------------------------------------------------------------------
  !
  ! Initialization by specifying band limits and possibly g-point/band mapping
  !
  ! ---------------------------------------------------------------------------
  function init_and_alloc_1rescl(this, ncol, nlay, band_lims_wvn, band_lims_gpt, name) result(err_message)
    class(ty_optical_props_1rescl)             :: this
    integer,                      intent(in) :: ncol, nlay
    real(wp), dimension(:,:),     intent(in) :: band_lims_wvn
    integer,  dimension(:,:), &
                      optional,   intent(in) :: band_lims_gpt
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message

    err_message = this%ty_optical_props%init(band_lims_wvn, &
                                             band_lims_gpt, name)
    if(err_message /= "") return
    err_message = this%alloc_only_1rescl(ncol, nlay)
  end function init_and_alloc_1rescl

  !-------------------------------------------------------------------------------------------------
  !
  ! Initialization from an existing spectral discretization/ty_optical_props
  !
  !-------------------------------------------------------------------------------------------------
  function copy_and_alloc_1rescl(this, ncol, nlay, spectral_desc, name) result(err_message)
    use mo_optical_props, only: ty_optical_props
    class(ty_optical_props_1rescl)           :: this
    integer,                      intent(in) :: ncol, nlay
    class(ty_optical_props     ), intent(in) :: spectral_desc
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message

    err_message = ""
    if(this%ty_optical_props%is_initialized()) call this%ty_optical_props%finalize()
    err_message = this%init_and_alloc_1rescl(ncol, nlay, &
                                             spectral_desc%get_band_lims_wavenumber(), &
                                             spectral_desc%get_band_lims_gpoint(), name)

  end function copy_and_alloc_1rescl

! ------------------------------------------------------------------------------------------
  !
  !  Routines for array classes: subsetting of optical properties arrays along x (col) direction
  !
  ! Allocate class, then arrays; copy. Could probably be more efficient if
  !   classes used pointers internally.
  !
  ! This set takes start position and number as scalars
  !
  ! ------------------------------------------------------------------------------------------

  function subset_1rescl_range(full, start, n, subset) result(err_message)
    use mo_optical_props_kernels, only: extract_subset
    use mo_optical_props,         only: ty_optical_props_arry, &
                                        ty_optical_props_2str, &
                                        ty_optical_props_nstr

    class(ty_optical_props_1rescl), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end if
    ncol = full%get_ncol()
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full)
    ! Seems like the deallocation statements should be needed under Fortran 2003
    !   but Intel compiler doesn't run without them
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return

      class is (ty_optical_props_1rescl)
        err_message = subset%alloc_1rescl(n, nlay)
        if(err_message /= "") return
        call extract_subset(ncol, nlay, ngpt, full%scaling, start, start+n-1, subset%scaling)

      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = 0._wp
        subset%g  (1:n,:,:) = 0._wp

      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = subset%get_nmom()
          deallocate(subset%p  )
        else
          nmom = 1
        end if
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = 0._wp
        subset%p(:,1:n,:,:) = 0._wp
    end select
    call extract_subset(ncol, nlay, ngpt, full%tau, start, start+n-1, subset%tau)

  end function subset_1rescl_range
end module mo_optical_props_add
