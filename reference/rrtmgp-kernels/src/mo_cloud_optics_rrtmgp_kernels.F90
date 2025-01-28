module mo_cloud_optics_rrtmgp_kernels
  use mo_rte_kind,      only : wp, wl
  implicit none
  private
  public :: compute_cld_from_table
  interface
    !---------------------------------------------------------------------------
    !
    ! Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
    !   elements starting at "offset." The table's second dimension is spectral -
    !   ngpt == nbnd if the tables are defined on bands.
    ! Returns 0 where the mask is false.
    ! We could also try gather/scatter for efficiency
    !
    subroutine compute_cld_from_table(ncol, nlay, ngpt, mask, lwp, re, &
                                      nsteps, step_size, offset,       &
                                      tau_table, ssa_table, asy_table, &
                                      tau, taussa, taussag) bind(C, name="rrtmgp_compute_cld_from_table")
      use mo_rte_kind, only : wp, wl
      integer,                                intent(in) :: ncol, nlay, ngpt, nsteps
      logical(wl), dimension(ncol,nlay),      intent(in) :: mask
      real(wp),    dimension(ncol,nlay),      intent(in) :: lwp, re
      real(wp),                               intent(in) :: step_size, offset
      real(wp),    dimension(nsteps,   ngpt), intent(in) :: tau_table, ssa_table, asy_table
      real(wp),    dimension(ncol,nlay,ngpt)             :: tau, taussa, taussag
    end subroutine compute_cld_from_table
  end interface
end module mo_cloud_optics_rrtmgp_kernels
