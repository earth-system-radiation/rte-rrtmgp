module mo_cloud_optics_rrtmgp_kernels
  use mo_rte_kind,      only : wp, wl
  implicit none
  private
  public :: compute_all_from_table, compute_all_from_pade
contains
  interface 
    !---------------------------------------------------------------------------
    !
    ! Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
    !   elements starting at "offset." The table's second dimension is band.
    ! Returns 0 where the mask is false.
    ! We could also try gather/scatter for efficiency
    !
    subroutine compute_all_from_table(ncol, nlay, nbnd, mask, lwp, re, &
                                      nsteps, step_size, offset,       &
                                      tau_table, ssa_table, asy_table, &
                                      tau, taussa, taussag) bind(C, name="rrtmgp_compute_all_from_table")
      integer,                                intent(in) :: ncol, nlay, nbnd, nsteps
      logical(wl), dimension(ncol,nlay),      intent(in) :: mask
      real(wp),    dimension(ncol,nlay),      intent(in) :: lwp, re
      real(wp),                               intent(in) :: step_size, offset
      real(wp),    dimension(nsteps,   nbnd), intent(in) :: tau_table, ssa_table, asy_table
      real(wp),    dimension(ncol,nlay,nbnd)             :: tau, taussa, taussag
    end subroutine compute_all_from_table

    !---------------------------------------------------------------------------
    !
    ! Pade functions
    !
    !---------------------------------------------------------------------------
    subroutine compute_all_from_pade(ncol, nlay, nbnd, nsizes, &
                                     mask, lwp, re,            &
                                     m_ext, n_ext, re_bounds_ext, coeffs_ext, &
                                     m_ssa, n_ssa, re_bounds_ssa, coeffs_ssa, &
                                     m_asy, n_asy, re_bounds_asy, coeffs_asy, &
                                     tau, taussa, taussag) bind(C, name="rrtmgp_compute_all_from_pade")
      integer,                        intent(in) :: ncol, nlay, nbnd, nsizes
      logical(wl),  &
                dimension(ncol,nlay), intent(in) :: mask
      real(wp), dimension(ncol,nlay), intent(in) :: lwp, re
      real(wp), dimension(nsizes+1),  intent(in) :: re_bounds_ext, re_bounds_ssa, re_bounds_asy
      integer,                        intent(in) :: m_ext, n_ext
      real(wp), dimension(nbnd,nsizes,0:m_ext+n_ext), &
                                      intent(in) :: coeffs_ext
      integer,                        intent(in) :: m_ssa, n_ssa
      real(wp), dimension(nbnd,nsizes,0:m_ssa+n_ssa), &
                                      intent(in) :: coeffs_ssa
      integer,                        intent(in) :: m_asy, n_asy
      real(wp), dimension(nbnd,nsizes,0:m_asy+n_asy), &
                                      intent(in) :: coeffs_asy
      real(wp), dimension(ncol,nlay,nbnd)        :: tau, taussa, taussag
    end subroutine compute_all_from_pade
  end interface
end module module mo_cloud_optics_rrtmgp_kernels

