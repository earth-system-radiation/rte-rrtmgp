! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Copyright 2024-,  Atmospheric and Environmental Research,
!    Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
module mo_cloud_optics_rrtmgp_kernels
  use mo_rte_kind,      only : wp, wl
  implicit none
  private
  public :: compute_cld_from_table
contains
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
    integer,                                 intent(in) :: ncol, nlay, ngpt, nsteps
    logical(wl), dimension(ncol,nlay),       intent(in) :: mask
    real(wp),    dimension(ncol,nlay),       intent(in) :: lwp, re
    real(wp),                                intent(in) :: step_size, offset
    real(wp),    dimension(nsteps,   ngpt), intent(in) :: tau_table, ssa_table, asy_table
    real(wp),    dimension(ncol,nlay,ngpt), intent(out):: tau, taussa, taussag
    ! ---------------------------
    integer  :: icol, ilay, igpt
    integer  :: index
    real(wp) :: fint
    real(wp) :: t, ts  ! tau, tau*ssa, tau*ssa*g
    ! ---------------------------
    !$acc parallel loop gang vector default(present) collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do igpt = 1, ngpt
      do ilay = 1,nlay
        do icol = 1, ncol
          if(mask(icol,ilay)) then
            index = min(floor((re(icol,ilay) - offset)/step_size)+1, nsteps-1)
            fint = (re(icol,ilay) - offset)/step_size - (index-1)
            t   = lwp(icol,ilay) * &
                  (tau_table(index,  igpt) + fint * (tau_table(index+1,igpt) - tau_table(index,igpt)))
            ts  = t              * &
                  (ssa_table(index,  igpt) + fint * (ssa_table(index+1,igpt) - ssa_table(index,igpt)))
            taussag(icol,ilay,igpt) =  &
                  ts             * &
                  (asy_table(index,  igpt) + fint * (asy_table(index+1,igpt) - asy_table(index,igpt)))
            taussa (icol,ilay,igpt) = ts
            tau    (icol,ilay,igpt) = t
          else
            tau    (icol,ilay,igpt) = 0._wp
            taussa (icol,ilay,igpt) = 0._wp
            taussag(icol,ilay,igpt) = 0._wp
          end if
        end do
      end do
    end do
  end subroutine compute_cld_from_table

end module mo_cloud_optics_rrtmgp_kernels
