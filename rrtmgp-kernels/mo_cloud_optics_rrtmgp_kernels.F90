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
  !   elements starting at "offset." The table's second dimension is band.
  ! Returns 0 where the mask is false.
  ! We could also try gather/scatter for efficiency
  !
  subroutine compute_cld_from_table(ncol, nlay, nspec, mask, lwp, re, &
                                    nsteps, step_size, offset,       &
                                    tau_table, ssa_table, asy_table, &
                                    tau, taussa, taussag) bind(C, name="rrtmgp_compute_cld_from_table")
    integer,                                 intent(in) :: ncol, nlay, nspec, nsteps
    logical(wl), dimension(ncol,nlay),       intent(in) :: mask
    real(wp),    dimension(ncol,nlay),       intent(in) :: lwp, re
    real(wp),                                intent(in) :: step_size, offset
    real(wp),    dimension(nsteps,   nspec), intent(in) :: tau_table, ssa_table, asy_table
    real(wp),    dimension(ncol,nlay,nspec), intent(out):: tau, taussa, taussag
    ! ---------------------------
    integer  :: icol, ilay, ispec
    integer  :: index
    real(wp) :: fint
    real(wp) :: t, ts  ! tau, tau*ssa, tau*ssa*g
    ! ---------------------------
    !$acc parallel loop gang vector default(present) collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do ispec = 1, nspec
      do ilay = 1,nlay
        do icol = 1, ncol
          if(mask(icol,ilay)) then
            index = min(floor((re(icol,ilay) - offset)/step_size)+1, nsteps-1)
            fint = (re(icol,ilay) - offset)/step_size - (index-1)
            t   = lwp(icol,ilay) * &
                  (tau_table(index,  ispec) + fint * (tau_table(index+1,ispec) - tau_table(index,ispec)))
            ts  = t              * &
                  (ssa_table(index,  ispec) + fint * (ssa_table(index+1,ispec) - ssa_table(index,ispec)))
            taussag(icol,ilay,ispec) =  &
                  ts             * &
                  (asy_table(index,  ispec) + fint * (asy_table(index+1,ispec) - asy_table(index,ispec)))
            taussa (icol,ilay,ispec) = ts
            tau    (icol,ilay,ispec) = t
          else
            tau    (icol,ilay,ispec) = 0._wp
            taussa (icol,ilay,ispec) = 0._wp
            taussag(icol,ilay,ispec) = 0._wp
          end if
        end do
      end do
    end do
  end subroutine compute_cld_from_table

end module mo_cloud_optics_rrtmgp_kernels
