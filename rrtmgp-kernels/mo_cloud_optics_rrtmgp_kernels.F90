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
  public :: compute_cld_from_table, compute_cld_from_pade
  interface pade_eval
    module procedure pade_eval_nbnd, pade_eval_1
  end interface pade_eval
contains
  !---------------------------------------------------------------------------
  !
  ! Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
  !   elements starting at "offset." The table's second dimension is band.
  ! Returns 0 where the mask is false.
  ! We could also try gather/scatter for efficiency
  !
  subroutine compute_cld_from_table(ncol, nlay, nbnd, mask, lwp, re, &
                                    nsteps, step_size, offset,       &
                                    tau_table, ssa_table, asy_table, &
                                    tau, taussa, taussag) bind(C, name="rrtmgp_compute_cld_from_table")
    integer,                                intent(in) :: ncol, nlay, nbnd, nsteps
    logical(wl), dimension(ncol,nlay),      intent(in) :: mask
    real(wp),    dimension(ncol,nlay),      intent(in) :: lwp, re
    real(wp),                               intent(in) :: step_size, offset
    real(wp),    dimension(nsteps,   nbnd), intent(in) :: tau_table, ssa_table, asy_table
    real(wp),    dimension(ncol,nlay,nbnd)             :: tau, taussa, taussag
    ! ---------------------------
    integer  :: icol, ilay, ibnd
    integer  :: index
    real(wp) :: fint
    real(wp) :: t, ts  ! tau, tau*ssa, tau*ssa*g
    ! ---------------------------
    !$acc parallel loop gang vector default(present) collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do ibnd = 1, nbnd
      do ilay = 1,nlay
        do icol = 1, ncol
          if(mask(icol,ilay)) then
            index = min(floor((re(icol,ilay) - offset)/step_size)+1, nsteps-1)
            fint = (re(icol,ilay) - offset)/step_size - (index-1)
            t   = lwp(icol,ilay) * &
                  (tau_table(index,  ibnd) + fint * (tau_table(index+1,ibnd) - tau_table(index,ibnd)))
            ts  = t              * &
                  (ssa_table(index,  ibnd) + fint * (ssa_table(index+1,ibnd) - ssa_table(index,ibnd)))
            taussag(icol,ilay,ibnd) =  &
                  ts             * &
                  (asy_table(index,  ibnd) + fint * (asy_table(index+1,ibnd) - asy_table(index,ibnd)))
            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          else
            tau    (icol,ilay,ibnd) = 0._wp
            taussa (icol,ilay,ibnd) = 0._wp
            taussag(icol,ilay,ibnd) = 0._wp
          end if
        end do
      end do
    end do
  end subroutine compute_cld_from_table
  !---------------------------------------------------------------------------
  !
  ! Pade functions
  !
  !---------------------------------------------------------------------------
  subroutine compute_cld_from_pade(ncol, nlay, nbnd, nsizes, &
                                   mask, lwp, re,            &
                                   m_ext, n_ext, re_bounds_ext, coeffs_ext, &
                                   m_ssa, n_ssa, re_bounds_ssa, coeffs_ssa, &
                                   m_asy, n_asy, re_bounds_asy, coeffs_asy, &
                                   tau, taussa, taussag) bind(C, name="rrtmgp_compute_cld_from_pade")
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
    ! ---------------------------
    integer  :: icol, ilay, ibnd, irad
    real(wp) :: t, ts

    !$acc parallel loop gang vector default(present) collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do ibnd = 1, nbnd
      do ilay = 1, nlay
        do icol = 1, ncol
          if(mask(icol,ilay)) then
            !
            ! Finds index into size regime table
            ! This works only if there are precisely three size regimes (four bounds) and it's
            !   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
            !
            irad = min(floor((re(icol,ilay) - re_bounds_ext(2))/re_bounds_ext(3))+2, 3)
            t   = lwp(icol,ilay) *     &
                  pade_eval(ibnd, nbnd, nsizes, m_ext, n_ext, irad, re(icol,ilay), coeffs_ext)

            irad = min(floor((re(icol,ilay) - re_bounds_ssa(2))/re_bounds_ssa(3))+2, 3)
            ! Pade approximants for co-albedo can sometimes be negative
            ts  = t              * (1._wp - max(0._wp, &
                  pade_eval(ibnd, nbnd, nsizes, m_ssa, n_ssa, irad, re(icol,ilay), coeffs_ssa)))

            irad = min(floor((re(icol,ilay) - re_bounds_asy(2))/re_bounds_asy(3))+2, 3)
            taussag(icol,ilay,ibnd) =  &
                  ts             *     &
                  pade_eval(ibnd, nbnd, nsizes, m_asy, n_asy, irad, re(icol,ilay), coeffs_asy)

            taussa (icol,ilay,ibnd) = ts
            tau    (icol,ilay,ibnd) = t
          else
            tau    (icol,ilay,ibnd) = 0._wp
            taussa (icol,ilay,ibnd) = 0._wp
            taussag(icol,ilay,ibnd) = 0._wp
          end if
        end do
      end do
    end do

  end subroutine compute_cld_from_pade
  !---------------------------------------------------------------------------
  !
  ! Evaluate Pade approximant of order [m/n]
  !
  function pade_eval_nbnd(nbnd, nrads, m, n, irad, re, pade_coeffs)
    integer,                intent(in) :: nbnd, nrads, m, n, irad
    real(wp), dimension(nbnd, nrads, 0:m+n), &
                            intent(in) :: pade_coeffs
    real(wp),               intent(in) :: re
    real(wp), dimension(nbnd)          :: pade_eval_nbnd

    integer :: iband
    real(wp) :: numer, denom
    integer  :: i

    do iband = 1, nbnd
      denom = pade_coeffs(iband,irad,n+m)
      do i = n-1+m, 1+m, -1
        denom = pade_coeffs(iband,irad,i)+re*denom
      end do
      denom =  1._wp                     +re*denom

      numer = pade_coeffs(iband,irad,m)
      do i = m-1, 1, -1
        numer = pade_coeffs(iband,irad,i)+re*numer
      end do
      numer = pade_coeffs(iband,irad,0)  +re*numer

      pade_eval_nbnd(iband) = numer/denom
    end do
  end function pade_eval_nbnd
  !---------------------------------------------------------------------------
  !
  ! Evaluate Pade approximant of order [m/n]
  !
  function pade_eval_1(iband, nbnd, nrads, m, n, irad, re, pade_coeffs)
    !$acc routine seq
    !$omp declare target
    !
    integer,                intent(in) :: iband, nbnd, nrads, m, n, irad
    real(wp), dimension(nbnd, nrads, 0:m+n), &
                            intent(in) :: pade_coeffs
    real(wp),               intent(in) :: re
    real(wp)                           :: pade_eval_1

    real(wp) :: numer, denom
    integer  :: i

    denom = pade_coeffs(iband,irad,n+m)
    do i = n-1+m, 1+m, -1
      denom = pade_coeffs(iband,irad,i)+re*denom
    end do
    denom =  1._wp                     +re*denom

    numer = pade_coeffs(iband,irad,m)
    do i = m-1, 1, -1
      numer = pade_coeffs(iband,irad,i)+re*numer
    end do
    numer = pade_coeffs(iband,irad,0)  +re*numer

    pade_eval_1 = numer/denom
  end function pade_eval_1

end module mo_cloud_optics_rrtmgp_kernels
