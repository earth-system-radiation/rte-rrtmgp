module mo_gas_optics_ddq_kernels
  use mo_rte_kind,             only: wp, wl
  use mo_gas_optics_constants, only: boltzmann_k, lightspeed, planck_h

  implicit none
  private
  public :: tau_absorption_from_fits, add_tau_rayleigh
  integer, parameter, public :: fax_norder = 2, fax_nterms = 2, xsec_nterms = 3

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute absorption optical depth from second-order polynomial approximations
  !    for absorption cross-section
  !
  subroutine tau_absorption_from_fits(ncol, nlay, nnu, ngas, &
                  nus, &
                  play, tlay, dry_num, vmrs, &
                  fax_ngas, fax_num_index, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, fax_vmr0, &
                  xsec_ngas, xsec_num_index, xsec_p, &
                  mtckd_ngas, mtckd_num_index, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0, &
                  tau) bind(C, name="ddq_compute_tau_absorption")
    integer,  intent(in)  :: ncol, nlay, nnu, ngas
    real(wp), intent(in)  :: nus(nnu)
    real(wp), dimension(     ncol, nlay), &
              intent(in)  :: play, tlay, dry_num
! VMRs start at ngas = 0, with element 0 having value 0, to stand in for gases users haven't supplied.
    real(wp), intent(in)  :: vmrs(ncol, nlay, 0:ngas)

    ! Functional approximations to cross-sections
    integer,  intent(in) :: fax_ngas
    integer,  intent(in) :: fax_num_index(fax_ngas)
    real(wp), dimension(0:2, fax_ngas, nnu), &
              intent(in)  :: fax_a, fax_b
    real(wp), intent(in)  :: fax_c(0:2, fax_ngas, nnu)
    real(wp), intent(in)  :: fax_sigma0(fax_ngas, nnu)
    real(wp), dimension(fax_ngas), &
              intent(in)  :: fax_S, fax_T0, fax_p0, fax_vmr0
    ! Cross-section fits
    integer,  intent(in) :: xsec_ngas
    integer,  intent(in) :: xsec_num_index(xsec_ngas)
    real(wp), intent(in) :: xsec_p(0:3, xsec_ngas, nnu)

    integer,  intent(in) :: mtckd_ngas
    integer,  intent(in) :: mtckd_num_index(mtckd_ngas)
    real(wp), dimension(mtckd_ngas, nnu), &
              intent(in) :: mtckd_cself, mtckd_cfrgn, mtckd_n
    real(wp), intent(in) :: mtckd_T0, mtckd_p0

    real(wp), intent(out) :: tau(ncol, nlay, nnu)
    ! -----------------
    integer  :: igas, icol, ilay, inu
    real(wp) :: vmr
    real(wp) :: x, w, q, num, den, delta_T, c0, lin, rs, acc
    real(wp) :: cself, cfrgn, R ! MT_CKD
    real(wp) :: pr, tr, wn      ! MT_CKD per-layer scalars
    real(wp) :: fax_x  (ncol, nlay, fax_ngas), fax_dT(ncol, nlay, fax_ngas), fax_w(ncol, nlay, fax_ngas)
    real(wp) :: xsec_w (ncol, nlay, xsec_ngas)
    real(wp) :: mt_logr(ncol, nlay, mtckd_ngas), &
                mt_sw  (ncol, nlay, mtckd_ngas), mt_fw  (ncol, nlay, mtckd_ngas)
    real(wp) :: inv_2kT(ncol, nlay)
    real(wp) :: fax_rs (fax_ngas, nnu)

    !$acc        data copyin(nus, play, tlay, dry_num, vmrs)                                                 &
    !$acc             copyin(fax_num_index, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, fax_vmr0) &
    !$acc             copyin(xsec_num_index, xsec_p)                                                         &
    !$acc             copyin(mtckd_num_index, mtckd_cself, mtckd_cfrgn, mtckd_n)                             &
    !$acc             copyout(tau)                                                                           &
    !$acc             create(fax_x, fax_dT, fax_w, fax_rs, xsec_w, mt_logr, mt_sw, mt_fw, inv_2kT)
    !$omp target data map(to:nus, play, tlay, dry_num, vmrs)                                                 &
    !$omp             map(to:fax_num_index, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, fax_vmr0) &
    !$omp             map(to:xsec_num_index, xsec_p)                                                         &
    !$omp             map(to:mtckd_num_index, mtckd_cself, mtckd_cfrgn, mtckd_n)                             &
    !$omp             map(from:tau)                                                                          &
    !$omp             map(alloc:fax_x, fax_dT, fax_w, fax_rs, xsec_w, mt_logr, mt_sw, mt_fw, inv_2kT)
    !
    ! Hoist all nu-independent quantities out of the spectral loop
    !
    !$acc                         parallel loop    collapse(3) private(vmr)
    !$omp target teams distribute parallel do simd collapse(3) private(vmr)
    do igas = 1, fax_ngas
      do ilay = 1, nlay
        do icol = 1, ncol
          vmr = vmrs(icol, ilay, fax_num_index(igas))
          ! Increase pressure to account for self-broadening
          fax_x (icol, ilay, igas) = play(icol, ilay) / fax_p0(igas)        &
                                   * (1._wp + vmr            * fax_S(igas)) &
                                   / (1._wp + fax_vmr0(igas) * fax_S(igas))
          fax_dT(icol, ilay, igas) = tlay(icol, ilay) - fax_T0(igas)
          fax_w (icol, ilay, igas) = vmr * dry_num(icol, ilay)   ! Integrated number density [mol/m**2]
        end do
      end do
    end do
    !$acc                         parallel loop    collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do inu = 1, nnu
      do igas = 1, fax_ngas
        fax_rs(igas, inu) = 1._wp / (1._wp + fax_c(1, igas, inu))
      end do
    end do
    !$acc                         parallel loop    collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
    do igas = 1, xsec_ngas
      do ilay = 1, nlay
        do icol = 1, ncol
          xsec_w(icol, ilay, igas) = vmrs(icol, ilay, xsec_num_index(igas)) * dry_num(icol, ilay)
        end do
      end do
    end do
    !$acc                         parallel loop    collapse(3) private(vmr, pr, tr, wn)
    !$omp target teams distribute parallel do simd collapse(3) private(vmr, pr, tr, wn)
    do igas = 1, mtckd_ngas
      do ilay = 1, nlay
        do icol = 1, ncol
          vmr = vmrs(icol, ilay, mtckd_num_index(igas))
          pr  = play(icol, ilay)/mtckd_p0
          tr  = mtckd_T0/tlay(icol, ilay)
          wn  = vmr * dry_num(icol, ilay)
          mt_logr(icol, ilay, igas) = log(tr)
          mt_sw  (icol, ilay, igas) = pr *          vmr  * wn
          mt_fw  (icol, ilay, igas) = pr * (1._wp - vmr) * wn * tr
        end do
      end do
    end do
    !$acc                         parallel loop    collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do ilay = 1, nlay
      do icol = 1, ncol
        inv_2kT(icol, ilay) = (planck_h * lightspeed * 100._wp) &
                            / (2._wp * boltzmann_k * tlay(icol, ilay))
      end do
    end do

    !$acc                         parallel loop    collapse(2) &
    !$acc   private(acc, x, delta_T, c0, lin, rs, w, q, num, den, cself, cfrgn, R)
    !$omp target teams distribute parallel do simd collapse(2) &
    !$omp   private(acc, x, delta_T, c0, lin, rs, w, q, num, den, cself, cfrgn, R)
    do ilay = 1, nlay
      do icol = 1, ncol
        !$acc loop seq
        do inu = 1, nnu
          acc = 0._wp
          !
          ! Functional approximation to cross-sections
          !
          !$acc loop seq
          do igas = 1, fax_ngas
            x       = fax_x (icol, ilay, igas)
            delta_T = fax_dT(icol, ilay, igas)
            c0      = fax_c(0, igas, inu)
            lin     = fax_c(2, igas, inu)
            rs      = fax_rs(igas, inu)
            w       = rs * x + (1._wp - rs)
            q       = c0 + (1._wp - c0) * w * w
            num     = (lin * x * q + (1._wp - lin) * w) &
                    * (fax_a(0, igas, inu) + delta_T * (fax_a(1, igas, inu) + fax_a(2, igas, inu) * delta_T))
            den     = q &
                    * (fax_b(0, igas, inu) + delta_T * (fax_b(1, igas, inu) + fax_b(2, igas, inu) * delta_T))
            acc = acc &
              + fax_sigma0(igas, inu) * max(num / den, 0._wp) & ! cross-section [m**2/mol]
              * fax_w(icol, ilay, igas)
          end do
          !
          ! Cross-sections pressure and temperature dependence following doi:10.1029/2022MS003239
          !
          !$acc loop seq
          do igas = 1, xsec_ngas
            acc = acc &
              + (xsec_p(0, igas, inu)                                                         &
                 + tlay(icol, ilay) * (xsec_p(1, igas, inu) + xsec_p(2, igas, inu) * tlay(icol, ilay)) &
                 + xsec_p(3, igas, inu) * play(icol, ilay))                                   &
              * xsec_w(icol, ilay, igas)
          end do
          !
          ! MT_CKD continuum
          !
          !$acc loop seq
          do igas = 1, mtckd_ngas
            ! (T0/T)**(1+n) written as exp((1+n)*log(T0/T)) with the log hoisted
            cself = exp((1._wp + mtckd_n(igas, inu)) * mt_logr(icol, ilay, igas)) &
                  * mt_sw(icol, ilay, igas) * mtckd_cself(igas, inu)
            cfrgn = mt_fw(icol, ilay, igas) * mtckd_cfrgn(igas, inu)
            ! nu supplied in kaysers (cm^-1); convert to MKS for tanh argument;
            ! R needs to be in units cm^-1 as cself and cfrgn are in of m^2/molecule cm
            ! R * (cself + cfrgn) is in units of m^2/molecule
            R = nus(inu) * tanh(nus(inu) * inv_2kT(icol, ilay))
            acc = acc + R * (cself + cfrgn)
          end do
          tau(icol, ilay, inu) = max(0._wp, acc)
        end do
      end do
    end do

    !$acc end data
    !$omp end target data
  end subroutine tau_absorption_from_fits
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Add Rayleigh scattering to absorption optical depth; compute single-scattering albedo
  !
  subroutine add_tau_rayleigh(ncol, nlay, nnu,  &
                              dry_num,          &
                              rayleigh_xsec,    &
                              tau, ssa) bind(C, name="ddq_add_tau_rayleigh")
    integer,  intent(in)    :: ncol, nlay, nnu
    real(wp), intent(in)    :: dry_num(ncol, nlay)
    real(wp), intent(in)    :: rayleigh_xsec(nnu)
    real(wp), intent(inout) :: tau(ncol, nlay, nnu), ssa(ncol, nlay, nnu)
    ! -----------------
    integer  :: icol, ilay, inu
    real(wp) :: t, t_r

    !$acc                         parallel loop    collapse(3) private(t, t_r) &
    !$acc   copyin(dry_num, rayleigh_xsec) copy(tau) copyout(ssa)
    !$omp target teams distribute parallel do simd collapse(3) private(t, t_r) &
    !$omp   map(to:dry_num, rayleigh_xsec) map(tofrom:tau) map(from:ssa)
    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          t_r = dry_num(icol, ilay) * rayleigh_xsec(inu)
          t = tau(icol, ilay, inu)
          tau(icol, ilay, inu) = t + t_r
          ssa(icol, ilay, inu) = t_r/(t + t_r)
        end do
      end do
    end do
  end subroutine add_tau_rayleigh
end module mo_gas_optics_ddq_kernels
