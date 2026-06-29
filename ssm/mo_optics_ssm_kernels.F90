! This is an implementation of a simple spectral model
!
! Contacts: Andrew Williams and Robert Pincus
! email:  andrewwilliams@ucsd.edu
!!
! Copyright 2025-,  ... and Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Kernels
!>
!> Details here
! -------------------------------------------------------------------------------------------------
module mo_optics_ssm_kernels
  use mo_rte_kind,      only : wp, wl
  use mo_gas_optics_constants, &
                         only: boltzmann_k, planck_h, lightspeed, &
                               m_h2o, m_dry, avogad, R_univ_gconst, grav
  implicit none
  interface compute_Planck_source_ssm
    module procedure compute_Planck_source_2D_ssm, compute_Planck_source_1D_ssm
  end interface
  private
  public :: compute_tau, compute_Planck_source_ssm, get_layer_mass
contains
  !
  ! Computes optical depth from atmospheric state and gas optics
  !
  subroutine compute_tau(ncol, nlay, nnu, ngas,   &
                         absorption_coeffs, play, pref, layer_mass, &
                         tau) bind(C, name="ssm_compute_tau_absorption")
    integer,  &
      intent(in ) :: ncol, nlay, nnu, ngas
    real(wp), dimension(ngas, nnu), &
      intent(in ) :: absorption_coeffs
    real(wp), dimension(ncol, nlay), &
      intent(in   ) :: play     !! layer pressures [Pa]; (ncol,nlay)
    real(wp), dimension(ngas, ncol, nlay), &
      intent(in ) :: layer_mass !! !! mass of atm layer [kg/m2]; (ngas, ncol, nlay)
    real(wp), &
      intent(in ) :: pref       !! reference pressure [Pa], do foreign broadening if this is non-zero
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: tau

    ! Local variables
    integer :: icol, ilay, inu, igas

    ! Apply pressure broadening if pref input is non-zero
    if (pref /= 0._wp) then
      !$acc                         parallel loop    collapse(3)
      !$omp target teams distribute parallel do simd collapse(3)
      do inu = 1, nnu
        do ilay = 1, nlay
          do icol = 1, ncol
            tau(icol, ilay, inu) = &
              sum( [(layer_mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas)] ) &
              * play(icol, ilay) / pref
          end do
        end do
      end do
    else
      !$acc                         parallel loop    collapse(3)
      !$omp target teams distribute parallel do simd collapse(3)
      do inu = 1, nnu
        do ilay = 1, nlay
          do icol = 1, ncol
            tau(icol, ilay, inu) = &
              sum( [(layer_mass(igas, icol, ilay) * absorption_coeffs(igas, inu), igas = 1, ngas)] )
          end do
        end do
      end do
    end if
  end subroutine compute_tau

  ! -------------------------------------------------------------------------------------------------
  ! Planck source functions
  ! -------------------------------------------------------------------------------------------------
  !
  ! Planck function (gets wrapped by 1D, 2D codes)
  !
  elemental function B_nu(T, nu)
    real(wp), intent(in) :: T, nu
    real(wp)             :: B_nu
    B_nu = 100._wp*2._wp*planck_h*((nu*100._wp)**3)*(lightspeed**2) / &
         (exp((planck_h * lightspeed * nu * 100._wp) / (boltzmann_k * T)) - 1._wp)
  end function
  ! -----------------------------------------
  subroutine compute_Planck_source_2D_ssm(&
      ncol, nlay, nnu, &
      nus, dnus, T, &
      source) bind(C, name="rte_compute_Planck_source_2D_ssm")
    integer,  &
      intent(in ) :: ncol, nlay, nnu
    real(wp), dimension(nnu), &
      intent(in ) :: nus, dnus
    real(wp), dimension(ncol, nlay), &
      intent(in ) :: T
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: source

     ! Local variables
    integer :: icol, ilay, inu

   !$acc                         parallel loop    collapse(3) copyin(nus,dnus,T) copyout(source)
   !$omp target teams distribute parallel do simd collapse(3)
   do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          source(icol, ilay, inu) = B_nu(T(icol, ilay), nus(inu)) * dnus(inu)
        end do
      end do
    end do

  end subroutine compute_Planck_source_2D_ssm
  ! -----------------------------------------
  subroutine compute_Planck_source_1D_ssm(&
      ncol, nnu, &
      nus, dnus, T, &
      source) bind(C, name="rte_compute_Planck_source_1D_ssm")
    integer,  &
      intent(in ) :: ncol, nnu
    real(wp), dimension(nnu), &
      intent(in ) :: nus, dnus
    real(wp), dimension(ncol), &
      intent(in ) :: T
    real(wp), dimension(ncol, nnu), &
      intent(out) :: source

     ! Local variables
     integer :: icol, ilay, inu

    !$acc                         parallel loop    collapse(2) copyin(nus,dnus,T)  copyout(source)
    !$omp target teams distribute parallel do simd collapse(2) map(to:nus,dnus,T) map(from:source)
    do inu = 1, nnu
      do icol = 1, ncol
        source(icol, inu) = B_nu(T(icol), nus(inu)) * dnus(inu)
      end do
    end do

  end subroutine compute_Planck_source_1D_ssm
  ! -------------------------------------------------------------------------------------------------
  ! Layer mass (kg), layer number density
  ! -------------------------------------------------------------------------------------------------
  subroutine get_layer_mass(ncol, nlay, ngas, vmr, plev, mol_weights, m_dry, layer_mass)
    !
    !> mass (kg m^-2) each gas in the layer
    !>
    integer, intent(in)                                  :: ncol, nlay, ngas
    real(wp), dimension(ngas, ncol, nlay  ), intent(in ) :: vmr
    real(wp), dimension(      ncol, nlay+1), intent(in ) :: plev
    real(wp), dimension(ngas),               intent(in ) :: mol_weights
    real(wp),                                intent(in ) :: m_dry
    real(wp), dimension(ngas, ncol, nlay),   intent(out) :: layer_mass

    integer :: igas, icol, ilay
    ! Convert pressures and vmr to layer masses (ngas, ncol, nlay)
    ! mmr = vmr * (Mgas/Mair)
    ! layer_mass = mmr * dp / g
    !$acc                         parallel loop    collapse(3) copyin(vmr,mol_weights,plev)  copyout(layer_mass)
    !$omp target teams distribute parallel do simd collapse(3) map(to:vmr,mol_weights,plev) map(from:layer_mass)
    do ilay = 1, nlay
      do icol = 1, ncol
        do igas = 1, ngas
          layer_mass(igas, icol, ilay) = vmr(igas, icol, ilay) * &
            (mol_weights(igas) / m_dry) * &
            abs(plev(icol, ilay+1) - plev(icol, ilay)) / grav
        end do
      end do
    end do
  end subroutine get_layer_mass

end module mo_optics_ssm_kernels
