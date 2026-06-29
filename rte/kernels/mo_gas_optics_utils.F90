!> Physical constants, planetary/atmospheric parameters, and utility functions for
!>    low-level gas optics calculations including Planck source functions.
!>
!> layer mass for each species
!> layer number density for each species (TK)
!> The latter two don't have C bindings
!
! Copyright 2026-, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_gas_optics_utils
  use mo_rte_kind,       only : wp, wl
  use mo_gas_optics_constants, &
                         only: boltzmann_k, planck_h, lightspeed, &
                               m_h2o, m_dry, avogad, R_univ_gconst, grav
  use mo_rte_util_array, only : zero_array

  implicit none

  interface compute_Planck_source
    module procedure compute_Planck_source_2D, compute_Planck_source_1D
  end interface

  private :: B_nu
  public  :: compute_Planck_source, get_layer_mass, get_layer_number

contains
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
  subroutine compute_Planck_source_2D(&
      ncol, nlay, nnu, &
      nus, dnus, T, &
      source) bind(C, name="rte_compute_Planck_source_2D")
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

  end subroutine compute_Planck_source_2D
  ! -----------------------------------------
  subroutine compute_Planck_source_1D(&
      ncol, nnu, &
      nus, dnus, T, &
      source) bind(C, name="rte_compute_Planck_source_1D")
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

  end subroutine compute_Planck_source_1D
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
  !--------------------------------------------------------------------------------------------------------------------
  function get_layer_number(ncol, nlay, vmr_h2o, plev) result(col_dry)
    !
    !> Number density (#/cm^-2) of dry air molecules
    !>    "col_dry" in RRTMGP
    ! input
    integer, intent(in)                           :: ncol, nlay
    real(wp), dimension(ncol, nlay  ), intent(in) :: vmr_h2o  ! volume mixing ratio of water vapor to dry air
    real(wp), dimension(ncol, nlay+1), intent(in) :: plev     ! Layer boundary pressures [Pa]
    ! output
    real(wp), dimension(ncol, nlay) :: col_dry ! Column dry amount
    ! ------------------------------------------------
    real(wp):: delta_plev, m_air, fact
    integer :: icol, ilev
    ! ------------------------------------------------
    !$acc                parallel loop gang vector collapse(2) copyin(plev,vmr_h2o)  copyout(col_dry)
    !$omp target teams distribute parallel do simd collapse(2) map(to:plev,vmr_h2o) map(from:col_dry)
    do ilev = 1, nlay
      do icol = 1, ncol
        delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1))
        ! Get average mass of moist air per mole of moist air
        fact = 1._wp / (1.+vmr_h2o(icol,ilev))
        m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact
        col_dry(icol,ilev) = 10._wp * delta_plev * avogad * fact/(1000._wp*m_air*100._wp*grav)
      end do
    end do
  end function get_layer_number
end module mo_gas_optics_utils
