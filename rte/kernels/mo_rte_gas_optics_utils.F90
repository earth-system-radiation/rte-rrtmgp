!> Physical constants, planetary/atmospheric parameters, and utility functions for 
!>    low-level gas optics calculations including Planck source functions. 
!>
!> layer mass for each species 
!> layer number density for each species (TK)
!
! Copyright 2026-, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_rte_gas_optics_utils
  use mo_rte_kind,      only : wp, wl
  use mo_rte_util_array,only : zero_array

  implicit none

  interface compute_Planck_source
    module procedure compute_Planck_source_2D, compute_Planck_source_1D
  end interface

  private :: B_nu
  public  :: compute_Planck_source

  ! -----------------------------------------
  ! Physical constants, 2018 SI defintion of metric system
  !   doi:10.1088/1681-7575/aa950a (see also https://www.nist.gov/si-redefinition/meet-constants)
  ! Boltzmann therrmodynamics constant [J/K] = [(kg m^2)/(K s^2)]
  real(wp), parameter :: boltzmann_k = 1.380649e-23_wp

  !  molecular weight of water [kg/mol]
  real(wp), parameter :: m_h2o =  0.018016_wp

  ! Avogadro's number [molec/mol]
  real(wp), parameter :: avogad = 6.02214076e23_wp

  ! Universal gas constant [J/(mol K)]
  real(wp), parameter :: R_univ_gconst = avogad * boltzmann_k

  ! Need to verify against NIST values 
  real(wp), parameter :: planck_h     = 6.626075540e-34_wp    ! Planck's constant
  real(wp), parameter :: lightspeed   = 2.99792458e8_wp       ! Speed of light

  ! -----------------------------------------
  !
  ! Constants specific to the earth's atmosphere -- changeable in init() because they
  !   might be different on e.g. other planets

  ! molecular weight of dry air [kg/mol]
  real(wp), protected :: m_dry = 0.028964_wp

  ! Gravity at Earth's surface [m/s2]
  real(wp), protected :: grav = 9.80665_wp

  ! Specific heat at constant pressure for dry air [J/(K kg)]
  real(wp), protected :: cp_dry = 1004.64_wp
! -------------------------------------------------------------------------------------------------
contains
  ! -----------------------------------------
  subroutine init_constants(gravity, mol_weight_dry_air, heat_capacity_dry_air)
    real(wp), optional, intent(in) :: gravity, mol_weight_dry_air, heat_capacity_dry_air
      !! Planetary and atmospheric values used by RRTMGP in computing gas optical properties
      !! Default values reflect modern Earth but these can be changed using this routine

    if(present(gravity))               grav   = gravity
    if(present(mol_weight_dry_air))    m_dry  = mol_weight_dry_air
    if(present(heat_capacity_dry_air)) cp_dry = heat_capacity_dry_air

  end subroutine init_constants
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

   !$acc                         parallel loop    collapse(3)
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

    !$acc                         parallel loop    collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do inu = 1, nnu
      do icol = 1, ncol
        source(icol, inu) = B_nu(T(icol), nus(inu)) * dnus(inu)
      end do
    end do

  end subroutine compute_Planck_source_1D
  ! -------------------------------------------------------------------------------------------------
  ! Layer mass (kg), layer number density 
  ! -------------------------------------------------------------------------------------------------
  subroutine get_layer_mass(ncol, nlay, ngas, vmr, plev, mol_weights, m_dry, layer_mass) &
     bind(C, name="rte_get_layer_mass")
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
    !$acc                         parallel loop    collapse(3)
    !$omp target teams distribute parallel do simd collapse(3)
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
  function get_layer_number(vmr_h2o, plev) result(col_dry)
    !
    !> Number density (#/cm^-2) of dry air molecules
    !>    "col_dry" in RRTMGP
    ! input
    real(wp), dimension(:,:), intent(in) :: vmr_h2o  ! volume mixing ratio of water vapor to dry air; (ncol,nlay)
    real(wp), dimension(:,:), intent(in) :: plev     ! Layer boundary pressures [Pa] (ncol,nlay+1)
    ! output
    real(wp), dimension(size(plev,dim=1),size(plev,dim=2)-1) :: col_dry ! Column dry amount (ncol,nlay)
    ! ------------------------------------------------
    real(wp):: delta_plev, m_air, fact
    integer :: ncol, nlev
    integer :: icol, ilev ! nlay = nlev-1
    ! ------------------------------------------------
    ncol = size(plev, dim=1)
    nlev = size(plev, dim=2)
    !$acc                parallel loop gang vector collapse(2) copyin(plev,vmr_h2o)  copyout(col_dry)
    !$omp target teams distribute parallel do simd collapse(2) map(to:plev,vmr_h2o) map(from:col_dry)
    do ilev = 1, nlev-1
      do icol = 1, ncol
        delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1))
        ! Get average mass of moist air per mole of moist air
        fact = 1._wp / (1.+vmr_h2o(icol,ilev))
        m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact
        col_dry(icol,ilev) = 10._wp * delta_plev * avogad * fact/(1000._wp*m_air*100._wp*grav)
      end do
    end do
  end function get_layer_number
end module mo_rte_gas_optics_utils