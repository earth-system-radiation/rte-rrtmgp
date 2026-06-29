module mo_gas_optics_constants
  use mo_rte_kind,      only : wp, wl
  use mo_rte_util_array,only : zero_array

  implicit none
  public
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
end module mo_gas_optics_constants
