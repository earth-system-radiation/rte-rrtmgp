! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2022-  Atmospheric and Environmental Research and
! Trustees of Columbia University in New York.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ------------------------------------------------------------------------------
!
! Provides temperature, pressure, height, and gas concentrations following
!    the protocol for the Radiative-Convective Equilibrium MIP
!    (Wing et al. 2018, https://doi.org/10.5194/gmd-11-793-2018.)
!
! ------------------------------------------------------------------------------
module mo_rcemip_profiles
  use mo_rte_kind,           only: wp, wl
  use mo_rte_util_array_validation, &  
                             only: extents_are
  use mo_gas_concentrations, only: ty_gas_concs
  implicit none
  private
  public :: make_rcemip_profile_given_p, make_rcemip_profiles

  character(len=3), dimension(6), parameter :: rcemip_gases = ["h2o", "co2", "o3 ", "ch4", "n2o", "o2 "]
  real(wp), parameter :: molar_mass_h2o = 18.015_wp, molar_mass_air = 28.9644_wp
  !
  ! RCEMIP parameters
  !
  real(wp), parameter :: g  = 9.79764, & ! m/s^2
                         Rd  = 287.04, & !/kgK
                         p0  = 101480, & !Pa surface pressure
                         qt  = 1.E-14, & !g/g specific humidity at tropopause. Prior to 02Sept2020, errorenously specified as 10^(-11)
                         zq1 = 4000,   & ! m
                         zq2 = 7500,   & ! m
                         zt  = 15000,  & ! m, tropopause height
                         gamma = 0.0067  ! K/m lapse rate
  real(wp), parameter :: chi_co2 = 348.e-6, chi_ch4 = 1650.e-9, chi_n2o = 306.e-9

  real(wp), parameter :: SST = 295, q0 = 0.01200 ! specific humudity/mass mixing ratio
  real(wp), parameter :: g1 = 3.6478, g2 = 0.83209, g3 = 11.3515 ! ozone profile parameters
  ! For fixed SST
  real(wp), parameter :: T0 = SST                ! surface air temperature is SST
  real(wp), parameter :: Tv0 = T0*(1 + 0.608*q0) ! virtual temperature at surface
  real(wp), parameter :: Tvt = Tv0 - gamma*zt    ! virtual temperature at tropopause
  ! Pressure at tropopause (p(zt))
  real(wp), parameter :: pt = p0*(Tvt/Tv0)**(g/(Rd*gamma))

contains
  ! ------------------------------------------------------------------------------
  !
  ! Given pressure in Pa, provide values of temperature in K, height in m,
  !   and (optionally) water vapor and ozone volume mixing ratios
  !
  elemental subroutine zt_given_p(p, z, temp, q, o3)
    real(wp), intent(in ) :: p
    real(wp), intent(out) :: z
    real(wp), optional, &
              intent(out) :: temp, q, o3 ! K, specific humidity/mmr, vmr
    ! ----------------------------------
    real(wp) :: Tv, q_local
    ! ----------------------------------
    ! Height, specific humidity, virtual temperature
    if (p > pt) then
      z = (Tv0/gamma)*(1-(p/p0)**((Rd*gamma)/g))
      q_local = q0*exp(-z/zq1)*exp(-(z/zq2)**2)
      Tv = Tv0 - gamma*z
    else
      z  = zt + (Rd*Tvt/g)*log(pt/p)
      q_local  = qt
      Tv = Tvt
    end if

    ! Temperature
    if(present(temp)) temp = Tv/(1 + 0.608*q_local)
    ! specific humidity/mass mixing ratio -> molar mixing ratio
    if(present(q   )) q  = q_local * molar_mass_air/molar_mass_h2o
    ! Eq 1, converting from Pa to hPa and from ppmv to molar mixing ratio
    if(present(o3  )) o3 = g1 * (p/100._wp)**g2 * exp(-p/(100._wp * g3)) * 1.e-6
  end subroutine zt_given_p
  ! ------------------------------------------------------------------------------
  function make_rcemip_profile_given_p(p, temp, z, gas_concs) result(error_msg)
    real(wp), dimension(:), intent(in   ) :: p
    real(wp), dimension(:), intent(  out) :: temp, z
    type(ty_gas_concs),     optional, &
                            intent(inout) :: gas_concs
    character(len=128)                    :: error_msg
    ! ----------------------------------
    integer :: nlay
    real(wp), dimension(size(p)) :: q, o3
    ! ----------------------------------
    error_msg = ""
    nlay = size(p)
    if(.not. extents_are(temp, nlay)) &
      error_msg = "make_rcemip_profile: p, temp have different lengths"
    if(.not. extents_are(z,    nlay)) &
      error_msg = "make_rcemip_profile: p, z    have different lengths"
    if(len_trim(error_msg) /= 0) return

    call zt_given_p(p, z, temp, q, o3)
    if(present(gas_concs)) then
      error_msg = gas_concs%init(rcemip_gases)
      ! Scalars
      error_msg = gas_concs%set_vmr("co2", chi_co2)
      error_msg = gas_concs%set_vmr("ch4", chi_ch4)
      error_msg = gas_concs%set_vmr("n2o", chi_n2o)
      ! Profiles
      error_msg = gas_concs%set_vmr("h2o", q)
      error_msg = gas_concs%set_vmr("o3" , o3)
      error_msg = gas_concs%set_vmr("o2" , .21_wp)
    end if
    print *, "tropopause p", p0*(Tvt/Tv0)**(g/(Rd*gamma))
  end function make_rcemip_profile_given_p
  ! ------------------------------------------------------------------------------
  function make_rcemip_profiles(p_lay, p_lev, t_lay, z_lay, gas_concs, p_min) result(error_msg)
    real(wp), dimension(:), intent(  out) :: p_lay, p_lev, t_lay, z_lay
    type(ty_gas_concs),     intent(inout) :: gas_concs
    real(wp),     optional, intent(in   )    :: p_min
    character(len=128)                    :: error_msg
    ! ----------------------------------
    real(wp) :: p_top, t_lev(size(p_lev))
    integer  :: i, nlay
    ! ----------------------------------
    error_msg = ""
    nlay = size(p_lay)
    if(.not. extents_are(t_lay, nlay)) &
      error_msg = "make_rcemip_profile: p_lay, t_lay have different lengths"
    if(.not. extents_are(z_lay, nlay)) &
      error_msg = "make_rcemip_profile: p_lay, z_lay have different lengths"
    if(.not. extents_are(p_lev, nlay+1)) &
      error_msg = "make_rcemip_profile: p_lev should be 1 element longer than p_lay"
    if(len_trim(error_msg) /= 0) return

    p_top = 1. ! Min RRTMGP pressure at this time is 1 Pa
    if(present(p_min)) p_top = p_min
    ! Equal spacing in pressure
    p_lev = [(p_top + (p0 - p_top)/nlay * i, i = 0, nlay)]
    ! Layer center is mean pressure
    p_lay = 0.5_wp * (p_lev(1:nlay) + p_lev(2:))

    error_msg = make_rcemip_profile_given_p(p_lay, t_lay, z_lay, gas_concs)
  end function make_rcemip_profiles
  ! ------------------------------------------------------------------------------
end module mo_rcemip_profiles
