module mo_rcemip_profiles
  use mo_rte_kind,           only: wp, wl
  use mo_gas_concentrations, only: ty_gas_concs
  implicit none
  private
  public :: make_rcemip_profile_given_p
!
! Provides temperature, pressure, height, and gas concentrations following
!    the protocol for the Radiative-Convective Equilibrium MIP
!    (Wing et al. 2018, https://doi.org/10.5194/gmd-11-793-2018.)
!
  real(wp), parameter :: molar_mass_h2o = 18.015_wp, molar_mass_air = 28.9644_wp
!
! RCEMIP parameters
!
  real(wp), parameter :: g  = 9.79764, & ! m/s^2
                         Rd  = 287.04, & !/kgK
                         p0  = 101480, & !Pa surface pressure
                         qt  = 1.E-14, & !g/g specific humidity at tropopause. Prior to 02Sept2020, errorenously specified as 10^(-11)
                         zq1 = 4000, & ! m
                         zq2 = 7500, & ! m
                         zt  = 15000, & ! m, tropopause height
                         gamma = 0.0067 ! K/m lapse rate
  real(wp), parameter :: chi_co2 = 348.e-6, chi_ch4 = 1650.e-9, chi_n2o = 306.e-9

  real(wp), parameter :: SST = 295, q0 = 0.01200 ! specific humudity/mass mixing ratio
  real(wp), parameter :: g1 = 3.6478, g2 = 0.83209, g3 = 11.3515 ! ozone profile parameters
contains
  ! ------------------------------------------------------------------------------
  elemental subroutine zt_given_p(p, z, temp, q, o3)
    real(wp), intent(in ) :: p
    real(wp), intent(out) :: z
    real(wp), optional, &
              intent(out) :: temp, q, o3 ! K, specific humidity/mmr, vmr

    real(wp) :: T0, Tv0, Tvt, pt
    real(wp) :: Tv, q_local

    T0 = SST ! surface air temperature is SST

    ! Virtual Temperature at surface and tropopause
    Tv0 = T0*(1 + 0.608*q0) ! virtual temperature at surface
    Tvt = Tv0 - gamma*zt    ! virtual temperature at tropopause

    ! Pressure at tropopause (p(zt))
    pt = p0*(Tvt/Tv0)**(g/(Rd*gamma))

    ! Height, specific humidity, virtual temperature
    if (z .le. zt) then
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


    if(present(gas_concs)) then
      call gas_concs%reset()
      error_msg = gas_concs%set_vmr("co2", chi_co2)
      error_msg = gas_concs%set_vmr("ch4", chi_ch4)
      error_msg = gas_concs%set_vmr("n2o", chi_n2o)
      ! Don't forget water vapor, ozone
    end if
  end function make_rcemip_profile_given_p
  ! ------------------------------------------------------------------------------
end module mo_rcemip_profiles
