module mo_optics_utils
  use mo_rte_kind,      only: wp, wl
  use mo_optics_ssm,    only: ty_optics_ssm
  use mo_testing_utils, only: stop_on_err
  use mo_gas_concentrations, only: ty_gas_concs
  ! --------------------------------------------------
  implicit none

  type(ty_optics_ssm) :: gas_optics
  private
  public :: gas_optics, load_and_init, configure

contains

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Configures the default values of the SSM
  !
  subroutine configure(kdist)
    class(ty_optics_ssm), intent(inout) :: kdist

    ! call stop_on_err(kdist%configure())
  end subroutine configure
  !--------------------------------------------------------------------------------------------------------------------
  ! Maintains interface compatability with RRTMGP version but just calls configure()
  !
  subroutine load_and_init(kdist, filename, available_gases)
    class(ty_optics_ssm), intent(inout) :: kdist
    character(len=*),     intent(in   ) :: filename
    class(ty_gas_concs),  intent(in   ) :: available_gases ! Which gases does the host model have available?

    call configure(kdist)
  end subroutine load_and_init
  !--------------------------------------------------------------------------------------------------------------------
end module mo_optics_utils
