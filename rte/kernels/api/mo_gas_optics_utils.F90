module mo_gas_optics_utils
  implicit none
  public  :: compute_Planck_source, get_layer_mass, get_layer_number
  ! ------------------------------------------
  interface compute_Planck_source
    subroutine compute_Planck_source_2D(&
        ncol, nlay, nnu, &
        nus, dnus, T, &
        source) bind(C, name="rte_compute_Planck_source_2D")
      use mo_rte_kind,      only : wp, wl
      integer,  &
        intent(in ) :: ncol, nlay, nnu
      real(wp), dimension(nnu), &
        intent(in ) :: nus, dnus
      real(wp), dimension(ncol, nlay), &
        intent(in ) :: T
      real(wp), dimension(ncol, nlay, nnu), &
        intent(out) :: source
     end subroutine compute_Planck_source_2D

    subroutine compute_Planck_source_1D(&
        ncol, nnu, &
        nus, dnus, T, &
        source) bind(C, name="rte_compute_Planck_source_1D")
      use mo_rte_kind,      only : wp, wl
      integer,  &
        intent(in ) :: ncol, nnu
      real(wp), dimension(nnu), &
        intent(in ) :: nus, dnus
      real(wp), dimension(ncol), &
        intent(in ) :: T
      real(wp), dimension(ncol, nnu), &
        intent(out) :: source
      end subroutine compute_Planck_source_1D
    end interface compute_Planck_source

  !--------------------------------------------------------------------------------------------------------------------
  interface
    subroutine get_layer_mass(ncol, nlay, ngas, vmr, plev, mol_weights, m_dry, layer_mass)
      !>
      !> mass (kg m^-2) each gas in the layer
      !>
      use mo_rte_kind,      only : wp, wl
      integer, intent(in)                                  :: ncol, nlay, ngas
      real(wp), dimension(ngas, ncol, nlay  ), intent(in ) :: vmr
      real(wp), dimension(      ncol, nlay+1), intent(in ) :: plev
      real(wp), dimension(ngas),               intent(in ) :: mol_weights
      real(wp),                                intent(in ) :: m_dry
      real(wp), dimension(ngas, ncol, nlay),   intent(out) :: layer_mass
    end subroutine get_layer_mass
  end interface
  !--------------------------------------------------------------------------------------------------------------------
  interface
    function get_layer_number(ncol, nlay, vmr_h2o, plev) result(col_dry)
      !>
      !> Number density (#/cm^-2) of dry air molecules
      !>    "col_dry" in RRTMGP
      ! input
      use mo_rte_kind,      only : wp, wl
      integer, intent(in) :: ncol, nlay
      real(wp), dimension(ncol, nlay  ), intent(in) :: vmr_h2o  ! volume mixing ratio of water vapor to dry air
      real(wp), dimension(ncol, nlay+1), intent(in) :: plev     ! Layer boundary pressures [Pa]
      ! output
      real(wp), dimension(ncol, nlay) :: col_dry ! Column dry amount
    end function get_layer_number
  end interface
end module mo_gas_optics_utils
