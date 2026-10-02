module mo_gas_optics_utils
  implicit none
  public  :: compute_Planck_source, get_layer_number, interp_tlev_from_tlay
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
    function get_layer_number(ncol, nlay, vmr_h2o, plev) result(col_dry)
      !>
      !> Number density (#/m^-2) of dry air molecules
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
  !--------------------------------------------------------------------------------------------------------------------
  interface
    function interp_tlev_from_tlay(ncol, nlay, tlay, play, plev) result(tlev)
      !>
      !> Temperature at layer boundaries, interpolated from layer centers
      !>
      use mo_rte_kind,      only : wp, wl
      integer,  intent(in) :: ncol, nlay
      real(wp), dimension(ncol, nlay  ), intent(in) :: tlay, play ! Layer temperatures [K], pressures [Pa]
      real(wp), dimension(ncol, nlay+1), intent(in) :: plev       ! Layer boundary pressures [Pa]
      ! output
      real(wp), dimension(ncol, nlay+1) :: tlev ! Level temperatures [K]
    end function interp_tlev_from_tlay
  end interface
end module mo_gas_optics_utils
