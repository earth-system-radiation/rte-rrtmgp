! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2025 - 
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ----------------------------------------------------------------------------
!
! This module is intended to support generic codes using RTE-compatibale gas optics types.
!   It defines variable `gas_optics` of class `ty_gas_optics` and subroutine
!   `load_and_init()` that loads data into the gas_optics variable from a single netCDF file
!    This module might be replaced by others that implement the same variable and procedure
!
! Gas optics classes need to be initialized with data; for RRTMGP data comes from a netCDF file.
!    The gas optics classes themselves don't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward serial implementation of reading the data
!    and calling gas_optics%load().
!
!
module mo_gas_optics_defs
  use mo_rte_kind,           only: wp, wl
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_testing_utils,      only: stop_on_err
  ! --------------------------------------------------
  use mo_simple_netcdf, only: read_field, read_char_vec, read_logical_vec, var_exists, get_dim_size
  use netcdf
  implicit none

  type(ty_gas_optics_ssm) :: gas_optics
  private
  public :: gas_optics, load_and_init

contains
  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file - or don't bother, with the SSM 
  subroutine load_and_init(kdist, filename, available_gases)
    class(ty_gas_optics_ssm), intent(inout) :: kdist
    character(len=*),         intent(in   ) :: filename
    class(ty_gas_concs),      intent(in   ) :: available_gases ! Which gases does the host model have available?
    ! --------------------------------------------------
  
    call stop_on_err(kdist.configure())
  end subroutine
end module mo_gas_optics_defs