! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2020,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Control over input sanitization in Fortan front-end
!   Module variables can be changed only by calling one of the included subroutine
!
! -------------------------------------------------------------------------------------------------
module mo_rte_config
  use mo_rte_kind, only: wl
  implicit none
  private

  logical(wl), protected, public :: check_extents = .true.
  logical(wl), protected, public :: check_values  = .true.

  interface rte_config_checks
    module procedure rte_config_checks_each, rte_config_checks_all
  end interface
  public :: rte_config_checks
contains
  ! --------------------------------------------------------------
  subroutine rte_config_checks_each(extents, values)
    logical(wl), intent(in) :: extents, values

    check_extents = extents
    check_values  = values
  end subroutine rte_config_checks_each
  ! --------------------------------------------------------------
  subroutine rte_config_checks_all(do_checks)
    logical(wl), intent(in) :: do_checks

    check_extents = do_checks
    check_values  = do_checks
  end subroutine rte_config_checks_all
  ! --------------------------------------------------------------
end module mo_rte_config
