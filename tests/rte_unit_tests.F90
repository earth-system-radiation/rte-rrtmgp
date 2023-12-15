! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2023-  Atmospheric and Environmental Research,
!    Regents of the University of Colorado,
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ----------------------------------------------------------------------------
program rte_unit_tests
  !
  ! Exercise various paths through RTE code including solvers, optical properties, fluxes
  !   Tests are run on idealized problems with analytic solutions (e.g. radiative equilibrium)
  !   Some sections test e.g. initialization and finalization 
  !   Others compare two sets of fluxes which should be the same, e.g. with respect to vertical ordering 
  !
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_rte_util_array,     only: zero_array
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_testing_utils,      only: allclose, stop_on_err, report_err
  implicit none
  ! ----------------------------------------------------------------------------------
  !
  ! Do we vary the number of columns? Layers? 
  !   For gray radiative equilibrium: 
  !   We want 
  !     several optical thickness values and several surface temperatures - let's say 8 in total 
  !     even and uneven discretizations in tau?
  !   Want to check correctness
  !     OLR, surface T, and optical depth have the relationships we expect 
  !     net fluxes are constant with height 
  !     net fluxes on- and off-line are the same (this tests ty_fluxes_broadband)
  !     heating rates? couldn't check correctness 
  !   Then use an example to test invariance 
  !     to vertical orientation 
  !     subsetting
  !     maybe not to incrementing since we're restrictied to a single solver 
  !     vertical discretization? Maybe just check boundary fluxes
  !     *Jacobian
  !     *Computing Jacobian shouldn't change net fluxes 
  !   Expand testing to optical properties 1_scl? 
  !     Increment with 3 variants of transparent media
  !     Divide by two and increment 
  !   Test the application of the boundary condition? 
end program rte_unit_tests