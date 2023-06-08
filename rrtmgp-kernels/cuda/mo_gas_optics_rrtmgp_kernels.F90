module mo_gas_optics_rrtmgp_kernels
    ! CvH: The compiler complains that wp and wl are non-constant in the interface so I propagated
    ! CvH: the ifdef, I would like to fall back to the use statement.

    use mo_rte_kind,       only : wp, wl
    use mo_rte_util_array, only : zero_array

    implicit none
    ! private
    ! public :: interpolation, compute_tau_absorption, compute_tau_rayleigh, compute_Planck_source

!     interface
!         subroutine interpolation(                               &
!                 ncol,nlay,ngas,nflav,neta, npres, ntemp,        &
!                 flavor,                                         &
!                 press_ref_log, temp_ref,press_ref_log_delta,    &
!                 temp_ref_min,temp_ref_delta,press_ref_trop_log, &
!                 vmr_ref,                                        &
!                 play,tlay,col_gas,                              &
!                 jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress) bind(C, name="rrtmgp_interpolation")
! 
!             use, intrinsic :: iso_c_binding, only : c_double, c_float
! 
!             #ifdef RTE_USE_SP
!                 integer, parameter :: wp = c_float
!             #else
!                 integer, parameter :: wp = c_double
!             #endif
! 
!             #ifdef RTE_USE_CBOOL
!                 integer, parameter :: wl = c_bool
!             #else
!                 integer, parameter :: wl = kind(.true.)
!             #endif
! 
!             ! input dimensions
!             integer,                            intent(in) :: ncol,nlay
!               !! physical domain size
!             integer,                            intent(in) :: ngas,nflav,neta,npres,ntemp
!               !! k-distribution table dimensions 
!             integer,     dimension(2,nflav),    intent(in) :: flavor
!               !! index into vmr_ref of major gases for each flavor
!             real(wp),    dimension(npres),      intent(in) :: press_ref_log
!               !! log of pressure dimension in RRTMGP tables 
!             real(wp),    dimension(ntemp),      intent(in) :: temp_ref
!               !! temperature dimension in RRTMGP tables 
!             real(wp),                           intent(in) :: press_ref_log_delta, &
!                                                               temp_ref_min, temp_ref_delta, &
!                                                               press_ref_trop_log
!               !! constants related to RRTMGP tables
!             real(wp),    dimension(2,0:ngas,ntemp), intent(in) :: vmr_ref
!               !! reference volume mixing ratios used in compute "binary species parameter" eta
!         
!             ! inputs from profile or parent function
!             real(wp),    dimension(ncol,nlay),        intent(in) :: play, tlay
!               !! input pressure (Pa?) and temperature (K)
!             real(wp),    dimension(ncol,nlay,0:ngas), intent(in) :: col_gas
!               !! input column gas amount - molecules/cm^2 
!             ! outputs
!             integer,     dimension(ncol,nlay), intent(out) :: jtemp, jpress
!               !! temperature and pressure interpolation indexes 
!             logical(wl), dimension(ncol,nlay), intent(out) :: tropo
!               !! use lower (or upper) atmosphere tables 
!             integer,     dimension(2,    ncol,nlay,nflav), intent(out) :: jeta
!               !! Index for binary species interpolation 
! #if !defined(__INTEL_LLVM_COMPILER) && __INTEL_COMPILER >= 2021
!             ! A performance-hitting workaround for the vectorization problem reported in
!             ! https://github.com/earth-system-radiation/rte-rrtmgp/issues/159
!             ! The known affected compilers are Intel Fortran Compiler Classic
!             ! 2021.4, 2021.5 and 2022.1. We do not limit the workaround to these
!             ! versions because it is not clear when the compiler bug will be fixed, see
!             ! https://community.intel.com/t5/Intel-Fortran-Compiler/Compiler-vectorization-bug/m-p/1362591.
!             ! We, however, limit the workaround to the Classic versions only since the
!             ! problem is not confirmed for the Intel Fortran Compiler oneAPI (a.k.a
!             ! 'ifx'), which does not mean there is none though.
!             real(wp),    dimension(:,       :,   :,    :), intent(out) :: col_mix
! #else
!             real(wp),    dimension(2,    ncol,nlay,nflav), intent(out) :: col_mix
!               !! combination of major species's column amounts (first index is strat/trop)
! #endif
!             real(wp),    dimension(2,2,2,ncol,nlay,nflav), intent(out) :: fmajor
!               !! Interpolation weights in pressure, eta, strat/trop 
!             real(wp),    dimension(2,2,  ncol,nlay,nflav), intent(out) :: fminor
!               !! Interpolation fraction in eta, strat/trop
!         end subroutine interpolation
!     end interface

end module mo_gas_optics_rrtmgp_kernels
 
