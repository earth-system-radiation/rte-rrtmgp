/* This code is part of Radiative Transfer for Energetics (RTE)

Contacts: Robert Pincus and Eli Mlawer
email:  rrtmgp@aer.com

Copyright 2024-  
   Trustees of Columbia University in the City of New York
   All right reserved.

Use and duplication is permitted under the terms of the
    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause

This header files defines the C bindings for the kernels used in RTE 
  Adapted from code written by Chiel van Heerwaarden at Wageningen University and Research 

*/

#include "rte_types.h"

extern "C"
{
    // SHORTWAVE SOLVERS
    void rte_sw_solver_noscat(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
              // [in] ncol, nlay, ngpt: constant scalars
              // [in] top_at_1: constant Bool 
            Float* tau,
              // [in] tau: vector of length ncol*nlay*ngpt 
            Float* mu0,
              // [in] mu0: vector of length ncol*nlay
            Float* inc_flux_dir,
              // [in] inc_flux_dir: vector of length ncol*ngpt 
            Float* flux_dir);
              // [out] flux_dir: vector of length ncol*(nlay+1)*ngpt 

    void rte_sw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
              // [in] ncol, nlay, ngpt: constant scalars
              // [in] top_at_1: constant Bool 
            Float* tau, Float* ssa, Float* g,
              // [in] tau, ssa, g: vectors of length ncol*nlay*ngpt 
            Float* mu0,
              // [in] mu0: vector of length ncol*nlay
            Float* sfc_alb_dir, Float* sfc_alb_dif,
              // [in] sfc_alb_dir, sfc_alb_dif: vector of length ncol*ngpt
            Float* inc_flux_dir,
              // [in] inc_flux_dir: vector of length ncol*ngpt 
            Float* flux_up, Float* flux_dn, Float* flux_dir,
              // [out] flux_dir: vector of length ncol*(nlay+1)*ngpt 
            Bool* has_dif_bc, Float* inc_flux_dif,
              // [in] has_dif_bc: constant Bool 
              // [in] inc_flux_dif: vector of length ncol*ngpt 
            Bool* do_broadband, Float* broadband_up, Float* broadband_dn, Float* broadband_dir);
              // [in] has_dif_bc: constant Bool 
              // [out] broadband_up, broadband_dn, broadband_dir: vector of length ncol*(nlay+1)

    void rte_lw_solver_noscat(
            int* ncol, int* nlay, int* ngpt, 
              // [in] ncol, nlay, ngpt: constant scalars
            Bool* top_at_1, int* nmus,
              // [in] top_at_1: constant Bool 
              // [in] nmus: constant scalar
            Float* secants, Float* weights,
              // [in] secants, weights: vectors of length nmus
            Float* tau, Float* lay_source,
              // [in] taus, lay_source: vectors of length ncol*nlay*ngpt 
            Float* lev_source,
              // [in] lev_source: vectors of length ncol*(nlay+1)*ngpt 
            Float* sfc_emis, Float* sfc_src,
              // [in] sfc_emis, sfc_src: vectors of length ncol*ngpt 
            Float* inc_flux,
              // [in] inc_flux: vector of length ncol*ngpt 
            Float* flux_up, Float* flux_dn,
              // [out] flux_up, flux_dn: vector of length ncol*(nlay+1)*ngpt 
            Bool* do_broadband, 
              // [in] do_broadband: constant Bool 
            Float* broadband_up, Float* broadband_dn,
              // [out] broadband_up, broadband_dn: vector of length ncol*(nlay+1)
            Bool* do_jacobians, 
              // [in] do_jacobians: constant Bool 
            Float* sfc_src_jac, 
              // [in] sfc_src_jac: vector of length ncol*ngpt 
            Float* flux_up_jac,
              // [out] flux_up_jac: vector of length ncol*(nlay+1)*ngpt 
            Bool* do_rescaling,
              // [in] do_rescaling: constant Bool 
            Float* ssa, g);
              // [in] ssa, g: vectors of length ncol*nlay*ngpt 
            

    void rte_lw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
              // [in] ncol, nlay, ngpt: constant scalars
              // [in] top_at_1: constant Bool 
            Float* tau, Float* ssa, Float* g,
               // [in] tau, ssa, g: vectors of length ncol*nlay*ngpt 
            Float* lay_source, Float* lev_source,
               // [in] lay_source: vector of length ncol*nlay*ngpt 
               // [in] lev_source: vector of length ncol*(nlay+1)*ngpt 
            Float* sfc_emis, Float* sfc_src,
              // [in] sfc_emis, sfc_src: vectors of length ncol*ngpt 
            Float* inc_flux,
              // [in] inc_flux: vector of length ncol*ngpt 
            Float* flux_up, Float* flux_dn);
              // [out] flux_up, flux_dn: vector of length ncol*(nlay+1)*ngpt 

    // OPTICAL PROPS - INCREMENT
    void rte_increment_1scalar_by_1scalar(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* tau_in);
 

    void rte_increment_1scalar_by_2stream(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout,
            Float* tau_in, Float* ssa_in);
 
    void rte_increment_1scalar_by_nstream(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout,
            Float* tau_in, Float* ssa_in);
 
    void rte_increment_2stream_by_1scalar(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout,
            Float* tau_in);
 
    void rte_increment_2stream_by_2stream(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            Float* tau_in, Float* ssa_in, Float* g_in);
 

    void rte_increment_2stream_by_nstream(
            int* ncol, int* nlay, int* ngpt, int* nmom,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            Float* tau_in, Float* ssa_in, Float* p_in);
 
    void rte_increment_nstream_by_1scalar(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout,
            Float* tau_in);
 
    void rte_increment_nstream_by_2stream(
            int* ncol, int* nlay, int* ngpt, int* nmom1,
            Float* tau_inout, Float* ssa_inout, Float* p_inout,
            Float* tau_in, Float* ssa_in, Float* g_in);
 
    void rte_increment_nstream_by_nstream(
            int* ncol, int* nlay, int* ngpt, int* nmom1, int* nmom2,
            Float* tau_inout, Float* ssa_inout, Float* p_inout,
            Float* tau_in, Float* ssa_in, Float* p_in);

    // OPTICAL PROPS - INCREMENT BYBND
    void rte_inc_1scalar_by_1scalar_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* tau_in,
            int* nbnd, int* band_lims_gpoint);

    void rte_inc_1scalar_by_2stream_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout,
            Float* tau_in, Float* ssa_in,
            int* nbnd, int* band_lims_gpoint);

    void rte_inc_1scalar_by_nstream_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout,
            int* nbnd, int* band_lims_gpoint);

    void rte_inc_2stream_by_1scalar_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout,
            Float* tau_in,
            int* nbnd, int* band_lims_gpoint);

    void rte_inc_2stream_by_2stream_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            Float* tau_in, Float* ssa_in, Float* g_in,
            int* nbnd, int* band_lims_gpoint);
  
    void rte_inc_2stream_by_nstream_bybnd(
            int* ncol, int* nlay, int* ngpt, int* nmom,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            Float* tau_in, Float* ssa_in, Float* p_in,
            int* nbnd, int* band_lims_gpoint);


    void rte_inc_nstream_by_1scalar_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout,
            Float* tau_in,
            int* nbnd, int* band_lims_gpoint);

    void rte_inc_nstream_by_2stream_bybnd(
            int* ncol, int* nlay, int* ngpt, int* nmom1,
            Float* tau_inout, Float* ssa_inout, Float* p_inout,
            Float* tau_in, Float* ssa_in, Float* g_in,
            int* nbnd, int* band_lims_gpoint);

    void rte_inc_nstream_by_nstream_bybnd(
            int* ncol, int* nlay, int* ngpt, int* nmom1, int* nmom2,
            Float* tau_inout, Float* ssa_inout, Float* p_inout,
            Float* tau_in, Float* ssa_in, Float* p_in,
            int* nbnd, int* band_lims_gpoint);

    // OPTICAL PROPS - DELTA SCALING
    void rte_delta_scale_2str_k(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout);

    void rte_delta_scale_2str_f_k(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout, Float* f);

    // OPTICAL PROPS - SUBSET
    void rte_extract_subset_dim1_3d(
            int* ncol, int* nlay, int* ngpt, Float* array_in, int* ncol_start, int* ncol_end, Float* array_out);

    void rte_extract_subset_dim2_4d(
            int* ncol, int* nlay, int* ngpt, Float* array_in, int* ncol_start, int* ncol_end, Float* array_out);

    void rte_extract_subset_absorption_tau(
            int* ncol, int* nlay, int* ngpt, Float* tau_in, Float* ssa_in, int* ncol_start, int* ncol_end, Float* tau_out);

   // Fluxes - reduction 
    void rte_sum_broadband(
            int* ncol, int* nlev, int* ngpt,
            Float* gpt_flux, Float* flux);

    void rte_net_broadband_full(
            int* ncol, int* nlev, int* ngpt,
            Float* gpt_flux_dn, Float* gpt_flux_up,
            Float* flux_net);

    void rte_net_broadband_precalc(
            int* ncol, int* nlev,
            Float* broadband_flux_dn, Float* broadband_flux_up,
            Float* broadband_flux_net);

    void rte_sum_byband(
            int* ncol, int* nlev, int* ngpt, int* nbnd,
            int* band_lims,
            Float* gpt_flux,
            Float* bnd_flux);

    void rte_net_byband_full(
            int* ncol, int* nlev, int* ngpt, int* nbnd, int* band_lims,
            Float* bnd_flux_dn, Float* bnd_flux_up, Float* bnd_flux_net);

    // Array utilities 
    void zero_array_1D(const int* ni, Float* array);
      // [in]  ni is a constant scalar 
      // [out] array is a vector of size ni
    void zero_array_2D(int* ni, int* nj, Float* array);
      // [in]  ni, nj are constant scalars
      // [out] array is a vector of size ni*nj
    void zero_array_3D(int* ni, int* nj, int* nk, Float* array);
      // [in]  ni, nj, nk are constant scalars
      // [out] array is a vector of size ni*nj*nk
    void zero_array_4D(int* ni, int* nj, int* nk, int* nl, Float* array);
      // [in]  ni, nj, nk, nl are constant scalars
      // [out] array is a vector of size ni*nj*nk*nl
 
}
