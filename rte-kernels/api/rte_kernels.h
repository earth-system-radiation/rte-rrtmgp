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

include "rte_types.h"

extern "C"
{
    // SHORTWAVE SOLVERS
    void rte_sw_solver_noscat(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
            Float* tau,
            Float* mu0,
            Float* inc_flux_dir,
            Float* flux_dir);

    void rte_sw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
            Float* tau, Float* ssa, Float* g,
            Float* mu0,
            Float* sfc_alb_dir, Float* sfc_alb_dif,
            Float* inc_flux_dir,
            Float* flux_up, Float* flux_dn, Float* flux_dir,
            Bool* has_dif_bc, Float* inc_flux_dif,
            Bool* do_broadband, Float* flux_up_loc, Float* flux_dn_loc, Float* flux_dir_loc);

    void rte_lw_solver_noscat(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1, int* nmus,
            Float* secants, Float* weights,
            Float* tau, Float* lay_source,
            Float* lev_source,
            Float* sfc_emis, Float* sfc_src,
            Float* inc_flux,
            Float* flux_up, Float* flux_dn,
            Bool* do_broadband, Float* flux_up_loc, Float* flux_dn_loc,
            Bool* do_jacobians, Float* sfc_src_jac, Float* flux_up_jac);

    void rte_lw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
            Float* tau, Float* ssa, Float* g,
            Float* lay_source, Float* lev_source,
            Float* sfc_emis, Float* sfc_src,
            Float* inc_flux,
            Float* flux_up, Float* flux_dn);

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
            Float* tau_in, Float* ssa_in, Float* p_in)
 
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
            int* nbnd, int* band_lims_gpoint)''

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
    void zero_array_1D(int* ni, Float* array);
    void zero_array_2D(int* ni, int* nj, Float* array);
    void zero_array_3D(int* ni, int* nj, int* nk, Float* array);
    void zero_array_4D(int* ni, int* nj, int* nk, int* nl, Float* array)
 
}
