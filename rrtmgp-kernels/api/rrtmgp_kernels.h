/* This code is part of RRTMGP

Contacts: Robert Pincus and Eli Mlawer
email:  rrtmgp@aer.com

Copyright 2024-  
   Trustees of Columbia University in the City of New York
   All right reserved.

Use and duplication is permitted under the terms of the
    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause

This header files defines the C bindings for the kernels used in RRTMGP 

*/

extern "C"
{
    void rrtmgp_interpolation(
            int* ncol, int* nlay,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* flavor,
            Float* press_ref_log,
            Float* temp_ref,
            Float* press_ref_log_delta,
            Float* temp_ref_min,
            Float* temp_ref_delta,
            Float* press_ref_trop_log,
            Float* vmr_ref,
            Float* play,
            Float* tlay,
            Float* col_gas,
            int* jtemp,
            Float* fmajor, Float* fminor,
            Float* col_mix,
            Bool* tropo,
            int* jeta,
            int* jpress);

    void rrtmgp_compute_tau_absorption(
            int* ncol, int* nlay, int* nband, int* ngpt,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* nminorlower, int* nminorklower,
            int* nminorupper, int* nminorkupper,
            int* idx_h2o,
            int* gpoint_flavor,
            int* band_lims_gpt,
            Float* kmajor,
            Float* kminor_lower,
            Float* kminor_upper,
            int* minor_limits_gpt_lower,
            int* minor_limits_gpt_upper,
            Bool* minor_scales_with_density_lower,
            Bool* minor_scales_with_density_upper,
            Bool* scale_by_complement_lower,
            Bool* scale_by_complement_upper,
            int* idx_minor_lower,
            int* idx_minor_upper,
            int* idx_minor_scaling_lower,
            int* idx_minor_scaling_upper,
            int* kminor_start_lower,
            int* kminor_start_upper,
            Bool* tropo,
            Float* col_mix, Float* fmajor,
            Float* fminor, Float* play,
            Float* tlay, Float* col_gas,
            int* jeta, int* jtemp,
            int* jpress, Float* tau);

    void rrtmgp_compute_tau_rayleigh(
            int* ncol, int* nlay, int* nband, int* ngpt,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* gpoint_flavor,
            int* band_lims_gpt,
            Float* krayl,
            int* idx_h2o, Float* col_dry, Float* col_gas,
            Float* fminor, int* jeta,
            Bool* tropo, int* jtemp,
            Float* tau_rayleigh);

    void rrtmgp_compute_Planck_source(
            int* ncol, int* nlay, int* nbnd, int* ngpt,
            int* nflav, int* neta, int* npres, int* ntemp,
            int* nPlanckTemp,
            Float* tlay,
            Float* tlev,
            Float* tsfc,
            int* sfc_lay,
            Float* fmajor,
            int* jeta,
            Bool* tropo,
            int* jtemp,
            int* jpress,
            int* gpoint_bands,
            int* band_lims_gpt,
            Float* pfracin,
            Float* temp_ref_min, Float* totplnk_delta,
            Float* totplnk,
            int* gpoint_flavor,
            Float* sfc_src,
            Float* lay_src,
            Float* lev_src,
            Float* sfc_src_jac)
}
