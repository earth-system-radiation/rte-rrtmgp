/* This code is part of RRTMGP

Contacts: Robert Pincus and Eli Mlawer
email:  rrtmgp@aer.com

Copyright 2024-  
   Trustees of Columbia University in the City of New York
   All right reserved.

Use and duplication is permitted under the terms of the
    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause

This header files defines the C bindings for the kernels used in RRTMGP 
  Adapted from code written by Chiel van Heerwaarden at Wageningen University and Research 

*/

include "rte_types.h"

extern "C"
{
    void rrtmgp_interpolation(
            const int ncol&, const int nlay&,
            const int ngas&, const int nflav&, const int neta&, 
            const int npres&, const int ntemp&,
            int* flavor,
              // [in] (2,nflav)
            Float* press_ref_log,
              // [in] (npres)
            Float* temp_ref,
              // [in] (ntemp)
            const Float press_ref_log_delta&,
            const Float temp_ref_min&,
            const Float temp_ref_delta&,
            const Float press_ref_trop_log&,
            Float* vmr_ref,
              // [in] (2,ngas+1,ntemp)
            Float* play,
              // [in] (ncol,nlay) 
            Float* tlay,
              // [in] (ncol,nlay) 
            Float* col_gas,
              // [in] (ncol,nlay,ngas+1)
            int* jtemp,
              // [out] (ncol*nlay) 
            Float* fmajor,
              // [out] (2,2,2,ncol,nlay,nflav) 
            Float* fminor,
              // [out[ (2,2,  ncol,nlay,nflav))
            Float* col_mix,
              // [out] (2,    ncol,nlay,nflav)
            Bool* tropo,
              // [out] size (ncol*nlay) 
            int* jeta,
              // [out] size (2*ncol*nlay*nflav) 
            int* jpress);
              // [out] size (ncol*nlay) 

    void rrtmgp_compute_tau_absorption(
            const int ncol&, const int nlay&, const int nband&, const int ngpt&,
            const int ngas&, const int nflav&, const int neta&, 
            const int npres&, const int ntemp&,
            const int nminorlower&, const int nminorklower&,
            const int nminorupper&, const int nminorkupper&,
            const int idx_h2o&,
            int* gpoint_flavor,
              // [in] (2,ngpt)
            int* band_lims_gpt,
              // [in] (2,nbnd)
            Float* kmajor,
              // [in] (ntemp,neta,npres+1,ngpt)
            Float* kminor_lower,
              // [in] (ntemp,neta,nminorklower)
            Float* kminor_upper,
              // [in] (ntemp,neta,nminorkupper)
            int* minor_limits_gpt_lower,
              // [in] (2,nminorlower)
            int* minor_limits_gpt_upper,
              // [in] (2,nminorupper)
            Bool* minor_scales_with_density_lower,
              // [in] (  nminorlower)
            Bool* minor_scales_with_density_upper,
              // [in] (  nminorupper)
            Bool* scale_by_complement_lower,
              // [in] (  nminorlower)
            Bool* scale_by_complement_upper,
              // [in] (  nminorupper)
            int* idx_minor_lower,
              // [in] (  nminorlower)
            int* idx_minor_upper,
              // [in] (  nminorupper)
            int* idx_minor_scaling_lower,
              // [in] (  nminorlower)
            int* idx_minor_scaling_upper,
              // [in] (  nminorupper)
            int* kminor_start_lower,
              // [in] (  nminorlower)
            int* kminor_start_upper,
              // [in] (  nminorupper)
            Bool* tropo,
              // [in] (ncol,nlay)
            Float* col_mix, 
              // [in] (2,    ncol,nlay,nflav       )
            Float* fmajor,
              // [in] (2,2,2,ncol,nlay,nflav       )
            Float* fminor, 
             // [in] (2,2,  ncol,nlay,nflav       )
            Float* play,
              // [in] (ncol,nlay)
            Float* tlay, 
              // [in] (ncol,nlay)
            Float* col_gas,
              // [in] (ncol,nlay,ngas+1)
            int* jeta, 
              // [in] (2,    ncol,nlay,nflav       )
            int* jtemp,
              // [in] (ncol,nlay)
            int* jpress, 
              // [in] (ncol,nlay)
            Float* tau);
              // [inout] (ncol,nlay.ngpt)

    void rrtmgp_compute_tau_rayleigh(
            const int ncol&, const int nlay&, const int nband&, const int ngpt&,
            const int ngas&, const int nflav&, const int neta&, const int npres&, const int ntemp&,
            int* gpoint_flavor,
            int* band_lims_gpt,
            Float* krayl,
            int* idx_h2o, Float* col_dry, Float* col_gas,
            Float* fminor, int* jeta,
            Bool* tropo, int* jtemp,
            Float* tau_rayleigh);

    void rrtmgp_compute_Planck_source(
            const int ncol&, const int nlay&, int* nbnd, const int ngpt&,
            const int nflav&, const int neta&, const int npres&, const int ntemp&,
            int nPlanckTemp&,
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
