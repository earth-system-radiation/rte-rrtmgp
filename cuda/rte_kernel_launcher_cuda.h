/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/RobertPincus/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2020,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/microhh/rte-rrtmgp-cpp
 *
 * Contact: Chiel van Heerwaarden
 * email: chiel.vanheerwaarden@wur.nl
 *
 * Copyright 2020, Wageningen University & Research.
 *
 * Use and duplication is permitted under the terms of the
 * BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
 *
 */

#ifndef RTE_KERNELS_CUDA_H
#define RTE_KERNELS_CUDA_H

#include "Array.h"
#include "Types.h"


namespace rte_kernel_launcher_cuda
{
    void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
                  const Float* inc_flux_dir, const Float* mu0, Float* gpt_flux_dir);
    void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, Float* gpt_flux_dn);
    void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const Float* inc_flux_dif, Float* gpt_flux_dn);


    void sw_solver_2stream(
            const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
            const Float* tau, const Float* ssa, const Float* g,
            const Float* mu0,
            const Float* sfc_alb_dir, const Float* sfc_alb_dif,
            const Float* inc_flux_dir,
            Float* flux_up, Float* flux_dn, Float* flux_dir,
            const Bool has_dif_bc, const Float* inc_flux_dif,
            const Bool do_broadband, Float* flux_up_loc, Float* flux_dn_loc, Float* flux_dir_loc);


    void lw_solver_noscat_gaussquad(
            const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const int nmus,
            const Float* secants, const Float* weights,
            const Float* tau, const Float* lay_source,
            const Float* lev_source_inc, const Float* lev_source_dec,
            const Float* sfc_emis, const Float* sfc_src,
            const Float* inc_flux,
            Float* flux_up, Float* flux_dn,
            const Bool do_broadband, Float* flux_up_loc, Float* flux_dn_loc,
            const Bool do_jacobians, const Float* sfc_src_jac, Float* flux_up_jac);


    void lw_secants_array(
            const int ncol, const int ngpt, const int n_quad_angs, const int max_gauss_pts,
            const Float* Gauss_Ds, Float* secants);
}
#endif
