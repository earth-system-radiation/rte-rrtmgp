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

#ifndef OPTICAL_PROPS_KERNELS_CUDA_H
#define OPTICAL_PROPS_KERNELS_CUDA_H

#include "types.h"


namespace Optical_props_kernels_cuda
{
    void increment_1scalar_by_1scalar(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, const Float* tau_in);

    void increment_2stream_by_2stream(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            const Float* tau_in, const Float* ssa_in, const Float* g_in);

    void inc_1scalar_by_1scalar_bybnd(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, const Float* tau_in,
            int nbnd, const int* band_lims_gpoint);

    void inc_2stream_by_2stream_bybnd(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            const Float* tau_in, const Float* ssa_in, const Float* g_in,
            int nbnd, const int* band_lims_gpoint);

    void delta_scale_2str_k(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout);
}
#endif
