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

#ifndef FLUXES_KERNELS_CUDA_H
#define FLUXES_KERNELS_CUDA_H

#include "types.h"


namespace Fluxes_kernels_cuda
{
    void sum_broadband(
            int ncol, int nlev, int ngpt,
            const Float* gpt_flux, Float* flux);

    void net_broadband_precalc(
            int ncol, int nlev,
            const Float* broadband_flux_dn, const Float* broadband_flux_up,
            Float* broadband_flux_net);

    void sum_byband(
            int ncol, int nlev, int ngpt, int nbnd,
            const int* band_lims,
            const Float* gpt_flux,
            Float* bnd_flux);

    void net_byband_full(
            int ncol, int nlev, int ngpt, int nbnd, const int* band_lims,
            const Float* bnd_flux_dn, const Float* bnd_flux_up, Float* bnd_flux_net);
}
#endif
