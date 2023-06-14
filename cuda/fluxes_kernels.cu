/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/earth-system-radiation/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2020,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/earth-system-radiation/rte-rrtmgp-cpp
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


__global__
void sum_broadband_kernel(
            const int ncol, const int nlev, const int ngpt,
            const Float* __restrict__ spectral_flux, Float* __restrict__ broadband_flux)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilev = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (ilev < nlev) )
    {
        const int idx_out = icol + ilev*ncol;
        Float bb_flux_s = 0;
        for (int igpt=0; igpt < ngpt; ++igpt)
        {
            const int idx_in = icol + ilev*ncol + igpt*nlev*ncol;
            bb_flux_s += spectral_flux[idx_in];
        }
        broadband_flux[idx_out] = bb_flux_s;
    }
}


__global__
void net_broadband_precalc_kernel(
            const int ncol, const int nlev,
            const Float* __restrict__ flux_dn, const Float* __restrict__ flux_up,
            Float* __restrict__ broadband_flux_net)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilev = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (ilev < nlev) )
    {
        const int idx = icol + ilev*ncol;
        broadband_flux_net[idx] = flux_dn[idx] - flux_up[idx];
    }
}


__global__
void sum_byband_kernel(
            const int ncol, const int nlev, const int ngpt, const int nbnd,
            const int* __restrict__ band_lims, const Float* __restrict__ spectral_flux,
            Float* __restrict__ byband_flux)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilev = blockIdx.y*blockDim.y + threadIdx.y;
    const int ibnd = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilev < nlev) && (ibnd < nbnd) )
    {
        const int idx_bnd = icol + ilev*ncol + ibnd*ncol*nlev;
        const int gpt_start = band_lims[2*ibnd];
        const int gpt_end = band_lims[2*ibnd+1];

        byband_flux[idx_bnd] = 0;

        for (int igpt = gpt_start; igpt < gpt_end; ++igpt)
        {
            const int idx_gpt = icol + ilev*ncol + igpt*ncol*nlev;
            byband_flux[idx_bnd] += spectral_flux[idx_gpt];
        }
    }
}


__global__
void net_byband_full_kernel(
            const int ncol, const int nlev, const int ngpt, const int nbnd,
            const int* __restrict__ band_lims, const Float* __restrict__ spectral_flux_dn,
            const Float* __restrict__ spectral_flux_up, Float* __restrict__ byband_flux_net)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilev = blockIdx.y*blockDim.y + threadIdx.y;
    const int ibnd = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilev < nlev) && (ibnd < nbnd) )
    {
        const int idx_bnd = icol + ilev*ncol + ibnd*ncol*nlev;
        const int gpt_start = band_lims[2*ibnd];
        const int gpt_end = band_lims[2*ibnd+1];
        byband_flux_net[idx_bnd] = 0;
        for (int igpt = gpt_start; igpt < gpt_end; ++igpt)
        {
            const int idx_gpt = icol + ilev*ncol + igpt*ncol*nlev;
            byband_flux_net[idx_bnd] += spectral_flux_dn[idx_gpt] - spectral_flux_up[idx_gpt];
        }
    }
}
