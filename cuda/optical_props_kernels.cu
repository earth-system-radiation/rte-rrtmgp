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

#ifndef kernel_tuner
const int loop_unroll_factor_nbnd = 2;
#endif


__global__
void increment_1scalar_by_1scalar_kernel(
            const int ncol, const int nlay, const int ngpt,
            Float* __restrict__ tau1, const Float* __restrict__ tau2)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx = icol + ilay*ncol + igpt*ncol*nlay;
        tau1[idx] = tau1[idx]+tau2[idx];
    }
}


__global__
void increment_2stream_by_2stream_kernel(
            const int ncol, const int nlay, const int ngpt, const Float eps,
            Float* __restrict__ tau1, Float* __restrict__ ssa1, Float* __restrict__ g1,
            const Float* __restrict__ tau2, const Float* __restrict__ ssa2, const Float* __restrict__ g2)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx = icol + ilay*ncol + igpt*ncol*nlay;
        const Float tau1_value = tau1[idx];
        const Float tau2_value = tau2[idx];
        const Float tau12 = tau1_value + tau2_value;
        const Float ssa1_value = ssa1[idx];
        const Float ssa2_value = ssa2[idx];
        const Float tauscat12 = (tau1_value * ssa1_value) + (tau2_value * ssa2_value);

        g1[idx] = ((tau1_value * ssa1_value * g1[idx]) + (tau2_value * ssa2[idx] * g2[idx])) / max(tauscat12, eps);
        ssa1[idx] = tauscat12 / max(eps, tau12);
        tau1[idx] = tau12;
    }
}


__global__
void inc_1scalar_by_1scalar_bybnd_kernel(
            const int ncol, const int nlay, const int ngpt,
            Float* __restrict__ tau1, const Float* __restrict__ tau2,
            const int nbnd, const int* __restrict__ band_lims_gpt)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx_gpt = icol + ilay*ncol + igpt*nlay*ncol;

        #pragma unroll loop_unroll_factor_nbnd
        for (int ibnd=0; ibnd<nbnd; ++ibnd)
        {
            if ( ((igpt+1) >= band_lims_gpt[ibnd*2]) && ((igpt+1) <= band_lims_gpt[ibnd*2+1]) )
            {
                const int idx_bnd = icol + ilay*ncol + ibnd*nlay*ncol;

                tau1[idx_gpt] = tau1[idx_gpt] + tau2[idx_bnd];
            }
        }
    }
}


__global__
void inc_2stream_by_2stream_bybnd_kernel(
            const int ncol, const int nlay, const int ngpt, const Float eps,
            Float* __restrict__ tau1, Float* __restrict__ ssa1, Float* __restrict__ g1,
            const Float* __restrict__ tau2, const Float* __restrict__ ssa2, const Float* __restrict__ g2,
            const int nbnd, const int* __restrict__ band_lims_gpt)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx_gpt = icol + ilay*ncol + igpt*nlay*ncol;

        #pragma unroll loop_unroll_factor_nbnd
        for (int ibnd=0; ibnd<nbnd; ++ibnd)
        {
            if ( ((igpt+1) >= band_lims_gpt[ibnd*2]) && ((igpt+1) <= band_lims_gpt[ibnd*2+1]) )
            {
                const int idx_bnd = icol + ilay*ncol + ibnd*nlay*ncol;
                const Float tau1_value = tau1[idx_gpt];
                const Float tau2_value = tau2[idx_bnd];
                const Float ssa1_value = ssa1[idx_gpt];
                const Float ssa2_value = ssa2[idx_bnd];
                const Float tau12 = tau1_value + tau2_value;
                const Float tauscat12 = (tau1_value * ssa1_value) + (tau2_value * ssa2_value);

                g1[idx_gpt] = ((tau1_value * ssa1_value * g1[idx_gpt]) + (tau2_value * ssa2_value * g2[idx_bnd])) / max(tauscat12, eps);
                ssa1[idx_gpt] = tauscat12 / max(eps, tau12);
                tau1[idx_gpt] = tau12;
            }
        }
    }
}


__global__
void delta_scale_2str_k_kernel(
            const int ncol, const int nlay, const int ngpt, const Float eps,
            Float* __restrict__ tau, Float* __restrict__ ssa, Float* __restrict__ g)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx = icol + ilay*ncol + igpt*nlay*ncol;
        const Float g_value = g[idx];
        const Float ssa_value = ssa[idx];
        const Float f = g_value * g_value;
        const Float wf = ssa_value * f;

        tau[idx] *= (Float(1.) - wf);
        ssa[idx] = (ssa_value - wf) / max(eps,(Float(1.)-wf));
        g[idx] = (g_value - f) / max(eps,(Float(1.)-f));

    }
}
