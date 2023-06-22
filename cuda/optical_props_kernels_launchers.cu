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

#include <limits>

#include "types.h"
#include "tuner.h"

namespace
{
    #include "optical_props_kernels.cu"
}

namespace Optical_props_kernels_cuda
{
    void increment_1scalar_by_1scalar(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, const Float* tau_in)

    {
        const int block_gpt = 32;
        const int block_lay = 16;
        const int block_col = 1;

        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);
        const int grid_lay = nlay/block_lay + (nlay%block_lay > 0);
        const int grid_col = ncol/block_col + (ncol%block_col > 0);

        dim3 grid_gpu(grid_col, grid_lay, grid_gpt);
        dim3 block_gpu(block_col, block_lay, block_gpt);

        increment_1scalar_by_1scalar_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlay, ngpt,
                tau_inout, tau_in);
    }


    void increment_2stream_by_2stream(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            const Float* tau_in, const Float* ssa_in, const Float* g_in)
    {
        const int block_gpt = 32;
        const int block_lay = 16;
        const int block_col = 1;

        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);
        const int grid_lay = nlay/block_lay + (nlay%block_lay > 0);
        const int grid_col = ncol/block_col + (ncol%block_col > 0);

        dim3 grid_gpu(grid_col, grid_lay, grid_gpt);
        dim3 block_gpu(block_col, block_lay, block_gpt);

        Float eps = std::numeric_limits<Float>::min() * Float(3.);

        increment_2stream_by_2stream_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlay, ngpt, eps,
                tau_inout, ssa_inout, g_inout,
                tau_in, ssa_in, g_in);
    }


    void inc_1scalar_by_1scalar_bybnd(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, const Float* tau_in,
            int nbnd, const int* band_lims_gpoint)

    {
        const int block_col = 16;
        const int block_lay = 4;
        const int block_gpt = 1;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_lay = nlay/block_lay + (nlay%block_lay > 0);
        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);

        dim3 grid_gpu(grid_col, grid_lay, grid_gpt);
        dim3 block_gpu(block_col, block_lay, block_gpt);

        inc_1scalar_by_1scalar_bybnd_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlay, ngpt,
                tau_inout, tau_in,
                nbnd, band_lims_gpoint);
    }


    void inc_2stream_by_2stream_bybnd(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            const Float* tau_in, const Float* ssa_in, const Float* g_in,
            int nbnd, const int* band_lims_gpoint)
    {
        Tuner_map& tunings = Tuner::get_map();

        dim3 grid(ncol, nlay, ngpt), block;

        Float eps = std::numeric_limits<Float>::min() * Float(3.);

        if (tunings.count("inc_2stream_by_2stream_bybnd_kernel") == 0)
        {
            std::tie(grid, block) = tune_kernel(
                    "inc_2stream_by_2stream_bybnd_kernel",
                    dim3(ncol, nlay, ngpt),
                    {1, 2, 3, 4, 8, 12, 16, 24},
                    {1, 2, 3, 4, 8, 12, 16, 24},
                    {1, 2, 3, 4, 8, 12, 16, 24},
                    inc_2stream_by_2stream_bybnd_kernel,
                    ncol, nlay, ngpt, eps,
                    tau_inout, ssa_inout, g_inout,
                    tau_in, ssa_in, g_in,
                    nbnd, band_lims_gpoint);

            tunings["inc_2stream_by_2stream_bybnd_kernel"].first = grid;
            tunings["inc_2stream_by_2stream_bybnd_kernel"].second = block;
        }
        else
        {
            grid = tunings["inc_2stream_by_2stream_bybnd_kernel"].first;
            block = tunings["inc_2stream_by_2stream_bybnd_kernel"].second;
        }

        inc_2stream_by_2stream_bybnd_kernel<<<grid, block>>>(
                ncol, nlay, ngpt, eps,
                tau_inout, ssa_inout, g_inout,
                tau_in, ssa_in, g_in,
                nbnd, band_lims_gpoint);
    }


    void delta_scale_2str_k(
            int ncol, int nlay, int ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout)
    {
        const int block_gpt = 32;
        const int block_lay = 16;
        const int block_col = 1;

        const int grid_gpt  = ngpt/block_gpt + (ngpt%block_gpt > 0);
        const int grid_lay  = nlay/block_lay + (nlay%block_lay > 0);
        const int grid_col  = ncol/block_col + (ncol%block_col > 0);

        dim3 grid_gpu(grid_col, grid_lay, grid_gpt);
        dim3 block_gpu(block_col, block_lay, block_gpt);

        Float eps = std::numeric_limits<Float>::min()*Float(3.);

        delta_scale_2str_k_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlay, ngpt, eps,
                tau_inout, ssa_inout, g_inout);
    }
}
