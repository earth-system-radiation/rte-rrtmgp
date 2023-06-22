#include "fluxes_kernels_cuda.h"
#include "tools_gpu.h"


namespace
{
    #include "fluxes_kernels.cu"
}


namespace Fluxes_kernels_cuda
{
    void sum_broadband(
            int ncol, int nlev, int ngpt,
            const Float* gpt_flux, Float* flux)
    {
        const int block_lev = 16;
        const int block_col = 16;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_lev = nlev/block_lev + (nlev%block_lev > 0);

        dim3 grid_gpu(grid_col, grid_lev);
        dim3 block_gpu(block_col, block_lev);

        sum_broadband_kernel<<<grid_gpu, block_gpu>>>(ncol, nlev, ngpt, gpt_flux, flux);
    }


    void net_broadband_precalc(
            int ncol, int nlev,
            const Float* flux_dn, const Float* flux_up,
            Float* flux_net)
    {
        const int block_lev = 16;
        const int block_col = 16;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_lev = nlev/block_lev + (nlev%block_lev > 0);

        dim3 grid_gpu(grid_col, grid_lev);
        dim3 block_gpu(block_col, block_lev);

        net_broadband_precalc_kernel<<<grid_gpu, block_gpu>>>(ncol, nlev, flux_dn, flux_up, flux_net);
    }


    void sum_byband(
            int ncol, int nlev, int ngpt, int nbnd,
            const int* band_lims,
            const Float* gpt_flux,
            Float* bnd_flux)
    {
        const int block_bnd = 1;
        const int block_lev = 16;
        const int block_col = 16;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_lev = nlev/block_lev + (nlev%block_lev > 0);
        const int grid_bnd = nbnd/block_bnd + (nbnd%block_bnd > 0);

        dim3 grid_gpu(grid_col, grid_lev, grid_bnd);
        dim3 block_gpu(block_col, block_lev, grid_bnd);

        sum_byband_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlev, ngpt, nbnd, band_lims,
                gpt_flux, bnd_flux);
    }


    void net_byband_full(
            int ncol, int nlev, int ngpt, int nbnd, const int* band_lims,
            const Float* bnd_flux_dn, const Float* bnd_flux_up, Float* bnd_flux_net)
    {
        const int block_bnd = 1;
        const int block_lev = 16;
        const int block_col = 16;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_lev = nlev/block_lev + (nlev%block_lev > 0);
        const int grid_bnd = nbnd/block_bnd + (nbnd%block_bnd > 0);

        dim3 grid_gpu(grid_col, grid_lev, grid_bnd);
        dim3 block_gpu(block_col, block_lev, grid_bnd);

        net_byband_full_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlev, ngpt, nbnd, band_lims,
                bnd_flux_dn, bnd_flux_up, bnd_flux_net);
    }
}
