#include "fluxes_kernels_cuda.h"
#include "tools_gpu.h"


namespace
{
    #include "fluxes_kernels.cu"

    using Tools_gpu::calc_grid_size;
}


namespace Fluxes_kernels_cuda
{
    void sum_broadband(
            int ncol, int nlev, int ngpt,
            const Float* gpt_flux, Float* flux)
    {
        dim3 block_gpu(16, 16);
        dim3 grid_gpu = calc_grid_size(block_gpu, dim3(ncol, nlev));

        sum_broadband_kernel<<<grid_gpu, block_gpu>>>(ncol, nlev, ngpt, gpt_flux, flux);
    }


    void net_broadband_precalc(
            int ncol, int nlev,
            const Float* flux_dn, const Float* flux_up,
            Float* flux_net)
    {
        dim3 block_gpu(16, 16);
        dim3 grid_gpu = calc_grid_size(block_gpu, dim3(ncol, nlev));

        net_broadband_precalc_kernel<<<grid_gpu, block_gpu>>>(ncol, nlev, flux_dn, flux_up, flux_net);
    }


    void sum_byband(
            int ncol, int nlev, int ngpt, int nbnd,
            const int* band_lims,
            const Float* gpt_flux,
            Float* bnd_flux)
    {
        dim3 block_gpu(16, 16, 1);
        dim3 grid_gpu = calc_grid_size(block_gpu, dim3(ncol, nlev, nbnd));

        sum_byband_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlev, ngpt, nbnd, band_lims,
                gpt_flux, bnd_flux);
    }


    void net_byband_full(
            int ncol, int nlev, int ngpt, int nbnd, const int* band_lims,
            const Float* bnd_flux_dn, const Float* bnd_flux_up, Float* bnd_flux_net)
    {
        dim3 block_gpu(16, 16, 1);
        dim3 grid_gpu = calc_grid_size(block_gpu, dim3(ncol, nlev, nbnd));

        net_byband_full_kernel<<<grid_gpu, block_gpu>>>(
                ncol, nlev, ngpt, nbnd, band_lims,
                bnd_flux_dn, bnd_flux_up, bnd_flux_net);
    }
}
