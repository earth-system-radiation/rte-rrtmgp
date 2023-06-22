#include <chrono>

#include "rte_solver_kernels_cuda.h"
#include "tools_gpu.h"
#include "tuner.h"

#include <iomanip>


namespace
{
    #include "rte_solver_kernels.cu"
}


namespace Rte_solver_kernels_cuda
{
    void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
                  const Float* inc_flux_dir, const Float* mu0, Float* gpt_flux_dir)
    {
        const int block_col = 32;
        const int block_gpt = 32;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);

        dim3 grid_gpu(grid_col, grid_gpt);
        dim3 block_gpu(block_col, block_gpt);

        apply_BC_kernel<<<grid_gpu, block_gpu>>>(ncol, nlay, ngpt, top_at_1, inc_flux_dir, mu0, gpt_flux_dir);
    }


    void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, Float* gpt_flux_dn)
    {
        const int block_col = 32;
        const int block_gpt = 32;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);

        dim3 grid_gpu(grid_col, grid_gpt);
        dim3 block_gpu(block_col, block_gpt);
        apply_BC_kernel<<<grid_gpu, block_gpu>>>(ncol, nlay, ngpt, top_at_1, gpt_flux_dn);
    }


    void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const Float* inc_flux_dif, Float* gpt_flux_dn)
    {
        const int block_col = 32;
        const int block_gpt = 32;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);

        dim3 grid_gpu(grid_col, grid_gpt);
        dim3 block_gpu(block_col, block_gpt);

        apply_BC_kernel<<<grid_gpu, block_gpu>>>(ncol, nlay, ngpt, top_at_1, inc_flux_dif, gpt_flux_dn);
    }


    void lw_secants_array(
            const int ncol, const int ngpt, const int n_gauss_quad, const int max_gauss_pts,
            const Float* gauss_Ds, Float* secants)
    {
        const int block_col = 32;
        const int block_gpt = 32;

        const int grid_col = ncol/block_col + (ncol%block_col > 0);
        const int grid_gpt = ngpt/block_gpt + (ngpt%block_gpt > 0);

        dim3 grid_gpu(grid_col, grid_gpt, n_gauss_quad);
        dim3 block_gpu(block_col, block_gpt, 1);

        lw_secants_array_kernel<<<grid_gpu, block_gpu>>>(
                ncol, ngpt, n_gauss_quad, max_gauss_pts,
                gauss_Ds, secants);
    }


    void lw_solver_noscat(
            const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const int nmus,
            const Float* secants, const Float* weights,
            const Float* tau, const Float* lay_source,
            const Float* lev_source_inc, const Float* lev_source_dec,
            const Float* sfc_emis, const Float* sfc_src,
            const Float* inc_flux,
            Float* flux_up, Float* flux_dn,
            const Bool do_broadband, Float* flux_up_loc, Float* flux_dn_loc,
            const Bool do_jacobians, const Float* sfc_src_jac, Float* flux_up_jac)
    {
        Float eps = std::numeric_limits<Float>::epsilon();

        const int flx_size = ncol*(nlay+1)*ngpt;
        const int opt_size = ncol*nlay*ngpt;
        const int sfc_size = ncol*ngpt;;

        Float* source_sfc = Tools_gpu::allocate_gpu<Float>(sfc_size);
        Float* source_sfc_jac = Tools_gpu::allocate_gpu<Float>(sfc_size);
        Float* sfc_albedo = Tools_gpu::allocate_gpu<Float>(sfc_size);
        Float* tau_loc = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* trans = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* source_dn = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* source_up = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* radn_dn = Tools_gpu::allocate_gpu<Float>(flx_size);
        Float* radn_up = Tools_gpu::allocate_gpu<Float>(flx_size);
        Float* radn_up_jac = Tools_gpu::allocate_gpu<Float>(flx_size);

        const int block_col2d = 64;
        const int block_gpt2d = 2;

        const int grid_col2d = ncol/block_col2d + (ncol%block_col2d > 0);
        const int grid_gpt2d = ngpt/block_gpt2d + (ngpt%block_gpt2d > 0);

        dim3 grid_gpu2d(grid_col2d, grid_gpt2d);
        dim3 block_gpu2d(block_col2d, block_gpt2d);

        const int block_col3d = 96;
        const int block_lay3d = 1;
        const int block_gpt3d = 1;

        const int grid_col3d = ncol/block_col3d + (ncol%block_col3d > 0);
        const int grid_lay3d = (nlay+1)/block_lay3d + ((nlay+1)%block_lay3d > 0);
        const int grid_gpt3d = ngpt/block_gpt3d + (ngpt%block_gpt3d > 0);

        dim3 grid_gpu3d(grid_col3d, grid_lay3d, grid_gpt3d);
        dim3 block_gpu3d(block_col3d, block_lay3d, block_gpt3d);

        const int top_level = top_at_1 ? 0 : nlay;


        // Upper boundary condition.
        if (inc_flux == nullptr)
            Rte_solver_kernels_cuda::apply_BC(ncol, nlay, ngpt, top_at_1, flux_dn);
        else
            Rte_solver_kernels_cuda::apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, flux_dn);


        // Step 1.
        Tuner_map& tunings = Tuner::get_map();

        dim3 grid_1, block_1;

        if (tunings.count("lw_step_1") == 0)
        {
                   Float* flux_up_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
                   Float* flux_dn_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
                   Float* flux_up_jac_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);

            std::tie(grid_1, block_1) = tune_kernel(
                    "lw_step_1",
                    dim3(ncol, nlay, ngpt),
                    {8, 16, 24, 32, 48, 64, 96, 128, 256, 512, 1024}, {1, 2, 4, 8}, {1},
                    lw_solver_noscat_step_1_kernel,
                    ncol, nlay, ngpt, eps, top_at_1,
                    secants, weights, tau, lay_source,
                    lev_source_inc, lev_source_dec,
                    sfc_emis, sfc_src, flux_up_tmp, flux_dn_tmp, sfc_src_jac,
                    flux_up_jac_tmp, tau_loc, trans,
                    source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);
                   
                    Tools_gpu::free_gpu<Float>(flux_up_tmp);
                    Tools_gpu::free_gpu<Float>(flux_dn_tmp);
                    Tools_gpu::free_gpu<Float>(flux_up_jac_tmp);

            tunings["lw_step_1"].first = grid_1;
            tunings["lw_step_1"].second = block_1;
        }
        else
        {
            grid_1 = tunings["lw_step_1"].first;
            block_1 = tunings["lw_step_1"].second;
        }

        lw_solver_noscat_step_1_kernel<<<grid_1, block_1>>>(
                ncol, nlay, ngpt, eps, top_at_1,
                secants, weights, tau, lay_source,
                lev_source_inc, lev_source_dec,
                sfc_emis, sfc_src, flux_up, flux_dn, sfc_src_jac,
                flux_up_jac, tau_loc, trans,
                source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);


        // Step 2.
        dim3 grid_2, block_2;

        if (tunings.count("lw_step_2") == 0)
        {
            Float* flux_up_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            Float* flux_dn_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            Float* flux_up_jac_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            
            std::tie(grid_2, block_2) = tune_kernel(
                    "lw_step_2",
                    dim3(ncol, ngpt),
                    {64, 128, 256, 384, 512, 768, 1024}, {1, 2, 4}, {1},
                    lw_solver_noscat_step_2_kernel,
                    ncol, nlay, ngpt, eps, top_at_1,
                    secants, weights, tau, lay_source,
                    lev_source_inc, lev_source_dec,
                    sfc_emis, sfc_src, flux_up_tmp, flux_dn_tmp, sfc_src_jac,
                    flux_up_jac_tmp, tau_loc, trans,
                    source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

            Tools_gpu::free_gpu<Float>(flux_up_tmp);
            Tools_gpu::free_gpu<Float>(flux_dn_tmp);
            Tools_gpu::free_gpu<Float>(flux_up_jac_tmp);
            
            tunings["lw_step_2"].first = grid_2;
            tunings["lw_step_2"].second = block_2;
        }
        else
        {
            grid_2 = tunings["lw_step_2"].first;
            block_2 = tunings["lw_step_2"].second;
        }

        lw_solver_noscat_step_2_kernel<<<grid_2, block_2>>>(
                ncol, nlay, ngpt, eps, top_at_1,
                secants, weights, tau, lay_source,
                lev_source_inc, lev_source_dec,
                sfc_emis, sfc_src, flux_up, flux_dn, sfc_src_jac,
                flux_up_jac, tau_loc, trans,
                source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);


        // Step 3.
        dim3 grid_3, block_3;

        if (tunings.count("lw_step_3") == 0)
        {
            Float* flux_up_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            Float* flux_dn_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            Float* flux_up_jac_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            
            std::tie(grid_3, block_3) = tune_kernel(
                    "lw_step_3",
                    dim3(ncol, nlay+1, ngpt),
                    {8, 16, 24, 32, 48, 64, 96, 128, 256}, {1, 2, 4, 8}, {1},
                    lw_solver_noscat_step_3_kernel,
                    ncol, nlay, ngpt, eps, top_at_1,
                    secants, weights, tau, lay_source,
                    lev_source_inc, lev_source_dec,
                    sfc_emis, sfc_src, flux_up_tmp, flux_dn_tmp, sfc_src_jac,
                    flux_up_jac_tmp, tau_loc, trans,
                    source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

            Tools_gpu::free_gpu<Float>(flux_up_tmp);
            Tools_gpu::free_gpu<Float>(flux_dn_tmp);
            Tools_gpu::free_gpu<Float>(flux_up_jac_tmp);

            tunings["lw_step_3"].first = grid_3;
            tunings["lw_step_3"].second = block_3;
        }
        else
        {
            grid_3 = tunings["lw_step_3"].first;
            block_3 = tunings["lw_step_3"].second;
        }

        lw_solver_noscat_step_3_kernel<<<grid_3, block_3>>>(
                ncol, nlay, ngpt, eps, top_at_1,
                secants, weights, tau, lay_source,
                lev_source_inc, lev_source_dec,
                sfc_emis, sfc_src, flux_up, flux_dn, sfc_src_jac,
                flux_up_jac, tau_loc, trans,
                source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

        apply_BC_kernel_lw<<<grid_gpu2d, block_gpu2d>>>(top_level, ncol, nlay, ngpt, top_at_1, flux_dn, radn_dn);

        if (nmus > 1)
        {
            for (int imu=1; imu<nmus; ++imu)
            {
                throw std::runtime_error("Not implemented due to lacking test case");
                /*
                lw_solver_noscat_step_1_kernel<<<grid_1, block_1>>>(
                        ncol, nlay, ngpt, eps, top_at_1,
                        secants+imu, weights+imu, tau, lay_source,
                        lev_source_inc, lev_source_dec,
                        sfc_emis, sfc_src, radn_up, radn_dn, sfc_src_jac,
                        radn_up_jac, tau_loc, trans,
                        source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

                lw_solver_noscat_step_2_kernel<<<grid_2, block_2>>>(
                        ncol, nlay, ngpt, eps, top_at_1,
                        secants+imu, weights+imu, tau, lay_source,
                        lev_source_inc, lev_source_dec,
                        sfc_emis, sfc_src,
                        radn_up, radn_dn, sfc_src_jac,
                        radn_up_jac, tau_loc, trans,
                        source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

                lw_solver_noscat_step_3_kernel<<<grid_3, block_3>>>(
                        ncol, nlay, ngpt, eps, top_at_1,
                        secants+imu, weights+imu, tau, lay_source,
                        lev_source_inc, lev_source_dec,
                        sfc_emis, sfc_src, radn_up, radn_dn, sfc_src_jac,
                        radn_up_jac, tau_loc, trans,
                        source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

                add_fluxes_kernel<<<grid_gpu3d, block_gpu3d>>>(
                        ncol, nlay+1, ngpt,
                        radn_up, radn_dn, radn_up_jac,
                        flux_up, flux_dn, flux_up_jac);
                        */
            }
        }

        Tools_gpu::free_gpu(source_sfc);
        Tools_gpu::free_gpu(source_sfc_jac);
        Tools_gpu::free_gpu(sfc_albedo);
        Tools_gpu::free_gpu(tau_loc);
        Tools_gpu::free_gpu(trans);
        Tools_gpu::free_gpu(source_dn);
        Tools_gpu::free_gpu(source_up);
        Tools_gpu::free_gpu(radn_dn);
        Tools_gpu::free_gpu(radn_up);
        Tools_gpu::free_gpu(radn_up_jac);
    }


    void sw_solver_2stream(
            const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
            const Float* tau, const Float* ssa, const Float* g,
            const Float* mu0,
            const Float* sfc_alb_dir, const Float* sfc_alb_dif,
            const Float* inc_flux_dir,
            Float* flux_up, Float* flux_dn, Float* flux_dir,
            const Bool has_dif_bc, const Float* inc_flux_dif,
            const Bool do_broadband, Float* flux_up_loc, Float* flux_dn_loc, Float* flux_dir_loc)
    {
        const int opt_size = ncol*nlay*ngpt;
        const int alb_size = ncol*ngpt;
        const int flx_size = ncol*(nlay+1)*ngpt;

        Float* r_dif = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* t_dif = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* r_dir = nullptr;
        Float* t_dir = nullptr;
        Float* t_noscat = nullptr;
        Float* source_up = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* source_dn = Tools_gpu::allocate_gpu<Float>(opt_size);
        Float* source_sfc = Tools_gpu::allocate_gpu<Float>(alb_size);
        Float* albedo = Tools_gpu::allocate_gpu<Float>(flx_size);
        Float* src = Tools_gpu::allocate_gpu<Float>(flx_size);
        Float* denom = Tools_gpu::allocate_gpu<Float>(opt_size);

        dim3 grid_source(ncol, ngpt), block_source;

        // Step0. Upper boundary condition. At this stage, flux_dn contains the diffuse radiation only.
        Rte_solver_kernels_cuda::apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dir, mu0, flux_dir);
        if (inc_flux_dif == nullptr)
            Rte_solver_kernels_cuda::apply_BC(ncol, nlay, ngpt, top_at_1, flux_dn);
        else
            Rte_solver_kernels_cuda::apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif, flux_dn);


        // Step 1.
        Tuner_map& tunings = Tuner::get_map();

        if (tunings.count("sw_source_2stream_kernel") == 0)
        {
            if (top_at_1)
            {
                std::tie(grid_source, block_source) = tune_kernel(
                        "sw_source_2stream_kernel",
                        dim3(ncol, ngpt),
                        {8, 16, 32, 64, 96, 128, 256, 384, 512, 768, 1024}, {1, 2, 4, 8, 16}, {1},
                        sw_source_2stream_kernel<1>,
                        ncol, nlay, ngpt, tau, ssa, g, mu0, r_dif, t_dif,
                        sfc_alb_dir, source_up, source_dn, source_sfc, flux_dir);
            }
            else
            {
                std::tie(grid_source, block_source) = tune_kernel(
                        "sw_source_2stream_kernel",
                        dim3(ncol, ngpt),
                        {8, 16, 32, 64, 96, 128, 256, 384, 512, 768, 1024}, {1, 2, 4, 8, 16}, {1},
                        sw_source_2stream_kernel<0>,
                        ncol, nlay, ngpt, tau, ssa, g, mu0, r_dif, t_dif,
                        sfc_alb_dir, source_up, source_dn, source_sfc, flux_dir);
            }

            tunings["sw_source_2stream_kernel"].first = grid_source;
            tunings["sw_source_2stream_kernel"].second = block_source;
        }
        else
        {
            grid_source = tunings["sw_source_2stream_kernel"].first;
            block_source = tunings["sw_source_2stream_kernel"].second;
        }

        if (top_at_1)
        {
            sw_source_2stream_kernel<1><<<grid_source, block_source>>>(
                    ncol, nlay, ngpt, tau, ssa, g, mu0, r_dif, t_dif,
                    sfc_alb_dir, source_up, source_dn, source_sfc, flux_dir);
        }
        else
        {
            sw_source_2stream_kernel<0><<<grid_source, block_source>>>(
                    ncol, nlay, ngpt, tau, ssa, g, mu0, r_dif, t_dif,
                    sfc_alb_dir, source_up, source_dn, source_sfc, flux_dir);
        }


        // Step 2.
        dim3 grid_adding, block_adding;

        if (tunings.count("sw_adding") == 0)
        {
            Float* flux_up_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            Float* flux_dn_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
            Float* flux_dir_tmp = Tools_gpu::allocate_gpu<Float>(ngpt*(nlay+1)*ncol);
                
            if (top_at_1)
            {
                std::tie(grid_adding, block_adding) = tune_kernel(
                        "sw_adding",
                        dim3(ncol, ngpt),
                        {8, 16, 32, 64, 96, 128, 256, 384, 512, 768, 1024}, {1, 2, 4, 8, 16}, {1},
                        sw_adding_kernel<1>,
                        ncol, nlay, ngpt, top_at_1,
                        sfc_alb_dif, r_dif, t_dif,
                        source_dn, source_up, source_sfc,
                        flux_up_tmp, flux_dn_tmp, flux_dir_tmp, albedo, src, denom);
            
            }
            else
            {
                std::tie(grid_adding, block_adding) = tune_kernel(
                        "sw_adding",
                        dim3(ncol, ngpt),
                        {8, 16, 32, 64, 96, 128, 256, 384, 512, 768, 1024}, {1, 2, 4, 8, 16}, {1},
                        sw_adding_kernel<0>,
                        ncol, nlay, ngpt, top_at_1,
                        sfc_alb_dif, r_dif, t_dif,
                        source_dn, source_up, source_sfc,
                        flux_up_tmp, flux_dn_tmp, flux_dir_tmp, albedo, src, denom);
            }
            
            Tools_gpu::free_gpu<Float>(flux_up_tmp);
            Tools_gpu::free_gpu<Float>(flux_dn_tmp);
            Tools_gpu::free_gpu<Float>(flux_dir_tmp);

            tunings["sw_adding"].first = grid_adding;
            tunings["sw_adding"].second = block_adding;
        }
        else
        {
            grid_adding = tunings["sw_adding"].first;
            block_adding = tunings["sw_adding"].second;
        }

        if (top_at_1)
        {
            sw_adding_kernel<1><<<grid_adding, block_adding>>>(
                ncol, nlay, ngpt, top_at_1,
                sfc_alb_dif, r_dif, t_dif,
                source_dn, source_up, source_sfc,
                flux_up, flux_dn, flux_dir, albedo, src, denom);
        }
        else
        {
            sw_adding_kernel<0><<<grid_adding, block_adding>>>(
                        ncol, nlay, ngpt, top_at_1,
                        sfc_alb_dif, r_dif, t_dif,
                        source_dn, source_up, source_sfc,
                        flux_up, flux_dn, flux_dir, albedo, src, denom);
        }

        Tools_gpu::free_gpu(r_dif);
        Tools_gpu::free_gpu(t_dif);
        Tools_gpu::free_gpu(source_up);
        Tools_gpu::free_gpu(source_dn);
        Tools_gpu::free_gpu(source_sfc);
        Tools_gpu::free_gpu(albedo);
        Tools_gpu::free_gpu(src);
        Tools_gpu::free_gpu(denom);
    }
}
