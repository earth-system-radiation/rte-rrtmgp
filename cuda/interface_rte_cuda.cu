#include <openacc.h>
#include <cstdio>
#include <exception>

#include "Types.h"
// CvH: zero_array should move to rte kernels.
#include "rrtmgp_kernel_launcher_cuda.h"
#include "rte_kernel_launcher_cuda.h"
#include "optical_props_kernel_launcher_cuda.h"


template<typename T> T* acc_to_cuda(T* ptr) { return static_cast<T*>(acc_deviceptr(ptr)); }


extern "C"
{
    // void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
    //               const Float* inc_flux_dir, const Float* mu0, Float* gpt_flux_dir);
    // void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, Float* gpt_flux_dn);
    // void apply_BC(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const Float* inc_flux_dif, Float* gpt_flux_dn);


    void rte_sw_solver_2stream(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1,
            Float* tau, Float* ssa, Float* g,
            Float* mu0,
            Float* sfc_alb_dir, Float* sfc_alb_dif,
            Float* inc_flux_dir,
            Float* flux_up, Float* flux_dn, Float* flux_dir,
            Bool* has_dif_bc, Float* inc_flux_dif,
            Bool* do_broadband, Float* flux_up_loc, Float* flux_dn_loc, Float* flux_dir_loc)
    {
        printf("CvH: sw_solver_2stream CUDA\n");
        rte_kernel_launcher_cuda::sw_solver_2stream(
                *ncol, *nlay, *ngpt, *top_at_1,
                acc_to_cuda(tau), acc_to_cuda(ssa), acc_to_cuda(g),
                acc_to_cuda(mu0),
                acc_to_cuda(sfc_alb_dir), acc_to_cuda(sfc_alb_dif),
                acc_to_cuda(inc_flux_dir),
                acc_to_cuda(flux_up), acc_to_cuda(flux_dn), acc_to_cuda(flux_dir),
                *has_dif_bc, acc_to_cuda(inc_flux_dif),
                *do_broadband, acc_to_cuda(flux_up_loc), acc_to_cuda(flux_dn_loc), acc_to_cuda(flux_dir_loc));
    }


    // void lw_solver_noscat_gaussquad(
    void rte_lw_solver_noscat(
            int* ncol, int* nlay, int* ngpt, Bool* top_at_1, int* nmus,
            Float* secants, Float* weights,
            Float* tau, Float* lay_source,
            Float* lev_source_inc, Float* lev_source_dec,
            Float* sfc_emis, Float* sfc_src,
            Float* inc_flux,
            Float* flux_up, Float* flux_dn,
            Bool* do_broadband, Float* flux_up_loc, Float* flux_dn_loc,
            Bool* do_jacobians, Float* sfc_src_jac, Float* flux_up_jac)
    {
        printf("CvH: lw_solver_noscat CUDA\n");
        rte_kernel_launcher_cuda::lw_solver_noscat_gaussquad(
                *ncol, *nlay, *ngpt, *top_at_1, *nmus,
                acc_to_cuda(secants), acc_to_cuda(weights),
                acc_to_cuda(tau), acc_to_cuda(lay_source),
                acc_to_cuda(lev_source_inc), acc_to_cuda(lev_source_dec),
                acc_to_cuda(sfc_emis), acc_to_cuda(sfc_src),
                acc_to_cuda(inc_flux),
                acc_to_cuda(flux_up), acc_to_cuda(flux_dn),
                *do_broadband, acc_to_cuda(flux_up_loc), acc_to_cuda(flux_dn_loc),
                *do_jacobians, acc_to_cuda(sfc_src_jac), acc_to_cuda(flux_up_jac));
    }

    // void lw_secants_array(
    //         const int ncol, const int ngpt, const int n_quad_angs, const int max_gauss_pts,
    //         const Float* Gauss_Ds, Float* secants);

    // OPTICAL PROPS.
    void rte_increment_1scalar_by_1scalar(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, const Float* tau_in)
    {
        optical_props_kernel_launcher_cuda::increment_1scalar_by_1scalar(
                *ncol, *nlay, *ngpt,
                tau_inout, tau_in);
    }


    void rte_increment_2stream_by_2stream(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            const Float* tau_in, const Float* ssa_in, const Float* g_in)
    {
        optical_props_kernel_launcher_cuda::increment_2stream_by_2stream(
                *ncol, *nlay, *ngpt,
                tau_inout, ssa_inout, g_inout,
                tau_in, ssa_in, g_in);
    }


    void rte_inc_1scalar_by_1scalar_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, const Float* tau_in,
            int* nbnd, int* band_lims_gpoint)
    {
        optical_props_kernel_launcher_cuda::inc_1scalar_by_1scalar_bybnd(
                *ncol, *nlay, *ngpt,
                tau_inout, tau_in,
                *nbnd, band_lims_gpoint);
    }


    void rte_inc_2stream_by_2stream_bybnd(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout,
            const Float* tau_in, const Float* ssa_in, const Float* g_in,
            int* nbnd, const int* band_lims_gpoint)
    {
        optical_props_kernel_launcher_cuda::inc_2stream_by_2stream_bybnd(
                *ncol, *nlay, *ngpt,
                tau_inout, ssa_inout, g_inout,
                tau_in, ssa_in, g_in,
                *nbnd, band_lims_gpoint);
    }

    void rte_delta_scale_2str_k(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout)
    {
        printf("CvH: delta_scale_2str_k CUDA\n");
        optical_props_kernel_launcher_cuda::delta_scale_2str_k(
            *ncol, *nlay, *ngpt,
            tau_inout, ssa_inout, g_inout);
    }


    void rte_delta_scale_2str_f_k(
            int* ncol, int* nlay, int* ngpt,
            Float* tau_inout, Float* ssa_inout, Float* g_inout, Float* f)
    {
        throw std::runtime_error("delta_scale_2str_f_k is not implemented in CUDA");
    }


}
