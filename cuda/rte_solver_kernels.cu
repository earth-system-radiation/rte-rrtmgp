#include <float.h>
#include "Types.h"

#ifndef kernel_tuner
const int loop_unroll_factor_nlay = 4;
#endif


template<typename TF> __device__ constexpr TF k_min();
template<> __device__ constexpr double k_min() { return 1.e-12; }
template<> __device__ constexpr float k_min() { return 1.e-4f; }


__global__
void lw_secants_array_kernel(
        const int ncol, const int ngpt, const int n_gauss_quad, const int max_gauss_pts,
        const Float* __restrict__ gauss_Ds, Float* __restrict__ secants)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;
    const int imu = blockIdx.z;

    if ( (icol < ncol) && (igpt < ngpt) && (imu < n_gauss_quad) )
    {
        const int idx_s = icol + igpt*ncol + imu*ncol*ngpt;
        const int idx_g = imu + (n_gauss_quad-1)*max_gauss_pts;

        secants[idx_s] = gauss_Ds[idx_g];
    }
}

 
__device__
void lw_transport_noscat_kernel(
        const int icol, const int igpt, const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
        const Float* __restrict__ tau, const Float* __restrict__ trans, const Float* __restrict__ sfc_albedo,
        const Float* __restrict__ source_dn, const Float* __restrict__ source_up, const Float* __restrict__ source_sfc,
        Float* __restrict__ radn_up, Float* __restrict__ radn_dn, const Float* __restrict__ source_sfc_jac, Float* __restrict__ radn_up_jac)
{
    if (top_at_1)
    {
        #pragma unroll loop_unroll_factor_nlay
        for (int ilev=1; ilev<(nlay+1); ++ilev)
        {
            const int idx1 = icol + ilev*ncol + igpt*ncol*(nlay+1);
            const int idx2 = icol + (ilev-1)*ncol + igpt*ncol*(nlay+1);
            const int idx3 = icol + (ilev-1)*ncol + igpt*ncol*nlay;
            radn_dn[idx1] = trans[idx3] * radn_dn[idx2] + source_dn[idx3];
        }

        const int idx_bot = icol + nlay*ncol + igpt*ncol*(nlay+1);
        const int idx2d = icol + igpt*ncol;
        radn_up[idx_bot] = radn_dn[idx_bot] * sfc_albedo[idx2d] + source_sfc[idx2d];
        radn_up_jac[idx_bot] = source_sfc_jac[idx2d];

        #pragma unroll loop_unroll_factor_nlay
        for (int ilev=nlay-1; ilev>=0; --ilev)
        {
            const int idx1 = icol + ilev*ncol + igpt*ncol*(nlay+1);
            const int idx2 = icol + (ilev+1)*ncol + igpt*ncol*(nlay+1);
            const int idx3 = icol + ilev*ncol + igpt*ncol*nlay;
            radn_up[idx1] = trans[idx3] * radn_up[idx2] + source_up[idx3];
            radn_up_jac[idx1] = trans[idx3] * radn_up_jac[idx2];
        }
    }
    else
    {
        #pragma unroll loop_unroll_factor_nlay
        for (int ilev=(nlay-1); ilev>=0; --ilev)
        {
            const int idx1 = icol + ilev*ncol + igpt*ncol*(nlay+1);
            const int idx2 = icol + (ilev+1)*ncol + igpt*ncol*(nlay+1);
            const int idx3 = icol + ilev*ncol + igpt*ncol*nlay;
            radn_dn[idx1] = trans[idx3] * radn_dn[idx2] + source_dn[idx3];
        }

        const int idx_bot = icol + igpt*ncol*(nlay+1);
        const int idx2d = icol + igpt*ncol;
        radn_up[idx_bot] = radn_dn[idx_bot] * sfc_albedo[idx2d] + source_sfc[idx2d];
        radn_up_jac[idx_bot] = source_sfc_jac[idx2d];

        #pragma unroll loop_unroll_factor_nlay
        for (int ilev=1; ilev<(nlay+1); ++ilev)
        {
            const int idx1 = icol + ilev*ncol + igpt*ncol*(nlay+1);
            const int idx2 = icol + (ilev-1)*ncol + igpt*ncol*(nlay+1);
            const int idx3 = icol + (ilev-1)*ncol + igpt*ncol*nlay;
            radn_up[idx1] = trans[idx3] * radn_up[idx2] + source_up[idx3];
            radn_up_jac[idx1] = trans[idx3] * radn_up_jac[idx2];;
        }
    }
}


__global__
void lw_solver_noscat_step_1_kernel(
        const int ncol, const int nlay, const int ngpt, const Float eps, const Bool top_at_1,
        const Float* __restrict__ D, const Float* __restrict__ weight, const Float* __restrict__ tau, const Float* __restrict__ lay_source,
        const Float* __restrict__ lev_source_inc, const Float* __restrict__ lev_source_dec, const Float* __restrict__ sfc_emis,
        const Float* __restrict__ sfc_src, Float* __restrict__ radn_up, Float* __restrict__ radn_dn,
        const Float* __restrict__ sfc_src_jac, Float* __restrict__ radn_up_jac, Float* __restrict__ tau_loc,
        Float* __restrict__ trans, Float* __restrict__ source_dn, Float* __restrict__ source_up,
        Float* __restrict__ source_sfc, Float* __restrict__ sfc_albedo, Float* __restrict__ source_sfc_jac)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx = icol + ilay*ncol + igpt*ncol*nlay;
        const int idx_D = icol + igpt*ncol;
        tau_loc[idx] = tau[idx] * D[idx_D];
        trans[idx] = exp(-tau_loc[idx]);

        const Float fact = (tau_loc[idx]>0. && (tau_loc[idx]*tau_loc[idx])>eps) ? (Float(1.) - trans[idx]) / tau_loc[idx] - trans[idx] : tau_loc[idx] * (Float(.5) - Float(1.)/Float(3.)*tau_loc[idx]);

        Float src_inc = (Float(1.) - trans[idx]) * lev_source_inc[idx] + Float(2.) * fact * (lay_source[idx]-lev_source_inc[idx]);
        Float src_dec = (Float(1.) - trans[idx]) * lev_source_dec[idx] + Float(2.) * fact * (lay_source[idx]-lev_source_dec[idx]);

        source_dn[idx] = top_at_1 ? src_inc : src_dec;
        source_up[idx] = top_at_1 ? src_dec : src_inc;
    }
}


__global__
void lw_solver_noscat_step_2_kernel(
        const int ncol, const int nlay, const int ngpt, const Float eps, const Bool top_at_1,
        const Float* __restrict__ D, const Float* __restrict__ weight, const Float* __restrict__ tau, const Float* __restrict__ lay_source,
        const Float* __restrict__ lev_source_inc, const Float* __restrict__ lev_source_dec, const Float* __restrict__ sfc_emis,
        const Float* __restrict__ sfc_src, Float* __restrict__ radn_up, Float* __restrict__ radn_dn,
        const Float* __restrict__ sfc_src_jac, Float* __restrict__ radn_up_jac, Float* __restrict__ tau_loc,
        Float* __restrict__ trans, Float* __restrict__ source_dn, Float* __restrict__ source_up,
        Float* __restrict__ source_sfc, Float* __restrict__ sfc_albedo, Float* __restrict__ source_sfc_jac)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        const int idx2d = icol + igpt*ncol;
        sfc_albedo[idx2d] = Float(1.) - sfc_emis[idx2d];
        source_sfc[idx2d] = sfc_emis[idx2d] * sfc_src[idx2d];
        source_sfc_jac[idx2d] = sfc_emis[idx2d] * sfc_src_jac[idx2d];

        const Float pi = acos(Float(-1.));
        const int idx_top = icol + (top_at_1 ? 0 : nlay)*ncol + igpt*ncol*(nlay+1);
        radn_dn[idx_top] = radn_dn[idx_top] / (Float(2.) * pi * weight[0]);


        lw_transport_noscat_kernel(
                icol, igpt, ncol, nlay, ngpt, top_at_1, tau, trans, sfc_albedo, source_dn,
                source_up, source_sfc, radn_up, radn_dn, source_sfc_jac, radn_up_jac);
    }
}


__global__
void lw_solver_noscat_step_3_kernel(
        const int ncol, const int nlay, const int ngpt, const Float eps, const Bool top_at_1,
        const Float* __restrict__ D, const Float* __restrict__ weight, const Float* __restrict__ tau, const Float* __restrict__ lay_source,
        const Float* __restrict__ lev_source_inc, const Float* __restrict__ lev_source_dec, const Float* __restrict__ sfc_emis,
        const Float* __restrict__ sfc_src, Float* __restrict__ radn_up, Float* __restrict__ radn_dn,
        const Float* __restrict__ sfc_src_jac, Float* __restrict__ radn_up_jac, Float* __restrict__ tau_loc,
        Float* __restrict__ trans, Float* __restrict__ source_dn, Float* __restrict__ source_up,
        Float* __restrict__ source_sfc, Float* __restrict__ sfc_albedo, Float* __restrict__ source_sfc_jac)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilev = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilev < (nlay+1)) && (igpt < ngpt) )
    {
        const Float pi = acos(Float(-1.));

        const int idx = icol + ilev*ncol + igpt*ncol*(nlay+1);
        radn_up[idx] *= Float(2.) * pi * weight[0];
        radn_dn[idx] *= Float(2.) * pi * weight[0];
        radn_up_jac[idx] *= Float(2.) * pi * weight[0];
    }
}


template<Bool top_at_1> __global__
void sw_adding_kernel(
        const int ncol, const int nlay, const int ngpt, const Bool _top_at_1,
        const Float* __restrict__ sfc_alb_dif, const Float* __restrict__ r_dif, const Float* __restrict__ t_dif,
        const Float* __restrict__ source_dn, const Float* __restrict__ source_up, const Float* __restrict__ source_sfc,
        Float* __restrict__ flux_up, Float* __restrict__ flux_dn, const Float* __restrict__ flux_dir,
        Float* __restrict__ albedo, Float* __restrict__ src, Float* __restrict__ denom)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        if (top_at_1)
        {
            const int sfc_idx_3d = icol + nlay*ncol + igpt*(nlay+1)*ncol;
            const int sfc_idx_2d = icol + igpt*ncol;
            albedo[sfc_idx_3d] = sfc_alb_dif[sfc_idx_2d];
            src[sfc_idx_3d] = source_sfc[sfc_idx_2d];

            #pragma unroll loop_unroll_factor_nlay
            for (int ilay=nlay-1; ilay >= 0; --ilay)
            {
                const int lay_idx  = icol + ilay*ncol + igpt*ncol*nlay;
                const int lev_idx1 = icol + ilay*ncol + igpt*ncol*(nlay+1);
                const int lev_idx2 = icol + (ilay+1)*ncol + igpt*ncol*(nlay+1);
                denom[lay_idx] = Float(1.)/(Float(1.) - r_dif[lay_idx] * albedo[lev_idx2]);
                albedo[lev_idx1] = r_dif[lay_idx] + t_dif[lay_idx] * t_dif[lay_idx]
                                                  * albedo[lev_idx2] * denom[lay_idx];
                src[lev_idx1] = source_up[lay_idx] + t_dif[lay_idx] * denom[lay_idx] *
                                (src[lev_idx2] + albedo[lev_idx2] * source_dn[lay_idx]);
            }
            const int top_idx = icol + igpt*(nlay+1)*ncol;
            flux_up[top_idx] = flux_dn[top_idx]*albedo[top_idx] + src[top_idx];

            for (int ilay=1; ilay < (nlay+2); ++ilay)
            {
                const int lev_idx1 = icol + ilay*ncol + igpt*(nlay+1)*ncol;
                const int lev_idx2 = icol + (ilay-1)*ncol + igpt*(nlay+1)*ncol;
                const int lay_idx = icol + (ilay-1)*ncol + igpt*(nlay)*ncol;
                if (ilay < (nlay+1)) {
                    flux_dn[lev_idx1] = (t_dif[lay_idx]*flux_dn[lev_idx2] +
                                         r_dif[lay_idx]*src[lev_idx1] +
                                         source_dn[lay_idx]) * denom[lay_idx];
                    flux_up[lev_idx1] = flux_dn[lev_idx1] * albedo[lev_idx1] + src[lev_idx1];
                }
                flux_dn[lev_idx2] += flux_dir[lev_idx2];
            }
        }
        else
        {
            const int sfc_idx_3d = icol + igpt*(nlay+1)*ncol;
            const int sfc_idx_2d = icol + igpt*ncol;
            albedo[sfc_idx_3d] = sfc_alb_dif[sfc_idx_2d];
            src[sfc_idx_3d] = source_sfc[sfc_idx_2d];

            #pragma unroll loop_unroll_factor_nlay
            for (int ilay=0; ilay<nlay; ++ilay)
            {
                const int lay_idx  = icol + ilay*ncol + igpt*ncol*nlay;
                const int lev_idx1 = icol + ilay*ncol + igpt*ncol*(nlay+1);
                const int lev_idx2 = icol + (ilay+1)*ncol + igpt*ncol*(nlay+1);
                denom[lay_idx] = Float(1.)/(Float(1.) - r_dif[lay_idx] * albedo[lev_idx1]);
                albedo[lev_idx2] = r_dif[lay_idx] + (t_dif[lay_idx] * t_dif[lay_idx] *
                                                     albedo[lev_idx1] * denom[lay_idx]);
                src[lev_idx2] = source_up[lay_idx] + t_dif[lay_idx]*denom[lay_idx]*
                                                     (src[lev_idx1]+albedo[lev_idx1]*source_dn[lay_idx]);
            }

            const int top_idx = icol + nlay*ncol + igpt*(nlay+1)*ncol;
            flux_up[top_idx] = flux_dn[top_idx] *albedo[top_idx] + src[top_idx];

            for (int ilay=nlay-1; ilay >= -1; --ilay)
            {
                const int lay_idx  = icol + ilay*ncol + igpt*nlay*ncol;
                const int lev_idx1 = icol + ilay*ncol + igpt*(nlay+1)*ncol;
                const int lev_idx2 = icol + (ilay+1)*ncol + igpt*(nlay+1)*ncol;

                    if (ilay >= 0)
                    {
                        flux_dn[lev_idx1] = (t_dif[lay_idx]*flux_dn[lev_idx2] +
                                             r_dif[lay_idx]*src[lev_idx1] +
                                             source_dn[lay_idx]) * denom[lay_idx];
                        flux_up[lev_idx1] = flux_dn[lev_idx1] * albedo[lev_idx1] + src[lev_idx1];
                    }

                    flux_dn[lev_idx2] += flux_dir[lev_idx2];
            }
        }
    }
}


template<Bool top_at_1> __global__
void sw_source_kernel(
        const int ncol, const int nlay, const int ngpt, const Bool _top_at_1,
        Float* __restrict__ r_dir, Float* __restrict__ t_dir, Float* __restrict__ t_noscat,
        const Float* __restrict__ sfc_alb_dir, Float* __restrict__ source_up, Float* __restrict__ source_dn,
        Float* __restrict__ source_sfc, Float* __restrict__ flux_dir)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        if (top_at_1)
        {
            for (int ilay=0; ilay<nlay; ++ilay)
            {
                const int idx_lay  = icol + ilay*ncol + igpt*nlay*ncol;
                const int idx_lev1 = icol + ilay*ncol + igpt*(nlay+1)*ncol;
                const int idx_lev2 = icol + (ilay+1)*ncol + igpt*(nlay+1)*ncol;
                source_up[idx_lay] = r_dir[idx_lay] * flux_dir[idx_lev1];
                source_dn[idx_lay] = t_dir[idx_lay] * flux_dir[idx_lev1];
                flux_dir[idx_lev2] = t_noscat[idx_lay] * flux_dir[idx_lev1];

            }
            const int sfc_idx = icol + igpt*ncol;
            const int flx_idx = icol + nlay*ncol + igpt*(nlay+1)*ncol;
            source_sfc[sfc_idx] = flux_dir[flx_idx] * sfc_alb_dir[icol];
        }
        else
        {
            for (int ilay=nlay-1; ilay>=0; --ilay)
            {
                const int idx_lay  = icol + ilay*ncol + igpt*nlay*ncol;
                const int idx_lev1 = icol + (ilay)*ncol + igpt*(nlay+1)*ncol;
                const int idx_lev2 = icol + (ilay+1)*ncol + igpt*(nlay+1)*ncol;
                source_up[idx_lay] = r_dir[idx_lay] * flux_dir[idx_lev2];   //uses updated flux_dir from previous iteration
                source_dn[idx_lay] = t_dir[idx_lay] * flux_dir[idx_lev2];   //uses updated flux_dir from previous
                flux_dir[idx_lev1] = t_noscat[idx_lay] * flux_dir[idx_lev2];//updates flux_dir for 0 to nlay-1

            }
            const int sfc_idx = icol + igpt*ncol;
            const int flx_idx = icol + igpt*(nlay+1)*ncol;
            source_sfc[sfc_idx] = flux_dir[flx_idx] * sfc_alb_dir[icol];
        }
    }
}

__global__
void apply_BC_kernel_lw(const int isfc, int ncol, const int nlay, const int ngpt, const Bool top_at_1, const Float* __restrict__ inc_flux, Float* __restrict__ flux_dn)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        const int idx_in  = icol + isfc*ncol + igpt*ncol*(nlay+1);
        const int idx_out = (top_at_1) ? icol + igpt*ncol*(nlay+1) : icol + nlay*ncol + igpt*ncol*(nlay+1);
        flux_dn[idx_out] = inc_flux[idx_in];
    }
}

__global__
void apply_BC_kernel(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const Float* __restrict__ inc_flux, Float* __restrict__ flux_dn)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;
    if ( (icol < ncol) && (igpt < ngpt) )
    {
        const int idx_out = icol + ((top_at_1 ? 0 : (nlay * ncol))) + (igpt * ncol * (nlay + 1));
        const int idx_in = icol + (igpt * ncol);
        flux_dn[idx_out] = inc_flux[idx_in];
    }
}

__global__
void apply_BC_kernel(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, const Float* __restrict__ inc_flux, const Float* __restrict__ factor, Float* __restrict__ flux_dn)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;
    if ( (icol < ncol) && (igpt < ngpt) )
    {
        const int idx_out = icol + ((top_at_1 ? 0 : (nlay * ncol))) + (igpt * ncol * (nlay + 1));
        const int idx_in = icol + (igpt * ncol);

        flux_dn[idx_out] = inc_flux[idx_in] * factor[icol];
    }
}

__global__
void apply_BC_kernel(const int ncol, const int nlay, const int ngpt, const Bool top_at_1, Float* __restrict__ flux_dn)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;
    if ( (icol < ncol) && (igpt < ngpt) )
    {
        const int idx_out = icol + ((top_at_1 ? 0 : (nlay * ncol))) + (igpt * ncol * (nlay + 1));
        flux_dn[idx_out] = Float(0.);
    }
}

__global__
void sw_2stream_kernel(
        const int ncol, const int nlay, const int ngpt, const Float tmin,
        const Float* __restrict__ tau, const Float* __restrict__ ssa,
        const Float* __restrict__ g, const Float* __restrict__ mu0,
        Float* __restrict__ r_dif, Float* __restrict__ t_dif,
        Float* __restrict__ r_dir, Float* __restrict__ t_dir,
        Float* __restrict__ t_noscat)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilay = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilay < nlay) && (igpt < ngpt) )
    {
        const int idx = icol + ilay*ncol + igpt*nlay*ncol;
        const Float mu0_inv = Float(1.)/mu0[icol];
        const Float gamma1 = (Float(8.) - ssa[idx] * (Float(5.) + Float(3.) * g[idx])) * Float(.25);
        const Float gamma2 = Float(3.) * (ssa[idx] * (Float(1.) -          g[idx])) * Float(.25);
        const Float gamma3 = (Float(2.) - Float(3.) * mu0[icol] *          g[idx])  * Float(.25);
        const Float gamma4 = Float(1.) - gamma3;

        const Float alpha1 = gamma1 * gamma4 + gamma2 * gamma3;
        const Float alpha2 = gamma1 * gamma3 + gamma2 * gamma4;

        const Float k = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), k_min<Float>()));
        const Float exp_minusktau = exp(-tau[idx] * k);
        const Float exp_minus2ktau = exp_minusktau * exp_minusktau;

        const Float rt_term = Float(1.) / (k      * (Float(1.) + exp_minus2ktau) +
                                     gamma1 * (Float(1.) - exp_minus2ktau));
        r_dif[idx] = rt_term * gamma2 * (Float(1.) - exp_minus2ktau);
        t_dif[idx] = rt_term * Float(2.) * k * exp_minusktau;
        t_noscat[idx] = exp(-tau[idx] * mu0_inv);

        const Float k_mu     = k * mu0[icol];
        const Float k_gamma3 = k * gamma3;
        const Float k_gamma4 = k * gamma4;

        const Float fact = (abs(Float(1.) - k_mu*k_mu) > tmin) ? Float(1.) - k_mu*k_mu : tmin;
        const Float rt_term2 = ssa[idx] * rt_term / fact;

        r_dir[idx] = rt_term2  * ((Float(1.) - k_mu) * (alpha2 + k_gamma3)   -
                                  (Float(1.) + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau -
                                   Float(2.) * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau * t_noscat[idx]);

        t_dir[idx] = -rt_term2 * ((Float(1.) + k_mu) * (alpha1 + k_gamma4) * t_noscat[idx]   -
                                  (Float(1.) - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * t_noscat[idx] -
                                   Float(2.) * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau);
    }
}

/*
template<typename Float>__global__
void sw_source_adding_kernel(const int ncol, const int nlay, const int ngpt, const Bool top_at_1,
                             const Float* __restrict__ sfc_alb_dir, const Float* __restrict__ sfc_alb_dif,
                             Float* __restrict__ r_dif, Float* __restrict__ t_dif,
                             Float* __restrict__ r_dir, Float* __restrict__ t_dir, Float* __restrict__ t_noscat,
                             Float* __restrict__ flux_up, Float* __restrict__ flux_dn, Float* __restrict__ flux_dir,
                             Float* __restrict__ source_up, Float* __restrict__ source_dn, Float* __restrict__ source_sfc,
                             Float* __restrict__ albedo, Float* __restrict__ src, Float* __restrict__ denom)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        sw_source_kernel(icol, igpt, ncol, nlay, top_at_1, r_dir, t_dir,
                         t_noscat, sfc_alb_dir, source_up, source_dn, source_sfc, flux_dir);

        sw_adding_kernel(icol, igpt, ncol, nlay, top_at_1, sfc_alb_dif,
                         r_dif, t_dif, source_dn, source_up, source_sfc,
                         flux_up, flux_dn, flux_dir, albedo, src, denom);
    }
}


__global__
void lw_solver_noscat_gaussquad_kernel(
        const int ncol, const int nlay, const int ngpt, const Float eps, const Bool top_at_1, const int nmus,
        const Float* __restrict__ secants, const Float* __restrict__ weights,
        const Float* __restrict__ tau, const Float* __restrict__ lay_source,
        const Float* __restrict__ lev_source_inc, const Float* __restrict__ lev_source_dec, const Float* __restrict__ sfc_emis,
        const Float* __restrict__ sfc_src, Float* __restrict__ radn_up, Float* __restrict__ radn_dn,
        const Float* __restrict__ sfc_src_jac, Float* __restrict__ radn_up_jac, Float* __restrict__ tau_loc,
        Float* __restrict__ trans, Float* __restrict__ source_dn, Float* __restrict__ source_up,
        Float* __restrict__ source_sfc, Float* __restrict__ sfc_albedo, Float* __restrict__ source_sfc_jac,
        Float* __restrict__ flux_up, Float* __restrict__ flux_dn, Float* __restrict__ flux_up_jac)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    // CvH ONLY TO MAKE IT COMPILE. REMOVE !!!!
    Float* ds = secants;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        lw_solver_noscat_kernel(
                icol, igpt, ncol, nlay, ngpt, eps, top_at_1, ds[0], weights[0], tau, lay_source,
                lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_up, flux_dn, sfc_src_jac,
                flux_up_jac, tau_loc, trans, source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

        const int top_level = top_at_1 ? 0 : nlay;
        apply_BC_kernel_lw(icol, igpt, top_level, ncol, nlay, ngpt, top_at_1, flux_dn, radn_dn);

        if (nmus > 1)
        {
            for (int imu=1; imu<nmus; ++imu)
            {
                lw_solver_noscat_kernel(
                        icol, igpt, ncol, nlay, ngpt, eps, top_at_1, ds[imu], weights[imu], tau, lay_source,
                        lev_source_inc, lev_source_dec, sfc_emis, sfc_src, radn_up, radn_dn, sfc_src_jac,
                        radn_up_jac, tau_loc, trans, source_dn, source_up, source_sfc, sfc_albedo, source_sfc_jac);

                for (int ilev=0; ilev<(nlay+1); ++ilev)
                {
                    const int idx = icol + ilev*ncol + igpt*ncol*(nlay+1);
                    flux_up[idx] += radn_up[idx];
                    flux_dn[idx] += radn_dn[idx];
                    flux_up_jac[idx] += radn_up_jac[idx];
                }
            }
        }
    }
}
*/


__global__
void add_fluxes_kernel(
        const int ncol, const int nlev, const int ngpt,
        const Float* __restrict__ radn_up, const Float* __restrict__ radn_dn, const Float* __restrict__ radn_up_jac,
        Float* __restrict__ flux_up, Float* __restrict__ flux_dn, Float* __restrict__ flux_up_jac)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int ilev = blockIdx.y*blockDim.y + threadIdx.y;
    const int igpt = blockIdx.z*blockDim.z + threadIdx.z;

    if ( (icol < ncol) && (ilev < nlev) && (igpt < ngpt) )
    {
        const int idx = icol + ilev*ncol + igpt*ncol*nlev;

        flux_up[idx] += radn_up[idx];
        flux_dn[idx] += radn_dn[idx];
        flux_up_jac[idx] += radn_up_jac[idx];
    }
}


template<typename TF> __device__ constexpr TF tmin();
template<> __forceinline__ __device__ constexpr double tmin() { return DBL_EPSILON; }
template<> __forceinline__ __device__ constexpr float tmin() { return FLT_EPSILON; }


__device__
void sw_2stream_function(
        const int icol, const int ilay, const int igpt,
        const int ncol, const int nlay, const int ngpt,
        const Float* __restrict__ tau, const Float* __restrict__ ssa,
        const Float* __restrict__ g, const Float* __restrict__ mu0,
        Float* __restrict__ r_dif, Float* __restrict__ t_dif,
        Float* __restrict__ r_dir, Float* __restrict__ t_dir,
        Float* __restrict__ t_noscat)
{
        const int idx = icol + ilay*ncol + igpt*nlay*ncol;

        const Float mu0_inv = Float(1.)/mu0[icol];
        const Float gamma1 = (Float(8.) - ssa[idx] * (Float(5.) + Float(3.) * g[idx])) * Float(.25);
        const Float gamma2 = Float(3.) * (ssa[idx] * (Float(1.) - g[idx])) * Float(.25);
        const Float gamma3 = (Float(2.) - Float(3.) * mu0[icol] * g[idx])  * Float(.25);
        const Float gamma4 = Float(1.) - gamma3;

        const Float alpha1 = gamma1 * gamma4 + gamma2 * gamma3;
        const Float alpha2 = gamma1 * gamma3 + gamma2 * gamma4;

        const Float k = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), k_min<Float>()));
        const Float exp_minusktau = exp(-tau[idx] * k);
        const Float exp_minus2ktau = exp_minusktau * exp_minusktau;

        const Float rt_term = Float(1.) / (k      * (Float(1.) + exp_minus2ktau) +
                                     gamma1 * (Float(1.) - exp_minus2ktau));
        r_dif[idx] = rt_term * gamma2 * (Float(1.) - exp_minus2ktau);
        t_dif[idx] = rt_term * Float(2.) * k * exp_minusktau;
        *t_noscat = exp(-tau[idx] * mu0_inv);

        const Float k_mu     = k * mu0[icol];
        const Float k_gamma3 = k * gamma3;
        const Float k_gamma4 = k * gamma4;

        const Float fact = (abs(Float(1.) - k_mu*k_mu) > tmin<Float>()) ? Float(1.) - k_mu*k_mu : tmin<Float>();
        const Float rt_term2 = ssa[idx] * rt_term / fact;

        *r_dir = rt_term2  * ((Float(1.) - k_mu) * (alpha2 + k_gamma3)   -
                                  (Float(1.) + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau -
                                   Float(2.) * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau * t_noscat[0]);

        *t_dir = -rt_term2 * ((Float(1.) + k_mu) * (alpha1 + k_gamma4) * t_noscat[0]   -
                                  (Float(1.) - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * t_noscat[0] -
                                   Float(2.) * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau);
        
        // fix thanks to peter ukkonen (see https://github.com/earth-system-radiation/rte-rrtmgp/pull/39#issuecomment-1026698541)
        *r_dir = max(tmin<Float>(), min(*r_dir, Float(1.0) - *t_noscat));
        *t_dir = max(tmin<Float>(), min(*t_dir, Float(1.0) - *t_noscat - *r_dir));
}


template<Bool top_at_1> __global__
void sw_source_2stream_kernel(
        const int ncol, const int nlay, const int ngpt,
        const Float* __restrict__ tau, const Float* __restrict__ ssa,
        const Float* __restrict__ g, const Float* __restrict__ mu0,
        Float* __restrict__ r_dif, Float* __restrict__ t_dif,
        const Float* __restrict__ sfc_alb_dir, Float* __restrict__ source_up, Float* __restrict__ source_dn,
        Float* __restrict__ source_sfc, Float* __restrict__ flux_dir)
{
    const int icol = blockIdx.x*blockDim.x + threadIdx.x;
    const int igpt = blockIdx.y*blockDim.y + threadIdx.y;

    if ( (icol < ncol) && (igpt < ngpt) )
    {
        if (top_at_1)
        {
            for (int ilay=0; ilay<nlay; ++ilay)
            {

                Float r_dir, t_dir, t_noscat;
                sw_2stream_function(icol, ilay, igpt,
                        ncol, nlay, ngpt,
                        tau, ssa, g, mu0,
                        r_dif, t_dif, &r_dir, &t_dir, &t_noscat);

                const int idx_lay  = icol + ilay*ncol + igpt*nlay*ncol;
                const int idx_lev1 = icol + ilay*ncol + igpt*(nlay+1)*ncol;
                const int idx_lev2 = icol + (ilay+1)*ncol + igpt*(nlay+1)*ncol;
                source_up[idx_lay] = r_dir * flux_dir[idx_lev1];
                source_dn[idx_lay] = t_dir * flux_dir[idx_lev1];
                flux_dir[idx_lev2] = t_noscat * flux_dir[idx_lev1];

            }
            const int sfc_idx = icol + igpt*ncol;
            const int flx_idx = icol + nlay*ncol + igpt*(nlay+1)*ncol;
            source_sfc[sfc_idx] = flux_dir[flx_idx] * sfc_alb_dir[icol];
        }
        else
        {
            for (int ilay=nlay-1; ilay>=0; --ilay)
            {
                Float r_dir, t_dir, t_noscat;
                sw_2stream_function(icol, ilay, igpt,
                        ncol, nlay, ngpt,
                        tau, ssa, g, mu0,
                        r_dif, t_dif, &r_dir, &t_dir, &t_noscat);

                const int idx_lay  = icol + ilay*ncol + igpt*nlay*ncol;
                const int idx_lev1 = icol + (ilay)*ncol + igpt*(nlay+1)*ncol;
                const int idx_lev2 = icol + (ilay+1)*ncol + igpt*(nlay+1)*ncol;
                source_up[idx_lay] = r_dir * flux_dir[idx_lev2];
                source_dn[idx_lay] = t_dir * flux_dir[idx_lev2];
                flux_dir[idx_lev1] = t_noscat * flux_dir[idx_lev2];

            }
            const int sfc_idx = icol + igpt*ncol;
            const int flx_idx = icol + igpt*(nlay+1)*ncol;
            source_sfc[sfc_idx] = flux_dir[flx_idx] * sfc_alb_dir[icol];
        }
    }
}

