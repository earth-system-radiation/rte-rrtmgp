#include <cstdio>
#include <openacc.h>
#include "Types.h"
#include "rrtmgp_kernel_launcher_cuda.h"

template<typename T> T* acc_to_cuda(T* ptr) { return static_cast<T*>(acc_deviceptr(ptr)); }

extern "C"
{
    void interpolation_cuda_(
            int* ncol, int* nlay,
            int* ngas, int* nflav, int* neta, int* npres, int* ntemp,
            int* flavor,
            Float* press_ref_log,
            Float* temp_ref,
            Float* press_ref_log_delta,
            Float* temp_ref_min,
            Float* temp_ref_delta,
            Float* press_ref_trop_log,
            Float* vmr_ref,
            Float* play,
            Float* tlay,
            Float* col_gas,
            int* jtemp,
            Float* fmajor, Float* fminor,
            Float* col_mix,
            Bool* tropo,
            int* jeta,
            int* jpress)
    {
        rrtmgp_kernel_launcher_cuda::interpolation(
                *ncol, *nlay,
                *ngas, *nflav, *neta, *npres, *ntemp,
                acc_to_cuda(flavor),
                acc_to_cuda(press_ref_log),
                acc_to_cuda(temp_ref),
                *press_ref_log_delta,
                *temp_ref_min,
                *temp_ref_delta,
                *press_ref_trop_log,
                acc_to_cuda(vmr_ref),
                acc_to_cuda(play),
                acc_to_cuda(tlay),
                acc_to_cuda(col_gas),
                acc_to_cuda(jtemp),
                acc_to_cuda(fmajor), acc_to_cuda(fminor),
                acc_to_cuda(col_mix),
                acc_to_cuda(tropo),
                acc_to_cuda(jeta),
                acc_to_cuda(jpress));
    }
}
