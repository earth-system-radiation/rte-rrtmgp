#ifndef TOOLS_GPU_H
#define TOOLS_GPU_H

#if defined(RTE_RRTMGP_GPU_MEMPOOL_OWN)
#include "mem_pool_gpu.h"
#endif
#include <cstdio>

#define cuda_safe_call(err) Tools_gpu::__cuda_safe_call(err, __FILE__, __LINE__)
#define cuda_check_error()  Tools_gpu::__cuda_check_error(__FILE__, __LINE__)
#define cuda_check_memory() Tools_gpu::__cuda_check_memory(__FILE__, __LINE__)

void prepare_cuda_mempool();

namespace Tools_gpu
{
    /* CUDA error checking.
       In debug mode, CUDACHECKS is defined and all kernel calls are checked with cudaCheckError().
       All CUDA api calls are always checked with cudaSafeCall() */

    // Wrapper to check for errors in CUDA api calls (e.g. cudaMalloc)
    inline void __cuda_safe_call(cudaError err, const char *file, const int line)
    {
        if (cudaSuccess != err)
        {
            printf("cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
            throw 1;
        }
    }

    // Function to check for errors in CUDA kernels. Call directly after kernel.
    inline void __cuda_check_error(const char *file, const int line)
    {
        #ifdef CUDACHECKS
        cudaError err = cudaGetLastError();
        if (cudaSuccess != err)
        {
            printf("cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString( err ) );
            throw 1;
        }

        err = cudaDeviceSynchronize();
        if (cudaSuccess != err)
        {
            printf("cudaCheckError() with sync failed at %s:%i : %s\n", file, line, cudaGetErrorString( err ) );
            throw 1;
        }
        #endif
    }

    // Check the memory usage.
    inline void __cuda_check_memory(const char *file, const int line)
    {
        #ifdef CUDACHECKS
        size_t free_byte, total_byte ;

        cudaError err = cudaMemGetInfo( &free_byte, &total_byte ) ;

        if ( cudaSuccess != err ){

            printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(err) );
            throw 1;

        }

        double used_db = (double)total_byte - (double)free_byte ;

        printf("GPU memory usage at %s:%i: %f MB\n", file, line, used_db/(1024.0*1024.0));
        #endif
    }

    template<typename T>
    T* allocate_gpu(int length)
    {
        T* data_ptr = nullptr;

        #if defined(RTE_RRTMGP_GPU_MEMPOOL_CUDA)
        prepare_cuda_mempool();
        cuda_safe_call(cudaMallocAsync((void **) &data_ptr, length*sizeof(T), 0));
        #elif defined(RTE_RRTMGP_GPU_MEMPOOL_OWN)
        data_ptr = (T*)(Memory_pool_gpu::get_instance().acquire(length*sizeof(T)));
        #else
        cuda_safe_call(cudaMalloc((void **) &data_ptr, length*sizeof(T)));
        #endif
        return data_ptr;
    }

    template<typename T>
    void free_gpu(T*& data_ptr)
    {
        #if defined(RTE_RRTMGP_GPU_MEMPOOL_CUDA)
        cuda_safe_call(cudaFreeAsync(data_ptr, 0));
        #elif defined(RTE_RRTMGP_GPU_MEMPOOL_OWN)
        Memory_pool_gpu::get_instance().release((void*)data_ptr);
        #else
        cuda_safe_call(cudaFree(data_ptr));
        #endif
        data_ptr = nullptr;
    }
}
#endif
