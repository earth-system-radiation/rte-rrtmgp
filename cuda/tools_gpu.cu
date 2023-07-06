#if defined(RTE_RRTMGP_GPU_MEMPOOL_CUDA)
#include <cstdint>
#include <cstdio>

static bool cuda_mempool_initialized = false;
void prepare_cuda_mempool() {
    if (cuda_mempool_initialized) return;

    printf("Setting up CUDA mempool.\n");
    cudaMemPool_t mempool;
    cudaDeviceGetDefaultMemPool(&mempool, 0);
    auto threshold = UINT64_MAX;
    cudaMemPoolSetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold);
    cuda_mempool_initialized = true;
}
#endif