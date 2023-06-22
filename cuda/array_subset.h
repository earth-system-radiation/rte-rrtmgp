#ifndef ARRAY_SUBSET_H
#define ARRAY_SUBSET_H


#ifdef __CUDACC__
template<int N>
struct Subset_data
{
    int sub_strides[N];
    int strides[N];
    int starts[N];
    int offsets[N];
    bool do_spread[N];
};


template<typename T, int N> __global__
void subset_kernel(
        T* __restrict__ a_sub,
        const T* __restrict__ a,
        Subset_data<N> subset_data,
        const int ncells)
{
    const int idx_out = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx_out < ncells)
    {
        int ic = idx_out;
        int idx_in = 0;

        #pragma unroll
        for (int n=N-1; n>=0; --n)
        {
            const int idx_dim = subset_data.do_spread[n] ? 1 : ic / subset_data.sub_strides[n] + subset_data.starts[n];
            ic %= subset_data.sub_strides[n];

            idx_in += (idx_dim - subset_data.offsets[n] - 1) * subset_data.strides[n];
        }

        a_sub[idx_out] = a[idx_in];
    }
}


template<typename T> __global__
void fill_kernel(
        T* __restrict__ a,
        const T v,
        const int ncells)
{
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx < ncells)
    {
        a[idx] = v;
    }
}
#endif

#endif
