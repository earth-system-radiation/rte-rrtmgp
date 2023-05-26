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

#ifndef ARRAY_H
#define ARRAY_H

#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <utility>

#ifdef __CUDACC__
#include "tools_gpu.h"
template<typename T, int N> class Array_gpu;
#endif


template<int N>
inline std::array<int, N> calc_strides(const std::array<int, N>& dims)
{
    std::array<int, N> strides;
    strides[0] = 1;
    for (int i=1; i<N; ++i)
        strides[i] = strides[i-1]*dims[i-1];

    return strides;
}

template<int N>
inline int calc_index(
        const std::array<int, N>& indices,
        const std::array<int, N>& strides,
        const std::array<int, N>& offsets)
{
    int sum = 0;
    for (int i=0; i<N; ++i)
        sum += (indices[i]-offsets[i]-1)*strides[i];

    return sum;
}

template<int N>
inline std::array<int, N> calc_indices(
        int index, const std::array<int, N>& strides, const std::array<int, N>& offsets)
{
    std::array<int, N> indices;

    for (int i=N-1; i>=1; --i)
    {
        indices[i] = index / strides[i];
        index %= strides[i];
    }
    indices[0] = index;

    for (int i=0; i<N; ++i)
        indices[i] += offsets[i] + 1;

    return indices;
}

template<int N>
inline int product(const std::array<int, N>& array)
{
    int product = array[0];
    for (int i=1; i<N; ++i)
        product *= array[i];

    return product;
}

template<typename T, int N>
class Array
{
    public:
        // Create an empty array, without dimensions.
        Array() :
            dims({}),
            ncells(0)
        {}

        // Create an array of zeros with given dimensions.
        Array(const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(ncells),
            strides(calc_strides<N>(dims)),
            offsets({})
        {}

        // Create an array from copying the contents of an std::vector.
        Array(const std::vector<T>& data, const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(data.begin(), data.begin() + ncells), // Do not copy beyond the end.
            strides(calc_strides<N>(dims)),
            offsets({})
        {} // CvH Do we need to size check data?

        // Create an array from moving the contents of an std::vector.
        Array(std::vector<T>&& data, const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(std::move(data)),
            strides(calc_strides<N>(dims)),
            offsets({})
        {} // CvH Do we need to size check data?

        // Define the default copy constructor and assignment operator.
        Array(const Array<T, N>&) = default;
        Array<T,N>& operator=(const Array<T, N>&) = default; // CvH does this one need empty checking?

        // Implement the move constructor to set ncells back to 0.
        Array(Array<T, N>&& array) :
            dims(std::exchange(array.dims, {})),
            ncells(std::exchange(array.ncells, 0)),
            data(std::move(array.data)),
            strides(std::exchange(array.strides, {})),
            offsets(std::exchange(array.offsets, {}))
        {}

        #ifdef __CUDACC__
        Array(const Array_gpu<T, N>& array_gpu) :
            dims(array_gpu.dims),
            ncells(array_gpu.ncells),
            data(ncells),
            strides(array_gpu.strides),
            offsets(array_gpu.offsets)
        {
            cuda_safe_call(cudaMemcpy(data.data(), array_gpu.ptr(), ncells*sizeof(T), cudaMemcpyDeviceToHost));
        }
        #endif

        inline void set_offsets(const std::array<int, N>& offsets)
        {
            this->offsets = offsets;
        }

        inline std::array<int, N> get_dims() const { return dims; }

        inline void set_dims(const std::array<int, N>& dims)
        {
            if (this->ncells != 0)
                throw std::runtime_error("Only arrays of size 0 can be resized");

            this->dims = dims;
            ncells = product<N>(dims);
            data.resize(ncells);
            strides = calc_strides<N>(dims);
            offsets = {};
        }

        inline std::vector<T>& v() { return data; }
        inline const std::vector<T>& v() const { return data; }

        inline T* ptr() { return data.data(); }
        inline const T* ptr() const { return data.data(); }

        inline int size() const { return ncells; }

        // inline std::array<int, N> find_indices(const T& value) const
        // {
        //     int pos = std::find(data.begin(), data.end(), value) - data.begin();
        //     return calc_indices<N>(pos, strides, offsets);
        // }

        inline T max() const
        {
            return *std::max_element(data.begin(), data.end());
        }

        inline T min() const
        {
            return *std::min_element(data.begin(), data.end());
        }

        inline void operator=(std::vector<T>&& data)
        {
            // CvH check size.
            this->data = data;
        }

        inline T& operator()(const std::array<int, N>& indices)
        {
            const int index = calc_index<N>(indices, strides, offsets);
            return data[index];
        }

        inline T operator()(const std::array<int, N>& indices) const
        {
            const int index = calc_index<N>(indices, strides, offsets);
            return data[index];
        }

        inline int dim(const int i) const { return dims[i-1]; }
        inline bool is_empty() const { return ncells == 0; }

        inline Array<T, N> subset(
                const std::array<std::pair<int, int>, N> ranges) const
        {
            // Calculate the dimension sizes based on the range.
            std::array<int, N> subdims;
            std::array<bool, N> do_spread;

            for (int i=0; i<N; ++i)
            {
                subdims[i] = ranges[i].second - ranges[i].first + 1;
                // CvH how flexible / tolerant are we?
                do_spread[i] = (dims[i] == 1);
            }

            // Create the array and fill it with the subset.
            Array<T, N> a_sub(subdims);
            for (int i=0; i<a_sub.ncells; ++i)
            {
                std::array<int, N> indices;
                int ic = i;
                for (int n=N-1; n>0; --n)
                {
                    indices[n] = do_spread[n] ? 1 : ic / a_sub.strides[n] + ranges[n].first;
                    ic %= a_sub.strides[n];
                }
                indices[0] = do_spread[0] ? 1 : ic + ranges[0].first;
                a_sub.data[i] = (*this)(indices);
            }

            return a_sub;
        }

        inline void fill(const T value)
        {
            std::fill(data.begin(), data.end(), value);
        }

        inline void dump(const std::string& name) const
        {
            std::string file_name = name + ".bin";
            std::ofstream binary_file(file_name, std::ios::out | std::ios::trunc | std::ios::binary);

            if (binary_file)
                binary_file.write(reinterpret_cast<const char*>(data.data()), ncells*sizeof(T));
            else
            {
                std::string error = "Cannot write file \"" + file_name + "\"";
                throw std::runtime_error(error);
            }
        }

    private:
        std::array<int, N> dims;
        int ncells;
        std::vector<T> data;
        std::array<int, N> strides;
        std::array<int, N> offsets;

        #ifdef __CUDACC__
        friend class Array_gpu<T, N>;
        #endif
};


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


template<typename T, int N>
class Array_gpu
{
    public:
        // Create an empty array, without dimensions.
        Array_gpu() :
            dims({}),
            ncells(0),
            data_ptr(nullptr),
            is_view(false)
        {}

        #ifdef __CUDACC__
        ~Array_gpu()
        {
            if (is_view)
                data_ptr = nullptr;
            else
                Tools_gpu::free_gpu(data_ptr);
        }
        #endif

        #ifdef __CUDACC__
        Array_gpu& operator=(const Array_gpu<T, N>& array)
        {
            if ( !(this->ncells == array.size() || (this->ncells == 0 && data_ptr == nullptr)) )
                throw std::runtime_error("Initialised arrays can not be resized");

            dims = array.dims;
            ncells = array.ncells;
            strides = array.strides;
            offsets = array.offsets;

            if (array.is_view)
            {
                is_view = true;
                data_ptr = array.data_ptr;
            }
            else if (this->ncells == 0)
            {
                is_view = false;
                data_ptr = Tools_gpu::allocate_gpu<T>(ncells);
                cuda_safe_call(cudaMemcpy(data_ptr, array.ptr(), ncells*sizeof(T), cudaMemcpyDeviceToDevice));
            }
            else
            {
                cuda_safe_call(cudaMemcpy(data_ptr, array.ptr(), ncells*sizeof(T), cudaMemcpyDeviceToDevice));
            }
            return (*this);
        }
        #endif

        #ifdef __CUDACC__
        Array_gpu& operator=(Array_gpu<T, N>&& array)
        {
            if ( !(this->ncells == array.size() || (this->ncells == 0 && data_ptr == nullptr)) )
                throw std::runtime_error("initialised arrays can not be resized");

            if (this->ncells > 0)
                Tools_gpu::free_gpu(data_ptr);

            dims = std::exchange(array.dims, {});
            ncells = std::exchange(array.ncells, 0);
            data_ptr = std::exchange(array.data_ptr, nullptr);
            strides = std::exchange(array.strides, {});
            offsets = std::exchange(array.offsets, {});
            is_view = std::exchange(array.is_view, false);

            return (*this);
        }
        #endif

        #ifdef __CUDACC__
        Array_gpu(const Array_gpu<T, N>& array) :
            dims(array.dims),
            ncells(array.ncells),
            data_ptr(nullptr),
            strides(array.strides),
            offsets(array.offsets),
            is_view(false)
        {
            if (array.is_view)
            {
                is_view = true;
                data_ptr = array.data_ptr;
            }
            else
            {
                data_ptr = Tools_gpu::allocate_gpu<T>(ncells);
                cuda_safe_call(cudaMemcpy(data_ptr, array.ptr(), ncells*sizeof(T), cudaMemcpyDeviceToDevice));
            }
        }
        #endif

        #ifdef __CUDACC__
        Array_gpu(Array_gpu<T, N>&& array) :
            dims(std::exchange(array.dims, {})),
            ncells(std::exchange(array.ncells, 0)),
            data_ptr(std::exchange(array.data_ptr, nullptr)),
            strides(std::exchange(array.strides, {})),
            offsets(std::exchange(array.offsets, {})),
            is_view(std::exchange(array.is_view, false))
        {
        }
        #endif

        #ifdef __CUDACC__
        // Create an array of zeros with given dimensions.
        Array_gpu(const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data_ptr(nullptr),
            strides(calc_strides<N>(dims)),
            offsets({}),
            is_view(false)
        {
            data_ptr = Tools_gpu::allocate_gpu<T>(ncells);
        }
        #endif

        #ifdef __CUDACC__
        // Create an array that is a view.
        Array_gpu(T* ptr, const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data_ptr(ptr),
            strides(calc_strides<N>(dims)),
            offsets({}),
            is_view(true)
        {
        }
        #endif

        #ifdef __CUDACC__
        Array_gpu(const Array<T, N>& array) :
            dims(array.dims),
            ncells(array.ncells),
            data_ptr(nullptr),
            strides(array.strides),
            offsets(array.offsets),
            is_view(false)
        {
            data_ptr = Tools_gpu::allocate_gpu<T>(ncells);
            cuda_safe_call(cudaMemcpy(data_ptr, array.ptr(), ncells*sizeof(T), cudaMemcpyHostToDevice));
        }
        #endif

        inline void set_offsets(const std::array<int, N>& offsets)
        {
            this->offsets = offsets;
        }

        inline std::array<int, N> get_dims() const { return dims; }

        #ifdef __CUDACC__
        inline void fill(const T value)
        {
            constexpr int block_ncells = 64;
            const int grid_ncells = this->ncells/block_ncells + (this->ncells%block_ncells > 0);

            dim3 block_gpu(block_ncells);
            dim3 grid_gpu(grid_ncells);

            fill_kernel<<<grid_gpu, block_gpu>>>(data_ptr, value, ncells);
        }
        #endif

        #ifdef __CUDACC__
        inline void set_data(const Array<T, N>& array)
        {
            data_ptr = Tools_gpu::allocate_gpu<T>(ncells);
            cuda_safe_call(cudaMemcpy(data_ptr, array.ptr(), ncells*sizeof(T), cudaMemcpyHostToDevice));
        }
        #endif

        #ifdef __CUDACC__
        inline void set_dims(const std::array<int, N>& dims)
        {
            if ( !(this->ncells == 0 && data_ptr == nullptr) )
                throw std::runtime_error("Only arrays with uninitialized pointers can be resized");

            this->dims = dims;
            ncells = product<N>(dims);
            data_ptr = Tools_gpu::allocate_gpu<T>(ncells);
            strides = calc_strides<N>(dims);
            offsets = {};
        }
        #endif

        inline void copy(const std::array<int, N>& indices, Array_gpu<T, N>& input, const std::array<int, N>& indices_input) const
        {
            #ifdef __CUDACC__
            const int index = calc_index<N>(indices, strides, offsets);
            const int index_in =  calc_index<N>(indices_input, input.strides, input.offsets);
            cuda_safe_call(cudaMemcpy(data_ptr + index, input.ptr() + index_in, sizeof(T), cudaMemcpyDeviceToDevice));
            #endif
        }

        inline void insert(const std::array<int, N>& indices, const T value) const
        {
            #ifdef __CUDACC__
            const int index = calc_index<N>(indices, strides, offsets);
            cuda_safe_call(cudaMemcpy(data_ptr + index, &value, sizeof(T), cudaMemcpyHostToDevice));
            #endif
        }

        inline T* ptr() { return data_ptr; }
        inline const T* ptr() const { return data_ptr; }

        inline int size() const { return ncells; }

        #ifdef __CUDACC__
        inline T operator()(const std::array<int, N>& indices) const
        {
            const int index = calc_index<N>(indices, strides, offsets);
            T value;
            cuda_safe_call(cudaMemcpy(&value, data_ptr + index, sizeof(T), cudaMemcpyDeviceToHost));
            return value;
        }
        #endif

        inline int dim(const int i) const { return dims[i-1]; }

        #ifdef __CUDACC__
        inline Array_gpu<T, N> subset(
                const std::array<std::pair<int, int>, N> ranges) const
        {
            // Calculate the dimension sizes based on the range.
            std::array<int, N> subdims;
            std::array<int, N> block_corners;

            for (int i=0; i<N; ++i)
            {
                subdims[i] = ranges[i].second - ranges[i].first + 1;
                // CvH how flexible / tolerant are we?
                block_corners[i] = ranges[i].first;
            }

            // Create the array and fill it with the subset.
            Array_gpu<T, N> a_sub(subdims);

            return subset_copy(a_sub, block_corners);
        }

        inline Array_gpu<T, N> subset_copy(Array_gpu<T, N>& a_sub,
                const std::array<int, N>& block_corners) const
        {
            Subset_data<N> subset_data;

            for (int i=0; i<N; ++i)
            {
                subset_data.sub_strides[i] = a_sub.strides[i];
                subset_data.strides[i] = strides[i];
                subset_data.starts[i] = block_corners[i];
                subset_data.offsets[i] = offsets[i];
                subset_data.do_spread[i] = (dims[i] == 1);
            }

            constexpr int block_ncells = 64;
            const int grid_ncells = a_sub.ncells/block_ncells + (a_sub.ncells%block_ncells > 0);

            dim3 block_gpu(block_ncells);
            dim3 grid_gpu(grid_ncells);

            subset_kernel<<<grid_gpu, block_gpu>>>(a_sub.data_ptr, data_ptr, subset_data, a_sub.ncells);

            return a_sub;
        }
        #endif

        inline void dump(const std::string& name) const
        {
            std::string file_name = name + ".bin";
            std::ofstream binary_file(file_name, std::ios::out | std::ios::trunc | std::ios::binary);

            Array<T, N> array_dump(*this);

            if (binary_file)
                binary_file.write(reinterpret_cast<char*>(array_dump.ptr()), ncells*sizeof(T));
            else
            {
                std::string error = "Cannot write file \"" + file_name + "\"";
                throw std::runtime_error(error);
            }
        }

    private:
        std::array<int, N> dims;
        int ncells;
        T* data_ptr;
        std::array<int, N> strides;
        std::array<int, N> offsets;
        bool is_view;

        friend class Array<T, N>;
};

template<typename T, int N>
bool any_vals_outside(const Array<T, N>& array, const T lower_limit, const T upper_limit)
{
    return std::any_of(
            array.v().begin(),
            array.v().end(),
            [lower_limit, upper_limit](T val){ return (val < lower_limit) || (val > upper_limit); });
}

template<typename T, int N>
bool any_vals_less_than(const Array<T, N>& array, const T lower_limit)
{
    return std::any_of(
            array.v().begin(),
            array.v().end(),
            [lower_limit](T val){ return (val < lower_limit); });
}
#endif
