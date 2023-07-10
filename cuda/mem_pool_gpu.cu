#include "mem_pool_gpu.h"


Memory_pool_gpu::Memory_pool_gpu(const std::map<std::size_t, std::size_t>& block_sizes_) : alloc_counter(0), alloc_bytes(0)
{
    float min_queue_split = 0.25;
    
    // Preprocess mapping: ensure large enough splitting between block sizes to enhance memory re-use
    std::map<std::size_t, std::size_t> split_block_sizes;
    if (min_queue_split > 0)
    {
        std::size_t num_ptrs = 0;
        for (auto it = block_sizes_.begin(); it != block_sizes_.end(); ++it)
        {
            num_ptrs += it->second;
            auto it2 = it;
            ++it2;
            if (it2 == block_sizes_.end() or it->first < (1. - min_queue_split) * (it2->first))
            {
                split_block_sizes[it->first] = num_ptrs;
                num_ptrs = 0;
            }
        }
    }
    else
    {
        split_block_sizes = block_sizes_;
    }

    std::size_t tot_size = 0;

    // First pass: determine big root block size
    for (auto it = split_block_sizes.begin(); it != split_block_sizes.end(); ++it)
    {
        tot_size += (it->first * it->second);
    }

    auto block = (char*)allocate(tot_size);
    
    // Second pass: split block into sub-block pointers
    char* ptr = block;
    for (auto it = split_block_sizes.begin(); it != split_block_sizes.end(); ++it)
    {
        auto block_size = it->first;
        std::list<void*> pointers;
        for (std::size_t i = 0; i < it->second ; ++i)
        {
            pointers.push_back((void*)ptr);
            size_registry[(void*)ptr] = block_size;
            ptr = &(*(ptr + block_size));
        }
        blocks[block_size] = pointers;
    }
}


Memory_pool_gpu::~Memory_pool_gpu()
{
    // Free raw pointers
    for (auto it = raw_pointers.begin(); it != raw_pointers.end(); ++it)
    {
        cudaFree(*it);
    }
    raw_pointers.clear();
    blocks.clear();
    size_registry.clear();
}


void* Memory_pool_gpu::allocate(std::size_t nbytes_)
{
    void* data_ptr = nullptr;
    int err = cudaMalloc((void **) &data_ptr, nbytes_);
    if (cudaSuccess != err)
    {
        printf("cudaMalloc failed attempting to allocate %lu bytes\n", nbytes_);
        throw 1;
    }
    alloc_counter++;
    alloc_bytes += nbytes_;
    raw_pointers.push_back(data_ptr);
    return data_ptr;
}


void* Memory_pool_gpu::acquire(std::size_t nbytes_)
{
    if (nbytes_ == 0)
    {
        return nullptr;
    }
    auto it = blocks.begin();
    // Find smallest block size large enough to fulfill request
    while(it != blocks.end() and it->first < nbytes_)
    {
        ++it;
    }
    void* ptr;
    // If the requested block size is too large, add it
    if (it == blocks.end())
    {
        ptr = allocate(nbytes_);
        size_registry[ptr] = nbytes_;
        blocks[nbytes_] = std::list<void*>();
    }
    else
    {
        // If the queue is exhausted, look forward for larger available blocks
        if (it->second.size() == 0)
        {
            auto it2 = it;
            while(it2 != blocks.end() and it2->second.size() == 0)
            {
                ++it2;
            }
            if (it2 == blocks.end())
            {
                ptr = allocate(it->first);
                size_registry[ptr] = it->first;
            }
            else
            {
                ptr = it2->second.back();
                it2->second.pop_back();
            }
        }
        else
        {
            ptr = it->second.back();
            it->second.pop_back();
        }
    }
    return ptr;
}


void Memory_pool_gpu::release(void* ptr_)
{
    std::size_t block_size = size_registry[ptr_];
    blocks[block_size].push_back(ptr_);
}
