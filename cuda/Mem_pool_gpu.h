#ifndef MEM_POOL_GPU_H
#define MEM_POOL_GPU_H

#include <vector>
#include <list>
#include <map>
#include <tuple>
#include <iostream>


// Multi-queue memory pool for device pointers. 
class Memory_pool_gpu
{
    public:
        typedef std::map<std::size_t, std::list<void*>> memory_storage_type;

        static Memory_pool_gpu& get_instance()
        {
            std::map<std::size_t, std::size_t> empty_map;
            return init_instance(empty_map);
        }

        static Memory_pool_gpu& init_instance(const std::map<std::size_t, std::size_t>& block_sizes_)
        {
            static Memory_pool_gpu memory_pool(block_sizes_);
            return memory_pool;
        }

        Memory_pool_gpu(const std::map<std::size_t, std::size_t>& block_sizes_);

        ~Memory_pool_gpu();

        void* acquire(std::size_t nbytes_);

        void release(void* ptr_);

    private:
        void* allocate(std::size_t nbytes);

        memory_storage_type blocks;
        std::list<void*> raw_pointers;
        std::map<void*, std::size_t> size_registry; 

        std::size_t alloc_counter;
        std::size_t alloc_bytes;
};
#endif
