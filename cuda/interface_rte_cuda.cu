#include <openacc.h>
#include <cstdio>

#include "Types.h"
#include "rte_kernel_launcher_cuda.h"

template<typename T> T* acc_to_cuda(T* ptr) { return static_cast<T*>(acc_deviceptr(ptr)); }

extern "C"
{
}
