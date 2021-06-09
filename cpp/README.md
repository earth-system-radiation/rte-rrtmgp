# Welcome to the C++ RRTMGP

RRTMGP++ is ported to C++ using the Yet Another Kernel Launcher (YAKL) tool, and it was ported with readibility in mind, keeping it as similar to the original Fortran code as possible, including:
* Multi-dimensional arrays
* Column-major index ordering
* Array lower bounds that default to one but can be set to any value
* Fortran intrinsics
* Array slicing

All of this was retained while also achieve excellent performance on Nvidia and AMD GPUs (with support for Intel GPUs coming soon). For instance, on the Summit supercomputer, the six Nvidia Tesla V100 GPUs on a node run RRTMGP++'s longwave solver 45x faster than the two IBM Power9 CPUs on a node including all data transfer costs. The shortwave solver runs 30x faster. Even at very small workloads, the GPUs give sizeable speed-up because all fo the RRTMGP++ code (the `load_and_init` function is the only exception) runs continuously on the GPU with minimal traffic between device and host memory and with all kernels launched asynchronously to hide kernel launch latencies. The use of a pool allocator was critical for performance as well.

YAKL is a C++ performance portability libarary with nearly identical syntax to Kokkos, but it has better Fortran-like coding support. For more information, see https://github.com/mrnorman/YAKL .

To build and test the code, for example on Summit, please do the following:

* `cd /ccs/home/$USER`
* `git clone https://github.com/mrnorman/YAKL.git`
* `git clone https://github.com/NVlabs/cub.git`
* `cd /gpfs/alpine/[PROJECT_ID]/scratch/$USER`
* `git clone https://github.com/E3SM-Project/rte-rrtmgp.git`
* `cd rte-rrtmgp/cpp/test/build`
* `source summit_gpu.sh`
* `./get_data.sh`
* `./cmakescript.sh`
* `make -j8`
* `make test`

To build and test the code on your local GPU, please do the following:

* `cd /home/$USER`
* `git clone https://github.com/mrnorman/YAKL.git`
* `git clone https://github.com/NVlabs/cub.git`
* `git clone https://github.com/E3SM-Project/rte-rrtmgp.git`
* `cd rte-rrtmgp/cpp/test/build`
* edit `thatchroof_gpu.sh` to match your GPU's architecture
* `source thatchroof_gpu.sh`
* `./get_data.sh`
* `./cmakescript.sh`
* `make -j8`
* `make test`

To build and test the code on your local CPU, please do the following:

* `cd /home/$USER`
* `git clone https://github.com/mrnorman/YAKL.git`
* `git clone https://github.com/E3SM-Project/rte-rrtmgp.git`
* `cd rte-rrtmgp/cpp/test/build`
* `source thatchroof_cpu.sh`
* `./get_data.sh`
* `./cmakescript.sh`
* `make -j8`
* `make test`

