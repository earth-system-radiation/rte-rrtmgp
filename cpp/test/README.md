# RRTMGPXX testing
To run tests, execute the following:

    ./get_data.sh
    ./cmakescript.sh
    cd build && make && ctest

You will need to set some environment variables before running cmakescript.sh. Examples can be found in macbook.sh.

On summit, to test on the GPU, you will need to make sure the entire repo is contained within `/gpfs`. A full checkout and test workflow should look something like the following:
```
PROJECT=<project ID>
cd ${HOME}
git clone https://github.com/mrnorman/YAKL.git
git clone https://github.com/NVlabs/cub.git
cd /gpfs/alpine/${PROJECT}/scratch/$USER
git clone https://github.com/E3SM-Project/rte-rrtmgp.git
cd rte-rrtmgp/cpp/test
source summit_gpu.sh
./get_data.sh
./cmakescript.sh 
cd build 
make -j8 
bsub -I -W 00:05 -P ${PROJECT} -nnodes 1 make test
```
