rm -rf build
# conda --version
# conda env create -f environment-dev.yml

export FC=gfortran
export FCFLAGS="-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -fmodule-private -fimplicit-none -finit-real=nan"
export RRTMGP_DATA_VERSION=v1.8.2
export FP_MODEL=DP
export RTE_BOOL=C
export RTE_KERNELS=default
export FAILURE_THRESHOLD=7.e-4

cmake -S . -B build -G "Ninja" \
        -DCMAKE_Fortran_COMPILER=$FC \
        -DCMAKE_Fortran_FLAGS="$FCFLAGS" \
        -DRRTMGP_DATA_VERSION=$RRTMGP_DATA_VERSION \
        -DPRECISION=$FP_MODEL \
        -DBOOL_TYPE=$RTE_BOOL \
        -DKERNEL_MODE=$RTE_KERNELS \
        -DENABLE_TESTS=ON \
        -DFAILURE_THRESHOLD=$FAILURE_THRESHOLD

# cmake --build build --config Release --target all --parallel

# ctest --test-dir build/ -V
