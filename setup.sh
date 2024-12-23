rm -rf build
# conda --version
# conda env create -f environment-dev.yml

FC=gfortran
FFLAGS='-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -fmodule-private -fimplicit-none -finit-real=nan'
RTE_ENABLE_SP=OFF
KERNEL_MODE=default
FAILURE_THRESHOLD='7.e-4'

cmake -S . -B build -G "Ninja" \
        -DCMAKE_Fortran_COMPILER=$FC \
        -DCMAKE_Fortran_FLAGS="$FFLAGS" \
        -DRTE_ENABLE_SP=$RTE_ENABLE_SP \
        -DKERNEL_MODE=$KERNEL_MODE \
        -DBUILD_TESTING=ON \
        -DFAILURE_THRESHOLD=$FAILURE_THRESHOLD \
        -DCMAKE_BUILD_TYPE=Release

# cmake --build build --target all --parallel

# The --test-dir option is available only starting CMake 3.20:
# ctest -V --test-dir build
