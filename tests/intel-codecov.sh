#!/bin/bash
ulimit -s hard
export FC=ifort
export FCFLAGS="-m64 -prof-gen=srcpos -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"
#
# Intel specific - where will the profiling files be generated?
#
export RRTMGP_ROOT=`cd ..;pwd`
export RRTMGP_BUILD=${RRTMGP_ROOT}/build
export PROF_DIR=$PWD
#
# Environment variables for netCDF Fortran and C installations
#
export NFHOME=${HOME}/Applications/${FC}
export NCHOME=/opt/local

#
# An Anaconda environent with modules needed for other python scripts
#   (xarray, matplotlib, ...)
#
source activate pangeo
cd ${PROF_DIR}
rm -rf *.dyn pgopti.* CODE_COVERAGE.html CodeCoverage/
#
# Build RTE+RRTMGP librarues
#
make -C ${RRTMGP_BUILD} clean
make -C ${RRTMGP_BUILD} -j 4 || exit 1

#
# Build and run RFMIP examples
#
cd ${RRTMGP_ROOT}/examples/rfmip-clear-sky || exit 1
export FCFLAGS+=" -I${RRTMGP_BUILD} -I${NFHOME}/include"
export LDFLAGS+=" -L${RRTMGP_BUILD} -L${NFHOME}/lib -L${NCHOME}/lib -lrte -lrrtmgp -lnetcdff -lnetcdf"
make clean || exit 1
make -j 4  || exit 1
python ./stage_files.py
python ./run-rfmip-examples.py
python ./compare-to-reference.py --fail=7.e-4
make clean
#
# Build and run all-sky examples
#
cd ${RRTMGP_ROOT}/examples/all-sky || exit 1
make clean || exit 1
make -j 4  || exit 1
python ./run-allsky-example.py
python ./compare-to-reference.py
make clean
#
# Build and run regression tests
#
cd ${PROF_DIR} || exit 1
make clean || exit 1
make -j 4  || exit 1
cp ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc test_atmospheres.nc
./clear_sky_regression test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
./clear_sky_regression test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
# Need to repeat for openacc-kernels

#
# Merge
#
cd ${PROF_DIR}
profmerge -a
echo "
mo_
~mo_simple_netcdf
~mo_rfmip_io
~mo_testing_io
~mo_garand_atmos_io
~mo_load" > intel_codecov_filter.txt
codecov -prj rte-rrtmgp -comp intel_codecov_filter.txt
rm intel_codecov_filter.txt
