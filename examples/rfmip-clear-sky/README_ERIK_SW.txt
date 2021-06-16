cp ../../rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc coefficients_sw.nc

module load arm-forge/19.1.4

export NFHOME="/sw/spack-rhel6/netcdf-fortran-4.5.3-nqwtm5"
export NCHOME="/sw/spack-rhel6/netcdf-c-4.7.4-hlp365"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NFHOME/lib:$NCHOME/lib"

ddt rrtmgp_rfmip_lw
