cp ../../extensions/cloud_optics/*.nc .
cp ../../rrtmgp/data/*.nc 

module load arm-forge/19.1.4

export NFHOME="/sw/spack-rhel6/netcdf-fortran-4.5.3-nqwtm5"
export NCHOME="/sw/spack-rhel6/netcdf-c-4.7.4-hlp365"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NFHOME/lib:$NCHOME/lib"

./rrtmgp_allsky rrtmgp-data-lw-g256-2018-12-04.nc rrtmgp-cloud-optics-coeffs-lw.nc  128 1
