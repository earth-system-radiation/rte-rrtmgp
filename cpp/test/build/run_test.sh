#!/bin/bash

printf "Rebuilding\n\n"

make -j8  || exit -1

printf "Running allsky longwave tests\n\n"

cd allsky
cp ../*.nc .
./allsky rrtmgp-allsky.nc rrtmgp-data-lw-g256-2018-12-04.nc rrtmgp-cloud-optics-coeffs-lw.nc 1 1  || exit -1


printf "Running allsky shortwave tests\n\n"
./allsky rrtmgp-allsky.nc rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-cloud-optics-coeffs-sw.nc 1 1  || exit -1

