#!/bin/bash

cd allsky
cp ../*.nc .
ncpdq -O -a -lay,-lev rrtmgp-allsky.nc rrtmgp-allsky-r.nc
$RUNCMD ./allsky rrtmgp-allsky-r.nc rrtmgp-data-lw-g256-2018-12-04.nc rrtmgp-cloud-optics-coeffs-lw.nc 1 1 || exit -1
exit 0
