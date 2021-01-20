#!/bin/bash

cd allsky
cp ../*.nc .
./allsky rrtmgp-allsky.nc rrtmgp-data-lw-g256-2018-12-04.nc rrtmgp-cloud-optics-coeffs-lw.nc 1 1 || exit -1

exit 0

