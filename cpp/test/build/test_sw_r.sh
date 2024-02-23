#!/bin/bash -ex

ncpdq -O -a -lay,-lev input_files/rrtmgp-allsky.nc input_files/rrtmgp-allsky-r.nc
$RUNCMD ./allsky/allsky input_files/rrtmgp-allsky-r.nc                \
                        input_files/rrtmgp-data-sw-g224-2018-12-04.nc \
                        input_files/rrtmgp-cloud-optics-coeffs-sw.nc  \
                        1                                             \
                        1
