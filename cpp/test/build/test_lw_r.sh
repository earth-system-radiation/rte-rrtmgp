#!/bin/bash -ex

ncpdq -O -a -lay,-lev input_files/rrtmgp-allsky.nc input_files/rrtmgp-allsky-r.nc
$RUNCMD ./allsky/allsky input_files/rrtmgp-allsky-r.nc                \
                        input_files/rrtmgp-data-lw-g256-2018-12-04.nc \
                        input_files/rrtmgp-cloud-optics-coeffs-lw.nc  \
                        1                                             \
                        1
