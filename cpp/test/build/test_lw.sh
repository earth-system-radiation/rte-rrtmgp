#!/bin/bash


$RUNCMD ./allsky/allsky input_files/rrtmgp-allsky.nc                  \
                        input_files/rrtmgp-data-lw-g256-2018-12-04.nc \
                        input_files/rrtmgp-cloud-optics-coeffs-lw.nc  \
                        1000                                          \
                        1                                             || exit -1

$RUNCMD ./allsky_fortran/allsky_fortran                               \
                        input_files/rrtmgp-allsky.nc                  \
                        input_files/rrtmgp-data-lw-g256-2018-12-04.nc \
                        input_files/rrtmgp-cloud-optics-coeffs-lw.nc  \
                        1000                                          \
                        1                                             || exit -1

$RUNCMD ./allsky_fortran_openacc/allsky_fortran_openacc               \
                        input_files/rrtmgp-allsky.nc                  \
                        input_files/rrtmgp-data-lw-g256-2018-12-04.nc \
                        input_files/rrtmgp-cloud-optics-coeffs-lw.nc  \
                        1000                                          \
                        1                                             || exit -1

