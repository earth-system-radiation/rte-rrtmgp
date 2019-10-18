#!/usr/bin/env python
from __future__ import print_function

# J. Delamere, October 2019
# Stopgap program to add the optimal single angle fit coefficients to the
# longwave absorption coefficient file. Will use NCO tools to add the variable.
# >python ./create_optimal_angle_nc.py
# >cp rrtmgp-data-lw-g256-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04-a.nc
# >ncks -A -v optimal_angle_fit sec_t_fit_coeffs.nc
#  rrtmgp-data-lw-g256-2018-12-04-a.nc

# sec_t_fit_coeffs.csv was copied from
#
# standard libraries
import pandas as pd
import netCDF4 as nc

dat = pd.read_csv('sec_t_fit_coeffs.csv', header=0)

dataset = nc.Dataset("sec_t_fit_coeffs.nc", "w")
band = dataset.createDimension("bnd", 16)
fit_coeffs = dataset.createDimension("fit_coeffs", 2)
optimal_angle_fit = \
    dataset.createVariable("optimal_angle_fit", "f4",
                           ("bnd", "fit_coeffs"))
optimal_angle_fit[:] = dat

optimal_angle_fit.description = \
    'Coefficients for linear fit used in longwave \
     optimal angle RT calculation'
dataset.close()
