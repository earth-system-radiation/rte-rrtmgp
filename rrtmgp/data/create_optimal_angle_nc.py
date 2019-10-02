#!/usr/bin/env python
from __future__ import print_function

# standard libraries
import os, sys, argparse, glob, shutil
import numpy as np
import pandas as pd
import configparser as ConfigParser
import netCDF4 as nc
import subprocess as sub

dat = pd.read_csv('sec_t_fit_coeffs.csv', header=None)

dataset = nc.Dataset("sec_t_fit_coeffs.nc", "w")
band = dataset.createDimension("bnd", 16)
fit_coeffs = dataset.createDimension("fit_coeffs", 2)
optimal_single_angle_fit = dataset.createVariable("optimal_single_angle_fit","f4",("bnd","fit_coeffs"))
optimal_single_angle_fit[:] = dat
optimal_single_angle_fit.description = 'Optimized angle coefficients for linear fit'
dataset.close()
