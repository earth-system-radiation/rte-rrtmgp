#! /usr/bin/env python
#
# This script runs runs RRTMGP on the RFMIP off-line test cases
#
import os, subprocess, glob
from shutil import copy2
import urllib.request
# Will be needed by scripts to generate output file templates
from netCDF4 import Dataset
import time, uuid, argparse
import json

#
# Download and/or create input files and output template files
#
rte_rrtmgp_dir  = os.path.join("..", "..")
rfmip_dir       = os.path.join(rte_rrtmgp_dir, "examples", "rfmip-clear-sky")
conds_file      = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
conds_url       = "http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/" + \
                  "atmos/fx/multiple/none/v20190401/" + conds_file
# This branch expect conditions file version 1.1; will need to update when the drivers move to v1.2 by removing the next 3 lines
conds_file      = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-1_none.nc"
conds_url       = "http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-1/" + \
                  "atmos/fx/multiple/none/v20190204/" + conds_file
templ_scr       = "generate-output-file-templates.py"
templ_scr_url   = "https://raw.githubusercontent.com/RobertPincus/RFMIP-IRF-Scripts/master/" + templ_scr
#
# Remove previous versions of files
#
for f in glob.glob("r??_Efx*.nc"): os.remove(f)
for f in glob.glob("multiple_input4MIPs_radiation_RFMIP*.nc"): os.remove(f)
#
# Download the profiles for RFMIP; make the empty output files
#
print("Dowloading RFMIP input files")
urllib.request.urlretrieve(conds_url,     conds_file)
print("Dowloading scripts for generating output templates")
urllib.request.urlretrieve(templ_scr_url, templ_scr)
#%run -i generate-output-file-templates.py --source_id RTE-RRTMGP-181204
subprocess.run(["python3", templ_scr, "--source_id", "RTE-RRTMGP-181204"])

#
# Build and run the RFMIP example programs that computes fluxes from netCDF Garand atmosphere files
#
print("Building RFMIP drivers")
subprocess.run(["export RRTMGP_ROOT=" + rte_rrtmgp_dir + "; " + \
                "cd " + rfmip_dir + "; make "], shell=True)
rfmip_lw_exe_name = "rrtmgp_rfmip_lw"
rfmip_sw_exe_name = "rrtmgp_rfmip_sw"

print("Running RFMIP drivers")
# arguments are block size, input conditions, coefficient files, forcing index, physics index
subprocess.run([os.path.join(rfmip_dir, rfmip_lw_exe_name), "8", conds_file, os.path.join(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-lw-g256-2018-12-04.nc")])
subprocess.run([os.path.join(rfmip_dir, rfmip_sw_exe_name), "8", conds_file, os.path.join(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-sw-g224-2018-12-04.nc")])
