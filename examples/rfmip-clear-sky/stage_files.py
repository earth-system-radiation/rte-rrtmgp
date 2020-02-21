#! /usr/bin/env python
#
# This script downloads and creates files needed for the RFMIP off-line test cases
#
import os
import subprocess
import glob
import urllib.request
# This script invokes templ_scr_url, which uses Python modules time, uuid, argparse, netCDF4

#
# Download and/or create input files and output template files
#
rte_rrtmgp_dir  = os.path.join("..", "..")
rfmip_dir       = os.path.join(rte_rrtmgp_dir, "examples", "rfmip-clear-sky")
conds_file      = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
conds_url       = "http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/" + \
                  "atmos/fx/multiple/none/v20190401/" + conds_file
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
subprocess.run(["python3", templ_scr, "--source_id", "RTE-RRTMGP-181204"])
