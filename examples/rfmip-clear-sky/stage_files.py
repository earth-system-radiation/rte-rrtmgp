#! /usr/bin/env python
#
# This script downloads and creates files needed for the RFMIP off-line test
# cases
#
import glob
import os
import sys

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

try:
    from subprocess import run as subprocess_call
except ImportError:
    from subprocess import call as subprocess_call

#
# Download and/or create input files and output template files
#
conds_url_base = "http://aims3.llnl.gov/thredds/fileServer/user_pub_work/" \
                 "input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/" \
                 "atmos/fx/multiple/none/v20190401/"
conds_file = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"

templ_scr_url_base = "https://raw.githubusercontent.com/RobertPincus/" \
                     "RFMIP-IRF-Scripts/master/"
templ_scr = "generate-output-file-templates.py"

#
# Remove previous versions of files
#
for f in glob.glob("r??_Efx*.nc") + \
         glob.glob("multiple_input4MIPs_radiation_RFMIP*.nc"):
    os.remove(f)
#
# Download the profiles for RFMIP; make the empty output files
#
print("Dowloading RFMIP input files")
urlretrieve(conds_url_base + conds_file, conds_file)
print("Dowloading scripts for generating output templates")
urlretrieve(templ_scr_url_base + templ_scr, templ_scr)
subprocess_call([sys.executable, templ_scr,
                 "--source_id", "RTE-RRTMGP-181204"])

#
# Reference results
#
print("Downloading reference results")
ref_dir = os.getcwd() + "/reference/"
if not os.path.exists(ref_dir):
    os.makedirs(ref_dir)
urlretrieve(
    "https://owncloud.gwdg.de/index.php/s/kbhl3JOSccGtR0m/download",
    os.path.join(ref_dir, "rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
urlretrieve(
    "https://owncloud.gwdg.de/index.php/s/iFa28GFxRaNGKU1/download",
    os.path.join(ref_dir, "rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
urlretrieve(
    "https://owncloud.gwdg.de/index.php/s/uCemCHlGxbGK0gJ/download",
    os.path.join(ref_dir, "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
urlretrieve(
    "https://owncloud.gwdg.de/index.php/s/l8ZG28j9ttZWD9r/download",
    os.path.join(ref_dir, "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
