#! /usr/bin/env python
#
# This script downloads and/or creates files needed for the RFMIP off-line test
# cases
#
import subprocess
import sys
import urllib.request
from pathlib import Path

#
# Download and/or create input files and output template files
#
rte_rrtmgp_dir = Path("..").joinpath("..")
rfmip_dir = rte_rrtmgp_dir.joinpath("examples", "rfmip-clear-sky")
conds_file = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
conds_url = ("http://aims3.llnl.gov/thredds/fileServer/user_pub_work/"
             "input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/" +
             "atmos/fx/multiple/none/v20190401/" + conds_file)
#
# The official RFMIP conditions are available from the ESFG, as above, but this
# fails from time to time, so we use the copy at Lamont-Doherty Earth
# Observatory, which we have to access via ftp(!)
#
conds_url = ("ftp://ftp.ldeo.columbia.edu/pub/robertp/rte-rrtmgp/"
             "continuous-integration/" + conds_file)
output_files = [f"r{wl}{d}_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"
                for wl in ['l', 's'] for d in ['u', 'd']]
#
# Remove previous versions of files
#
for f in output_files:
    Path(f).unlink(missing_ok=True)
Path(conds_file).unlink(missing_ok=True)

#
# Download the profiles for RFMIP; download or make the empty output files
#
print("Downloading RFMIP input files")
urllib.request.urlretrieve(conds_url, conds_file)

generate_templates = False
if generate_templates:
    print("Downloading scripts for generating output templates")
    urllib.request.urlretrieve(templ_scr_url, templ_scr)
    subprocess.run(
        [sys.executable, templ_scr, "--source_id", "RTE-RRTMGP-181204"])
else:
    print("Downloading output templates")
    for f in output_files:
        urllib.request.urlretrieve(conds_url.replace(conds_file, f), f)
