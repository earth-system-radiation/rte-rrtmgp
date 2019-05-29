#! /usr/bin/env python
#
# This script runs runs RRTMGP on the RFMIP off-line test cases
#
import os

try:
    from subprocess import run as subprocess_call
except ImportError:
    from subprocess import call as subprocess_call

#
# Run the RFMIP example programs that computes fluxes from netCDF Garand
# atmosphere files
#
rte_rrtmgp_srcdir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                 "..", "..")
rte_rrtmgp_builddir = os.path.join("..", "..")

rfmip_builddir = os.path.join(rte_rrtmgp_builddir,
                              "examples", "rfmip-clear-sky")
conds_file = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"

rfmip_lw_exe_name = "rrtmgp_rfmip_lw"
rfmip_sw_exe_name = "rrtmgp_rfmip_sw"
print("Running RFMIP drivers")
# arguments are block size, input conditions, coefficient files,
# forcing index, physics index
subprocess_call([os.path.join(rfmip_builddir, rfmip_lw_exe_name),
                 "8", conds_file,
                 os.path.join(rte_rrtmgp_srcdir, "rrtmgp", "data",
                              "rrtmgp-data-lw-g256-2018-12-04.nc")])
subprocess_call([os.path.join(rfmip_builddir, rfmip_sw_exe_name),
                 "8", conds_file,
                 os.path.join(rte_rrtmgp_srcdir, "rrtmgp", "data",
                              "rrtmgp-data-sw-g224-2018-12-04.nc")])
