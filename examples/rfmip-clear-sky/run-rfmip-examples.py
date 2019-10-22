#! /usr/bin/env python
#
# This script runs runs RRTMGP on the RFMIP off-line test cases
#
import os, subprocess, glob
import argparse

#
# Run the RFMIP example programs that computes fluxes from netCDF Garand atmosphere files
#
rte_rrtmgp_dir  = os.path.join("..", "..")
rfmip_dir       = "."
# Code should be run in the rfmip_dir directory

conds_file      = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
lw_gas_coeffs_file   = os.path.join(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-lw-g256-2018-12-04.nc")
sw_gas_coeffs_file   = os.path.join(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-sw-g224-2018-12-04.nc")

rfmip_lw_exe_name = "./rrtmgp_rfmip_lw"
rfmip_sw_exe_name = "./rrtmgp_rfmip_sw"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Runs all-sky examples, resetting output.")
    parser.add_argument("--run_command", type=str, default="",
                        help="Prefix ('jsrun' etc.) for running commands. Use quote marks to enclose multi-part commands.")
    parser.add_argument("--block_size", type=int, default=8,
                        help="Number of columns to compute at a time. Must be a factor of 1800 (ncol*nexp)")

    args = parser.parse_args()
    block_size_str   = '{0:4d}'.format(args.block_size)
    if args.run_command:
        print ("using the run command")
        rfmip_lw_exe_name = args.run_command + " " + rfmip_lw_exe_name
        rfmip_sw_exe_name = args.run_command + " " + rfmip_sw_exe_name

    print("Running " + rfmip_lw_exe_name + block_size_str + " " + conds_file + " " + lw_gas_coeffs_file)
    # arguments are block size, input conditions, coefficient files, forcing index, physics index
    subprocess.run([rfmip_lw_exe_name, block_size_str, conds_file, lw_gas_coeffs_file])
    subprocess.run([rfmip_sw_exe_name, block_size_str, conds_file, sw_gas_coeffs_file])
