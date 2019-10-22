#! /usr/bin/env python
#
# This script runs runs the RTE+RRTMGP all-sky examples
#
import os, subprocess, glob
import argparse

rte_rrtmgp_dir   = os.path.join("..", "..")
# This files lives in $RRTMGP_DIR/examples/all-sky/
all_sky_dir      = "."
# Code should be run in the all_sky_dir directory

lw_gas_coeffs_file   = os.path.join(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-lw-g256-2018-12-04.nc")
sw_gas_coeffs_file   = os.path.join(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-sw-g224-2018-12-04.nc")

lw_clouds_coeff_file = os.path.join(rte_rrtmgp_dir, "extensions", "cloud_optics", "rrtmgp-cloud-optics-coeffs-lw.nc")
sw_clouds_coeff_file = os.path.join(rte_rrtmgp_dir, "extensions", "cloud_optics", "rrtmgp-cloud-optics-coeffs-sw.nc")

# In the local directory
all_sky_exe_name = os.path.join(all_sky_dir, "rrtmgp_allsky")
atmos_file       = os.path.join(all_sky_dir, "rrtmgp-allsky.nc")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Runs all-sky examples, resetting output.")
    parser.add_argument("--run_command", type=str, default="",
                        help="Prefix ('jsrun' etc.) for running commands. Use quote marks to enclose multi-part commands.")
    parser.add_argument("--ncol", type=int, default=128,
                        help="Number of cloudy columns to compute (every one will have the same clouds)")
    parser.add_argument("--nloops", type=int, default=1,
                        help="Number of times to compute 'nloops' cloudy columns")

    args = parser.parse_args()
    ncol_str   = '{0:5d}'.format(args.ncol)
    nloops_str = '{0:5d}'.format(args.nloops)
    if args.run_command:
        print ("using the run command")
        all_sky_exe_name = args.run_command + " " + all_sky_exe_name

    os.chdir(all_sky_dir)
    # Remove cloudy-sky fluxes from the file containing the atmospheric profiles
    subprocess.run(["git", "checkout", atmos_file])
    subprocess.run([all_sky_exe_name, atmos_file, lw_gas_coeffs_file, lw_clouds_coeff_file, ncol_str, nloops_str])
    subprocess.run([all_sky_exe_name, atmos_file, sw_gas_coeffs_file, sw_clouds_coeff_file, ncol_str, nloops_str])

# end main
