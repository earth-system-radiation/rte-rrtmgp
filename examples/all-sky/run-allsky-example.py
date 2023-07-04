#! /usr/bin/env python
#
# This script runs runs the RTE+RRTMGP all-sky examples
#
import argparse
import glob
import os
import shutil
import subprocess

rte_rrtmgp_dir = os.environ["RRTMGP_DATA"]
# This files lives in $RRTMGP_ROOT/examples/all-sky/
all_sky_dir = "."
# Code should be run in the all_sky_dir directory

lw_gas_coeffs_file = os.path.join(rte_rrtmgp_dir, "rrtmgp-gas-lw-g256.nc")
sw_gas_coeffs_file = os.path.join(rte_rrtmgp_dir, "rrtmgp-gas-sw-g224.nc")

lw_clouds_coeff_file = os.path.join(rte_rrtmgp_dir, "rrtmgp-clouds-lw.nc")
sw_clouds_coeff_file = os.path.join(rte_rrtmgp_dir, "rrtmgp-clouds-sw.nc")

# In the local directory
all_sky_exe_name = os.path.join(all_sky_dir, "rrtmgp_allsky")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Runs all-sky examples, resetting output.")
    parser.add_argument("--run_command", type=str, default="",
                        help="Prefix ('jsrun' etc.) for running commands. "
                             "Use quote marks to enclose multi-part commands.")
    parser.add_argument("--ncol", type=int, default=24,
                        help="Number of columns to compute")
    parser.add_argument("--nlay", type=int, default=72,
                        help="Number of layers "
                             "(every one will have the same clouds)")
    parser.add_argument("--nloops", type=int, default=1,
                        help="Number of times to compute 'nloops' "
                             "cloudy columns")

    args = parser.parse_args()
    ncol_str   = '{0:5d}'.format(args.ncol)
    nlay_str   = '{0:5d}'.format(args.nlay)
    nloops_str = '{0:5d}'.format(args.nloops)
    if args.run_command:
        print("using the run command")
        all_sky_exe_name = args.run_command + " " + all_sky_exe_name
    os.chdir(all_sky_dir)
    # Remove cloudy-sky fluxes from the file containing the atmospheric profiles
    subprocess.run(
        [all_sky_exe_name, ncol_str, nlay_str, nloops_str, "rrtmgp-allsky-lw-no-aerosols.nc", lw_gas_coeffs_file, lw_clouds_coeff_file])
    subprocess.run(
        [all_sky_exe_name, ncol_str, nlay_str, nloops_str, "rrtmgp-allsky-sw-no-aerosols.nc", sw_gas_coeffs_file, sw_clouds_coeff_file])

# end main
