#! /usr/bin/env python
#
# This script compares RFMIP results from RTE+RRTMGP against a benchmark
#
import argparse
import os
import sys
import urllib.request

import numpy as np
import xarray as xr

tst_dir = "."
rrtmgp_suffix = "_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"


def construct_esgf_remote_name(var):
    #
    # For a given variable name, provide the OpenDAP URL for the RTE+RRTMGP
    # RFMIP results
    #   This doesn't seem to work on CSCS Piz Daint within the netcdf-python
    #   module
    #
    esgf_url_base = ("http://esgf3.dkrz.de/thredds/dodsC/"
                     "cmip6/RFMIP/RTE-RRTMGP-Consortium/RTE-RRTMGP-181204/"
                     "rad-irf/r1i1p1f1/Efx/")
    # DKRZ server has been unstable - better to try the other if one fails
    esgf_url_base = ("http://esgf-data1.llnl.gov/thredds/dodsC/css03_data/"
                     "CMIP6/RFMIP/RTE-RRTMGP-Consortium/RTE-RRTMGP-181204/"
                     "rad-irf/r1i1p1f1/Efx/")
    esgf_url_ver = "gn/v20191007/"
    return os.path.join(esgf_url_base, var, esgf_url_ver, var + rrtmgp_suffix)


#
# Construct URL for RTE+RRTMGP results for RFMIP from ESGF
#
def construct_esgf_file(var):
    esgf_url_base = ("http://esgf3.dkrz.de/thredds/fileServer/"
                     "cmip6/RFMIP/RTE-RRTMGP-Consortium/RTE-RRTMGP-181204/"
                     "rad-irf/r1i1p1f1/Efx/")
    # DKRZ node goes down frequently
    esgf_url_base = ("http://esgf-data1.llnl.gov/thredds/fileServer/css03_data/"
                     "CMIP6/RFMIP/RTE-RRTMGP-Consortium/RTE-RRTMGP-181204/"
                     "rad-irf/r1i1p1f1/Efx/")
    esgf_url_ver = "gn/v20191007/"
    return os.path.join(esgf_url_base, var, esgf_url_ver, var + rrtmgp_suffix)


#
# Comparing reference and test results
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compares all-sky example output to file in reference "
                    "directory")
    parser.add_argument("--ref_dir", type=str, default="reference",
                        help="Directory where reference values are")
    parser.add_argument("--download_reference", action="store_true",
                        help="Download reference files even if they exist")
    parser.add_argument("--report_threshold", type=float, default=0.,
                        help="Threshold for reporting differences")
    parser.add_argument("--failure_threshold", type=float, default=1.e-5,
                        help="Threshold at which differences cause failure "
                             "(for continuous integration)")
    args = parser.parse_args()

    variables = ['rlu', 'rld', 'rsu', 'rsd']
    # Download reference data
    #   If versions of all files exist in the reference directory, no need to
    #   download (can be over-ridden)
    if not all(
            [os.path.exists(os.path.join(args.ref_dir, v + rrtmgp_suffix))
             for v in variables]) or args.download_reference:
        print("Downloading reference data")
        os.makedirs(args.ref_dir, exist_ok=True)
        for v in variables:
            urllib.request.urlretrieve(construct_esgf_file(v),
                                       os.path.join(args.ref_dir,
                                                    v + rrtmgp_suffix))

    tst = xr.open_mfdataset(os.path.join(tst_dir, "r??" + rrtmgp_suffix),
                            combine='by_coords')
    ref = xr.open_mfdataset(os.path.join(args.ref_dir, "r??" + rrtmgp_suffix),
                            combine='by_coords')

    failed = False
    for v in variables:
        if np.all(np.isnan(tst.variables[v].values)):
            raise Exception(
                v + ": all test values are missing. Were the tests run?")
        if np.any(np.isnan(tst.variables[v].values)):
            raise Exception(
                v + ": some test values are missing. Now that is strange.")

        diff = abs((tst - ref).variables[v].values)
        avg = 0.5 * (tst + ref).variables[v].values
        # Division raises a runtime warning when we divide by zero even if the
        # values in those locations will be ignored.
        with np.errstate(divide='ignore', invalid='ignore'):
            frac_diff = np.where(
                (avg > 2. * np.finfo(float).eps), diff / avg, 0)

            if diff.max() > args.report_threshold:
                print(
                    'Variable %s differs (max abs difference: %e; '
                    'max percent difference: %e%%)' % (
                        v, diff.max(), 100.0 * frac_diff.max()))
            else:
                print('Variable %s: No diffs' % v)

            if diff.max() > args.failure_threshold:
                failed = True

    sys.exit(1) if failed else sys.exit(0)
