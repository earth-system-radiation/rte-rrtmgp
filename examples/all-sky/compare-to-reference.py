#! /usr/bin/env python
#
# This script compares RFMIP results from RTE+RRTMGP against a benchmark
#
import argparse
import os
import sys

import numpy as np
import xarray as xr

#
# Comparing reference and test results
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compares all-sky example output to file in reference "
                    "directory")
    parser.add_argument("--allsky_file", type=str, default="rrtmgp-allsky-lw.nc",
                        dest="file",
                        help="Name of file inputs and outputs for "
                             "all-sky problem (same for test and reference)")
    parser.add_argument("--ref_dir", type=str, default=os.path.join(os.environ["RRTMGP_DATA"], "examples", "all-sky", "reference"),
                        help="Directory where reference values are")
    parser.add_argument("--report_threshold", type=float, default=0.,
                        help="Threshold for reporting differences")
    parser.add_argument("--failure_threshold", type=float, default=1.e-5,
                        help="Threshold at which differences cause failure "
                             "(for continuous integration)")
    args = parser.parse_args()

    tst_file = args.file
    ref_file = os.path.join(args.ref_dir, tst_file)
    tst = xr.open_dataset(tst_file)
    ref = xr.open_dataset(ref_file)

    failed = False
    for v in tst.variables:
        if np.any(np.isnan(ref.variables[v].values)):
            raise Exception(v + ": some ref values are missing. Now that is strange.")
        if np.all(np.isnan(tst.variables[v].values)):
            raise Exception("All test values are missing. Were the tests run?")
        if np.any(np.isnan(tst.variables[v].values)):
            raise Exception(v + ":Some test values are missing. Now that is strange.")

        # express as (tst-ref).variables[v].values when replacing reference file
        # to have same number of columns
        diff = abs((tst.variables[v] - ref.variables[v]).values)
        avg = 0.5 * (tst.variables[v] + ref.variables[v]).values
        # Division raises a runtime warning when we divide by zero even if the
        # values in those locations will be ignored.
        with np.errstate(divide='ignore', invalid='ignore'):
            frac_diff = np.abs(
                np.where((avg > 2. * np.finfo(float).eps), diff / avg, 0))

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
