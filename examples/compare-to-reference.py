#! /usr/bin/env python
#
# General purpose comparison script -- compare all variables in a set of files,
#   write output if differences exceed some threshold,
#   error exit if differences exceed a different threshold
#
# Currently thresholds are specified as absolute differences
#    worth revising but could change development practice
# Thresholds come from environement variables when set?
#
#
import argparse
import os
import sys
import warnings

import numpy as np
import xarray as xr

#
# Comparing reference and test results
#
if __name__ == "__main__":
    warnings.simplefilter("ignore", xr.SerializationWarning)

    parser = argparse.ArgumentParser(
        description="Compares example output to file in reference directory"
    )
    parser.add_argument(
        "--ref_dir", type=str, help="Directory containing reference files"
    )
    parser.add_argument(
        "--tst_dir", type=str, default=".", help="Directory contining test values"
    )
    parser.add_argument(
        "--file_names",
        type=str,
        nargs="+",
        default=[],
        help="Name[s] of files to compare",
    )
    parser.add_argument(
        "--variables",
        type=str,
        nargs="+",
        default=None,
        help="Name[s] of files to compare",
    )
    parser.add_argument(
        "--report_threshold",
        type=float,
        default=os.getenv("REPORTING_THRESHOLD", 0.0),
        help="Threshold for reporting differences",
    )
    parser.add_argument(
        "--failure_threshold",
        type=float,
        default=os.getenv("FAILURE_THRESHOLD", 1.0e-5),
        help="Threshold at which differences cause failure "
        "(for continuous integration)",
    )
    args = parser.parse_args()

    tst = xr.open_mfdataset(
        [os.path.join(args.tst_dir, f) for f in args.file_names], combine="by_coords"
    )
    ref = xr.open_mfdataset(
        [os.path.join(args.ref_dir, f) for f in args.file_names], combine="by_coords"
    )
    variables = args.variables if args.variables is not None else tst.variables

    failed = False
    for v in variables:
        if np.any(np.isnan(ref.variables[v].values)):
            raise Exception(v + ": some ref values are missing. That's not right.")
        if np.all(np.isnan(tst.variables[v].values)):
            raise Exception(v + ": all test values are missing. Were the tests run?")
        if np.any(np.isnan(tst.variables[v].values)):
            raise Exception(v + ": some test values are missing. Now that is strange.")
        #
        # Reporting
        #
        if not np.allclose(tst[v], ref[v], atol=args.report_threshold, rtol=0):
            diff = abs((tst - ref)[v].values)
            avg = 0.5 * (tst + ref)[v].values
            # Division raises a runtime warning when we divide by zero even if
            # the values in those locations will be ignored.
            with np.errstate(divide="ignore", invalid="ignore"):
                frac_diff = np.where((avg > 2.0 * np.finfo(float).eps), diff / avg, 0)
                print(
                    "Variable %s differs (max abs difference: %e; "
                    "max percent difference: %e%%)"
                    % (v, diff.max(), 100.0 * frac_diff.max())
                )
        else:
            print("Variable %s: No diffs" % v)
        #
        # Failure
        #
        if not np.allclose(tst[v], ref[v], atol=args.failure_threshold, rtol=0):
            failed = True

    if failed:
        print("Tests failed")
    sys.exit(1) if failed else sys.exit(0)
