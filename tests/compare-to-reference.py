#! /usr/bin/env python
#
# This script compares RFMIP results from RTE+RRTMGP against a benchmark
#
import os
import sys
import numpy as np
import xarray as xr
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compares cloud scattering example output to reference")
    parser.add_argument("--tst_file", type=str, default="cloud_scattering.nc",
                        help="File name for reference values")
    parser.add_argument("--ref_file", type=str, default="cloud_scattering_ref.nc",
                        help="File name for reference values")
    parser.add_argument("--report_threshold", type=float, default=0.,
                        help="Threshold for reporting differences")
    parser.add_argument("--failure_threshold", type=float, default=1.e-5,
                        help="Threshold at which differences cause failure (for continuous integration)")
    args = parser.parse_args()

    ref = xr.open_dataset(args.ref_file)
    tst = xr.open_dataset(args.tst_file)
    vars = [v for v in ref.variables if "up" in v or "dn" in v]

    failed = False
    for v in vars:
      if np.all(np.isnan(tst[v].values)):
        raise Exception(var + ": all test values are missing. Were the tests run?")
      if np.any(np.isnan(tst[v].values)):
        raise Exception(var + ": some test values are missing. Now that is strange.")

      diff = abs((tst-ref)[v].values)
      avg  = 0.5*(tst+ref)[v].values
      # Division raises a runtime warning when we divide by zero even if the
      #   values in those locations will be ignored.
      with np.errstate(divide='ignore', invalid='ignore'):
        frac_diff = np.where((avg > 2.*np.finfo(float).eps), diff/avg, 0)

        if diff.max() > args.report_threshold:
            print('Variable %s differs (max abs difference: %e; max percent difference: %e%%)'%(v, diff.max(), 100.0 * frac_diff.max()))
        else:
            print('Variable %s: No diffs'%(v))

        if diff.max() > args.failure_threshold: failed = True

    sys.exit(1) if failed else sys.exit(0)
