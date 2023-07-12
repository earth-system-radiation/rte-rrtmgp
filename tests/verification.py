import sys
import os
import warnings
import argparse

import xarray as xr
import numpy  as np

REPORT  = os.getenv("REPORTING_THRESHOLD") if os.getenv("REPORTING_THRESHOLD") is not None else 1.e-10
FAILURE = os.getenv("FAILURE_THRESHOLD")   if os.getenv("FAILURE_THRESHOLD")   is not None else 1.e-5


def assert_equal(variants, reference):
    #
    # Computes difference between reference and each variant, reports if
    # differences exceed a threshold
    #
    if type(variants) is not list:
        return assert_equal([variants], reference)
    else:
        passed = True
        print('Using %s as reference:' % reference.description)
        for v in variants:
            diff = np.fabs(v - reference)
            print('    %s' % v.description)
            if diff.max() > report_threshold:
                print('      differs from reference by as much as %e' % (
                    diff.max()))
            passed = passed and diff.max() <= failure_threshold
        print('')
        return passed


########################################################################
if __name__ == '__main__':
    warnings.simplefilter("ignore", xr.SerializationWarning)
 
    parser = argparse.ArgumentParser(
        description="Compares all-sky example output to file in reference "
                    "directory")
    parser.add_argument("--report_threshold", type=float, default=REPORT,
                        help="Threshold for reporting differences")
    parser.add_argument("--failure_threshold", type=float, default=FAILURE,
                        help="Threshold at which differences cause failure")
    args = parser.parse_args()
    report_threshold = args.report_threshold
    failure_threshold = args.failure_threshold

    passed = True
    gp = xr.open_dataset("test_atmospheres.nc")
    ########################################################################
    #
    # Some variants we expect to be essentially equal to the default,
    # so we'll check for those
    ###############################
    #
    # Longwave
    #
    #
    # Does the flux plus the Jacobian equal a calculation with perturbed surface
    # temperature?
    #
    gp['lw_flux_up_from_deriv'] = gp.lw_flux_up_jaco + gp.lw_jaco_up
    gp.lw_flux_up_from_deriv.attrs = {
        "description": "LW flux up, surface T+1K, computed from Jacobian"}
    passed = (passed and
              assert_equal(gp.lw_flux_up_from_deriv, gp.lw_flux_up_stp1))

 
    if not passed:
        print("Something is terribly, terribly wrong")
    sys.exit(0) if passed else sys.exit(1)
