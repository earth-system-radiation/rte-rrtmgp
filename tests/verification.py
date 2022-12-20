import argparse
import sys

import xarray as xr
import numpy  as np


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
    parser = argparse.ArgumentParser(
        description="Compares all-sky example output to file in reference "
                    "directory")
    parser.add_argument("--report_threshold", type=float, default=1.e-10,
                        help="Threshold for reporting differences")
    parser.add_argument("--failure_threshold", type=float, default=1.e-5,
                        help="Threshold at which differences cause failure")
    args = parser.parse_args()
    report_threshold = args.report_threshold
    failure_threshold = args.failure_threshold

    gp = xr.open_dataset("test_atmospheres.nc")
    ########################################################################
    #
    # Some variants we expect to be essentially equal to the default,
    # so we'll check for those
    ###############################
    #
    # Longwave
    #
    gp['lw_flux_net_from_updn'] = gp.lw_flux_dn - gp.lw_flux_up
    gp.lw_flux_net_from_updn.attrs = {
        "description": "LW flux net, computed externally from dn-up"}

    passed = assert_equal(
        [gp.lw_flux_dn_vr, gp.lw_flux_dn_jaco, gp.lw_flux_dn_subset],
        gp.lw_flux_dn)
    passed = (passed and
              assert_equal(
                  [gp.lw_flux_up_vr, gp.lw_flux_up_jaco, gp.lw_flux_up_subset],
                  gp.lw_flux_up))
    passed = (passed and
              assert_equal([gp.lw_flux_net, gp.lw_flux_net_2],
                           gp.lw_flux_net_from_updn))
    #
    # Does the flux plus the Jacobian equal a calculation with perturbed surface
    # temperature?
    #
    gp['lw_flux_up_from_deriv'] = gp.lw_flux_up_jaco + gp.lw_jaco_up
    gp.lw_flux_up_from_deriv.attrs = {
        "description": "LW flux up, surface T+1K, computed from Jacobian"}
    passed = (passed and
              assert_equal(gp.lw_flux_up_from_deriv, gp.lw_flux_up_stp1))
    ###############################
    #
    # Shortwave
    #
    passed = (passed and
              assert_equal([gp.sw_flux_dn_vr, gp.sw_flux_dn_tsi],
                           gp.sw_flux_dn))
    passed = (passed and
              assert_equal([gp.sw_flux_up_vr, gp.sw_flux_up_tsi],
                           gp.sw_flux_up))

    ###############################
    #
    # Incrementing
    #
    print('Incrementing')
    passed = assert_equal(
        [gp.lw_flux_dn_inc_1scl_with_1scl, gp.lw_flux_dn_inc_1scl_with_2str,
         gp.lw_flux_dn_inc_1scl_with_nstr],
        gp.lw_flux_dn) and passed
    passed = assert_equal(
        [gp.lw_flux_up_inc_1scl_with_1scl, gp.lw_flux_up_inc_1scl_with_2str,
         gp.lw_flux_up_inc_1scl_with_nstr],
        gp.lw_flux_up) and passed
    # passed = assert_equal(gp.lw_flux_dn_inc_2str_with_1scl,
    #                       gp.lw_flux_dn_2str) and passed
    # passed = assert_equal(gp.lw_flux_up_inc_2str_with_1scl,
    #                       gp.lw_flux_up_2str) and passed
    passed = assert_equal([gp.sw_flux_dn_incr], gp.sw_flux_dn) and passed
    passed = assert_equal([gp.sw_flux_up_incr], gp.sw_flux_up) and passed

    if not passed:
        print("Something is terribly, terribly wrong")
    sys.exit(0) if passed else sys.exit(1)
