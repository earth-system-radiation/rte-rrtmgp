#! /usr/bin/env python
#
# This script compares RFMIP results from RTE+RRTMGP against a benchmark
#
import os
import numpy as np
import xarray as xr
import urllib.request

ref_dir  = "reference"
tst_dir = "."

#
# Download reference results
#
print("Downloading reference results")
ref_dir = "./reference/"
if not os.path.exists(ref_dir):
    os.makedirs(ref_dir)
urllib.request.urlretrieve("https://owncloud.gwdg.de/index.php/s/kbhl3JOSccGtR0m/download", \
                           os.path.join(ref_dir, "rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
urllib.request.urlretrieve("https://owncloud.gwdg.de/index.php/s/5DbhryVSfztioPG/download", \
                           os.path.join(ref_dir, "rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
urllib.request.urlretrieve("https://owncloud.gwdg.de/index.php/s/uCemCHlGxbGK0gJ/download", \
                           os.path.join(ref_dir, "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))
urllib.request.urlretrieve("https://owncloud.gwdg.de/index.php/s/l8ZG28j9ttZWD9r/download", \
                           os.path.join(ref_dir, "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"))

#
# Comparing reference and test results
#
rrtmg_suffix = "_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"
tst = xr.open_mfdataset(os.path.join(tst_dir, "r??" + rrtmg_suffix))
ref = xr.open_mfdataset(os.path.join(ref_dir, "r??" + rrtmg_suffix))

for v in ['rlu', 'rld', 'rsu', 'rsd']:
  if np.all(np.isnan(tst.variables[v].values)):
    raise Exception("All test values are missing. Were the tests run?")
  if np.any(np.isnan(tst.variables[v].values)):
    raise Exception("Some test values are missing. Now that is strange.")

  diff = abs((tst-ref).variables[v].values)
  avg  = 0.5*(tst+ref).variables[v].values
  # Division raises a runtime warning when we divide by zero even if the
  #   values in those locations will be ignored.
  with np.errstate(divide='ignore', invalid='ignore'):
    frac_diff = np.where((avg > 2.*np.finfo(float).eps), diff/avg, 0)

  if diff.max() > 0:
    print('Variable %s differs (max abs difference: %e; max frac. difference: %e%%)'%(v, diff.max(), 100.0 * frac_diff.max()))
  else:
    print('Variable %s: No diffs'%(v))
