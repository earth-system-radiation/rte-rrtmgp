#! /usr/bin/env python
#
# This script compares RFMIP results from RTE+RRTMGP against a benchmark
#
import os
import numpy as np
import xarray as xr

ref_dir  = "reference"
tst_dir = "."

rrtmg_suffix = "_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"

tst = xr.open_mfdataset(os.path.join(tst_dir, "r??" + rrtmg_suffix))
ref = xr.open_mfdataset(os.path.join(ref_dir, "r??" + rrtmg_suffix))

for v in ['rlu', 'rld', 'rsu', 'rsd']:
  if np.all(np.isnan(tst.variables[v].values)):
    raise Exception("All test values are missing. Were the tests run?")
  diff = abs((tst-ref).variables[v].values)
  avg  = 0.5 * (tst+ref).variables[v].values
  if np.any((diff > 0) & (avg > 0)):
    frac_dif = np.where((avg > 0) & (diff > 0), diff/avg, 0)
  else:
    frac_dif = np.zeros(1)

  if diff.max() > 0:
    print('Variable %s differs (max abs difference: %e; max frac. difference(%): %e%%)'% \
          (v, diff.max(), 100.0 * frac_diff.max()))
  else:
    print('Variable %s: No diffs'%(v))
