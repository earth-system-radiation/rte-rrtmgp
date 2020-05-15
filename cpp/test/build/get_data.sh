#!/bin/bash

rm -f *.nc

wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc 
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rrtmgp-data-lw-g256-2018-12-04.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rrtmgp-data-sw-g224-2018-12-04.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rrtmgp-cloud-optics-coeffs-lw.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rrtmgp-cloud-optics-coeffs-sw.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rrtmgp-allsky.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
wget https://code.ornl.gov/imn/data/raw/master/rrtmgp/rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
 
