set -eux
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw-bnd.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc 
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw-bnd.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc 
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw-bnd.nc 
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw-bnd.nc  
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols-gpt-clouds.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw-g256.nc  
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols-gpt-clouds.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw-g224.nc  
