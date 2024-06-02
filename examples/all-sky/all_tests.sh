set -eux
if [[ $(uname) == *"NT"* ]]; then
    # Rename files ending with .exe to the same name without the extension
    for file in *.exe; do
        mv "$file" "${file%.exe}"
    done
fi
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc 
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc 
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc 
./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc  