set -eux
if [[ $(uname) == *"NT"* ]]; then
    # Rename files ending with .exe to the same name without the extension
    for file in *.exe; do
        mv "$file" "${file%.exe}"
    done
fi
./rte_optic_prop_unit_tests
./rte_lw_solver_unit_tests
./rte_sw_solver_unit_tests
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc 
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
./check_variants    test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
./check_variants    test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc 