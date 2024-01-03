set -e
./optical_prop_unit_tests
./rte_unit_tests
./crash_reproducer
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
./check_variants    test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
./check_variants    test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc 