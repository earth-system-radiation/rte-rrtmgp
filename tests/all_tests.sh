./rte_unit_tests
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc 
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc 
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc 
./check_equivalence test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc 
./check_variants    test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
./check_variants    test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
