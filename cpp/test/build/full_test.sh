#!/bin/bash -ex

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

# Test both
$this_dir/cmakescript.sh -DRRTMGP_ENABLE_YAKL=On -DRRTMGP_ENABLE_KOKKOS=On -DCMAKE_BUILD_TYPE=Debug
make
$this_dir/test_lw.sh
$this_dir/test_sw.sh

# Just YAKL
$this_dir/cmakescript.sh -DRRTMGP_ENABLE_YAKL=On -DCMAKE_BUILD_TYPE=Debug
make
$this_dir/test_lw.sh
$this_dir/test_sw.sh

# Just Kokkos
$this_dir/cmakescript.sh -DRRTMGP_ENABLE_KOKKOS=On -DCMAKE_BUILD_TYPE=Debug
make
$this_dir/test_lw.sh
$this_dir/test_sw.sh
