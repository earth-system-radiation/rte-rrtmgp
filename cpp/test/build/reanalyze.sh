#!/bin/bash -ex

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

# rebuild
make -j20

# remove prior files
/bin/rm -f *JGF*

OMP_NUM_THREADS=64 $this_dir/test_lw_perf.sh 100 | grep TIMING > TIME_DATA_JGF

$this_dir/rrtmgp-perf-analysis TIME_DATA_JGF > ANALYSIS_JGF

$this_dir/rrtmgp-perf-analysis ANALYSIS_JGF KANAL_ORIG
