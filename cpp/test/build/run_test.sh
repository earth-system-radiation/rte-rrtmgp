#!/bin/bash -ex

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

printf "Running allsky longwave tests\n\n"
$this_dir/test_lw.sh

printf "Running allsky shortwave tests\n\n"
$this_dir/test_sw.sh

