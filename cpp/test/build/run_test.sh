#!/bin/bash -ex

printf "Running allsky longwave tests\n\n"
./test_lw.sh

printf "Running allsky shortwave tests\n\n"
./test_sw.sh

