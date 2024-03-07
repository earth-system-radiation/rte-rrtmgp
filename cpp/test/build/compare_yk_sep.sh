#! /bin/bash -x

# Run test_lw and grep for sentinel string that report values
# for arrays/views you are having problems with. To use this
# feature, add calls to conv::pNd(...) for both yakl and kokkos
# impls for views/arrays that are not getting equivalent results.
# This is an easy way to quickly isolate problems.

# Remove previous temp files
/bin/rm JGF*

# Build
../cmakescript.sh -DRRTMGP_ENABLE_KOKKOS=On -DCMAKE_BUILD_TYPE=Debug

make -j16
make -j16
make -j16
make -j16
make -j16 || exit 1

# Run lw and grab kokkos data, use sed to make sure it matches yakl sentinel
# Not sure if we want to sort since that hides iteration order changes.
../test_sw.sh | grep JGFK | sort | sed 's/JGFK/JGFY/g' > JGFK

# Build
../cmakescript.sh -DRRTMGP_ENABLE_YAKL=On -DCMAKE_BUILD_TYPE=Debug

make -j16
make -j16
make -j16
make -j16
make -j16 || exit 1

# Run lw and grab yakl data
../test_sw.sh | grep JGFY | sort > JGFY

# Check and report diffs
diff JGFK JGFY
wc -l JGFK
wc -l JGFY
diff JGFK JGFY | wc -l
