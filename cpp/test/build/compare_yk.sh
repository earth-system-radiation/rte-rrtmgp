#! /bin/bash -x

# Run test_lw and grep for sentinel string that report values
# for arrays/views you are having problems with. To use this
# feature, add calls to conv::pNd(...) for both yakl and kokkos
# impls for views/arrays that are not getting equivalent results.
# This is an easy way to quickly isolate problems.

# Build
make -j8 || exit 1

# Remove previous temp files
/bin/rm JGF*

# Run lw and grab kokkos data, use sed to make sure it matches yakl sentinel
# Not sure if we want to sort since that hides iteration order changes.
../test_sw.sh | grep JGFK | sort | sed 's/JGFK/JGFY/g' > JGFK

# Run lw and grab yakl data
../test_sw.sh | grep JGFY | sort > JGFY

# Check and report diffs
diff JGFK JGFY
wc -l JGFK
wc -l JGFY
diff JGFK JGFY | wc -l
