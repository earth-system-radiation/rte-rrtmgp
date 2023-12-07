#!/bin/sh

script_dir=`echo "$0" | sed 's@[^/]*$@@'`
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH
cd "$script_dir"

autoreconf -fvi || exit $?
