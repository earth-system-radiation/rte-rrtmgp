# ACX_LANG_OPENMP_FLAG([ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is a part of
# AC_OPENMP.
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to enable OpenMP support. The result is either
# "unknown", or the actual compiler flag required to enable OpenMP support,
# which may be an empty string.
#
# Known flags:
# Intel >= 16: -qopenmp
# NAG, Intel < 16: -openmp
# GNU: -fopenmp
# Cray: -homp
# SGI, PGI: -mp
# SunPRO: -xopenmp
# Tru64 Compaq C: -omp
# IBM XL: -qsmp=omp
# NEC SX: -Popenmp
# Lahey Fortran: --openmp
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_openmp_flag variable.
#
# Upon successful run, you can check for the version of the standard supported
# by the compiler by expanding:
#   AS_VAR_APPEND([_AC_LANG_PREFIX[]FLAGS],
#     [" $acx_cv_[]_AC_LANG_ABBREV[]_openmp_flag"])
#   ACX_LANG_MACRO_CHECK_VALUE([_OPENMP])
# and checking for the value of the
# acx_cv_[]_AC_LANG_ABBREV[]_macro__OPENMP_value shell variable. The possible
# (successful) values of the variable are dates, which map to the versions of
# the standard in the following way:
#   202011 5.1
#   201811 5.0
#   201307 4.0
#   201107 3.1
#   200805 3.0
#   200505 2.5
#   200203 C/C++ version 2.0
#   200011 Fortran version 2.0
#   199911 Fortran version 1.1
#   199810 C/C++ version 1.0
#   199710 Fortran version 1.0
#
AC_DEFUN([ACX_LANG_OPENMP_FLAG],
  [m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_openmp_flag])dnl
   AC_MSG_CHECKING([for _AC_LANG compiler flag needed to enable OpenMP dnl
support])
   AC_CACHE_VAL([acx_cache_var],
     [acx_cache_var=unknown
      acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      AC_LANG_CONFTEST([_ACX_LANG_OPENMP])
      for acx_lang_openmp_flag in '' -qopenmp -openmp -fopenmp -homp -mp dnl
-xopenmp -omp -qsmp=omp -Popenmp --openmp; do
        _AC_LANG_PREFIX[]FLAGS="${acx_save_[]_AC_LANG_PREFIX[]FLAGS} dnl
$acx_lang_openmp_flag"
        AC_LINK_IFELSE([], [acx_cache_var=$acx_lang_openmp_flag])
        test "x$acx_cache_var" != xunknown && break
      done
      rm -f conftest.$ac_ext
      _AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS])
   AS_IF([test -n "$acx_cache_var"],
     [AC_MSG_RESULT([$acx_cache_var])],
     [AC_MSG_RESULT([none needed])])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect _AC_LANG compiler flag needed to dnl
enable OpenMP support])])], [$1])
   m4_popdef([acx_cache_var])])

# _ACX_LANG_OPENMP()
# -----------------------------------------------------------------------------
# Expands into the source code of a program in the current language that is
# compiled successfully only when OpenMP support is enabled for the current
# compiler. By default, expands to _AC_LANG_OPENMP.
#
m4_define([_ACX_LANG_OPENMP],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [_AC_LANG_OPENMP])])])

# _ACX_LANG_OPENMP(Fortran)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_OPENMP for Fortran language. In addition to the
# standard implementation of _AC_LANG_OPENMP(Fortran), also checks whether the
# macro _OPENMP is set by the compiler (if AC_FC_PP_SRCEXT was expanded
# before).
#
m4_define([_ACX_LANG_OPENMP(Fortran)],
  [AC_LANG_SOURCE([AC_PROVIDE_IFELSE([AC_FC_PP_SRCEXT],
[m4_n([[#ifndef _OPENMP
      choke me
#endif]])])_AC_LANG_OPENMP])])
