# ACX_LANG_OPENACC_FLAG([ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to enable OpenACC support. The result is
# either "unknown", or the actual compiler flag required to enable OpenACC
# support, which may be an empty string.
#
# Known flags:
# GNU: -fopenacc
# Cray: -hacc
# PGI: -acc
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_openacc_flag variable.
#
# Upon successful run, you can check for the version of the standard supported
# by the compiler by expanding:
#   AS_VAR_APPEND([_AC_LANG_PREFIX[]FLAGS],
#     [" $acx_cv_[]_AC_LANG_ABBREV[]_openacc_flag"])
#   ACX_LANG_MACRO_CHECK_VALUE([_OPENACC])
# and checking for the value of the
# acx_cv_[]_AC_LANG_ABBREV[]_macro__OPENACC_value shell variable. The possible
# (successful) values of the variable are dates, which map to the versions of
# the standard in the following way:
#   201711 2.6
#   201510 2.5
#   201308 2.0 (the corrected version)
#   201306 2.0
#   201111 1.0
#
AC_DEFUN([ACX_LANG_OPENACC_FLAG],
  [m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_openacc_flag])dnl
   AC_MSG_CHECKING([for _AC_LANG compiler flag needed to enable OpenACC dnl
support])
   AS_VAR_SET_IF([acx_cache_var],
     [AS_ECHO_N(["(cached) "]) >&AS_MESSAGE_FD],
     [acx_cache_var=unknown
      acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      AC_LANG_CONFTEST([_ACX_LANG_OPENACC])
      for acx_lang_openacc_flag in '' -fopenacc -hacc -acc; do
        _AC_LANG_PREFIX[]FLAGS="${acx_save_[]_AC_LANG_PREFIX[]FLAGS} dnl
$acx_lang_openacc_flag"
        AC_LINK_IFELSE([],
          [acx_cache_var=$acx_lang_openacc_flag
           break])
      done
      rm -f conftest.$ac_ext
      _AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS])
   AS_IF([test -n "$acx_cache_var"],
     [AC_MSG_RESULT([$acx_cache_var])],
     [AC_MSG_RESULT([none needed])])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect _AC_LANG compiler flag needed to dnl
enable OpenACC support])])], [$1])
   m4_pushdef([acx_cache_var])])

# _ACX_LANG_OPENACC()
# -----------------------------------------------------------------------------
# Expands into the source code of a program in the current language that is
# compiled successfully only when OpenACC support is enabled for the current
# compiler. By default, expands to m4_fatal with the message saying that
# _AC_LANG is not supported.
#
m4_define([_ACX_LANG_OPENACC],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_fatal([the OpenACC test program is not defined for ]dnl
_AC_LANG[ language])])])])

# _ACX_LANG_OPENACC(C)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_OPENACC for C language.
#
m4_define([_ACX_LANG_OPENACC(C)],
  [AC_LANG_PROGRAM([[#include <openacc.h>]],
[[#ifndef _OPENACC
 choke me
#endif
return (int) acc_get_device_type()]])])

# _ACX_LANG_OPENACC(C++)
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_OPENACC for C++ language.
#
m4_copy([_ACX_LANG_OPENACC(C)], [_ACX_LANG_OPENACC(C++)])

# _ACX_LANG_OPENACC(Fortran)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_OPENACC for Fortran language.
#
# The implementation implies that the Fortran module should be available when
# OpenACC is enabled.
#
m4_define([_ACX_LANG_OPENACC(Fortran)],
  [AC_PROVIDE_IFELSE([AC_FC_PP_SRCEXT], [],
     [m4_warn([syntax],
        [ACX_LANG_OPENACC_FLAG requires calling the Fortran compiler with ]dnl
[a preprocessor but no call to AC_FC_PP_SRCEXT is detected])])dnl
   AC_LANG_PROGRAM([],
[[#ifndef _OPENACC
      choke me
#endif
      use openacc, only : acc_get_device_type
      implicit none
      integer :: t
      t = acc_get_device_type()]])])
