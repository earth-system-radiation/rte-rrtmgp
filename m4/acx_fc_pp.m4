# ACX_FC_PP_SRCEXT(EXTENSION,
#                  [ACTION-IF-SUCCESS],
#                  [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is known as
# AC_FC_PP_SRCEXT.
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to enable Fortran preprocessing for
# source files with extension EXTENSION. The result is either "unknown",
# or the actual compiler flag required to enable Fortran preprocessing, which
# may be an empty string.
#
# If successful, sets the output variable FCFLAGS_[]EXTENSION and the shell
# variable ac_fcflags_srcext to the result, sets the shell variables
# ac_fc_srcext and ac_ext to the EXTENSION, and runs ACTION-IF-SUCCESS;
# otherwise sets the output variable FCFLAGS_[]EXTENSION and the shell variable
# ac_fcflags_srcext to empty strings, keeps the shell variables ac_fc_srcext
# and ac_ext unchanged, and runs ACTION-IF-FAILURE (defaults to failing with an
# error message).
#
# The flag is cached in the acx_cv_fc_pp_srcext_[]EXTENSION, which may contain
# whitespaces.
#
# Known flags:
# gfortran: -cpp
# SGI: -ftpp
# SUN: -xpp={fpp,cpp}
# IBM: -qsuffix=cpp=EXTENSION
# HP: +cpp
# PGI: -Mpreprocess
# Absoft: -cpp
# Cray: -e T, -e Z
# Intel: -fpp (-Tf may also be needed right before the source file name)
# PathScale: -ftpp, -cpp
# Lahey: -Cpp
# NAGWare: -fpp
# Compaq/Tru64: -cpp
# f2c: -cpp
# g95: -cpp
#
AC_DEFUN([ACX_FC_PP_SRCEXT],
  [m4_provide([AC_FC_PP_SRCEXT])dnl
   AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_pp_srcext_$1])dnl
   AC_MSG_CHECKING([for Fortran compiler flag needed to compile dnl
preprocessed .$1 files])
   AC_CACHE_VAL([acx_cache_var],
     [acx_cache_var=unknown
      acx_ext_save=$ac_ext
      acx_fcflags_srcext_save=$ac_fcflags_srcext
      ac_ext=$1
      AS_CASE([$acx_fc_pp_srcext_ext],
        [[[fF]]77], [acx_fc_pp_srcext_try=f77-cpp-input],
        [acx_fc_pp_srcext_try=f95-cpp-input])
      for ac_fcflags_srcext in '' -ftpp -fpp -Tf '-fpp -Tf' -xpp=fpp dnl
-Mpreprocess '-e T' '-e Z' -cpp -xpp=cpp -qsuffix=cpp=$1 dnl
"-x $acx_fc_pp_srcext_try" +cpp -Cpp; do
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
[[#if 0
#include <ac_nonexistent.h>
      choke me
#endif]])],
          [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
[[#if 1
#include <ac_nonexistent.h>
      choke me
#endif]])],
             [],
             [acx_cache_var=$ac_fcflags_srcext])])
        test "x$acx_cache_var" != xunknown && break
      done
      ac_fcflags_srcext=$acx_fcflags_srcext_save
      ac_ext=$acx_ext_save])
   AS_IF([test -n "$acx_cache_var"],
     [AC_MSG_RESULT([$acx_cache_var])],
     [AC_MSG_RESULT([none needed])])
   AC_SUBST([FCFLAGS_][$1])
   AS_VAR_IF([acx_cache_var], [unknown],
     [ac_fcflags_srcext=
      AS_VAR_COPY([FCFLAGS_][$1], [ac_fcflags_srcext])
      m4_default([$3],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
compile preprocessed .$1 files])])],
     [ac_fcflags_srcext=$acx_cache_var
      AS_VAR_COPY([FCFLAGS_][$1], [ac_fcflags_srcext])
      ac_fc_srcext=$1
      AS_VAR_COPY([ac_ext], [ac_fc_srcext])
      $2])
   m4_popdef([acx_cache_var])])
