# This file contains simplified (but probably less robust) versions of macros
# from https://github.com/skosukhin/mkhelper/tree/master/m4/acx_compiler.m4

# ACX_COMPILER_FC_VENDOR()
# -----------------------------------------------------------------------------
# Detects the vendor of the Fortran compiler. The result is "intel", "nag",
# "portland", "cray", "gnu" or "unknown".
#
# The result is cached in the acx_cv_fc_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_FC_VENDOR],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler vendor], [acx_cache_var],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^ifort (IFORT)' >/dev/null 2>&1],
        [acx_cache_var=intel],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
grep '^NAG Fortran Compiler Release' >/dev/null 2>&1],
        [acx_cache_var=nag],
        [AS_VAR_GET([_AC_CC]) -V 2>/dev/null| dnl
grep '^Copyright.*\(The Protland Group\|NVIDIA CORPORATION\)' >/dev/null 2>&1],
        [acx_cache_var=portland],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | grep '^Cray Fortran' >/dev/null 2>&1],
        [acx_cache_var=cray],
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^GNU Fortran' >/dev/null 2>&1],
        [acx_cache_var=gnu],
        [acx_cache_var=unknown])])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_FC_VERSION()
# -----------------------------------------------------------------------------
# Detects the version of the Fortran compiler. The result is either "unknown"
# or a string in the form "epoch:major.minor.patchversion", where "epoch:" is
# an optional prefix used in order to have an increasing version number in case
# of marketing change.
#
# The result is cached in the acx_cv_fc_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_FC_VERSION],
  [AC_REQUIRE([ACX_COMPILER_FC_VENDOR])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_version])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler version], [acx_cache_var],
     [AS_CASE([AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])],
        [intel],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^ifort (IFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p'`]],
        [nag],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/^NAG Fortran Compiler Release \([0-9][0-9]*\.[0-9][0-9]*\).*]dnl
[Build \([0-9][0-9]*\)/\1.\2/p'`]],
        [portland],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V | dnl
[sed -n 's/pgfortran \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\1.\2/p'`]],
        [cray],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*ersion \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p'`]],
        [gnu],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpversion 2>/dev/null`],
        [acx_cache_var=unknown])])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_CC_VENDOR()
# -----------------------------------------------------------------------------
# Detects the vendor of the C compiler. The result is  "intel", "nag",
# "portland", "cray", "gnu" or "unknown".
#
# The result is cached in the acx_cv_c_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_CC_VENDOR],
  [AC_LANG_ASSERT([C])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler vendor], [acx_cache_var],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^icc (ICC)' >/dev/null 2>&1],
        [acx_cache_var=intel],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
grep '^NAG Fortran Compiler Release' >/dev/null 2>&1],
        [acx_cache_var=nag],
        [AS_VAR_GET([_AC_CC]) -V 2>/dev/null| dnl
grep '^Copyright.*\(The Protland Group\|NVIDIA CORPORATION\)' >/dev/null 2>&1],
        [acx_cache_var=portland],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | grep '^Cray C' >/dev/null 2>&1],
        [acx_cache_var=cray],
        [AS_VAR_GET([_AC_CC]) --version | grep '^gcc' >/dev/null 2>&1],
        [acx_cache_var=gnu],
        [acx_cache_var=unknown])])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_CC_VERSION()
# -----------------------------------------------------------------------------
# Detects the version of the C compiler. The result is either "unknown"
# or a string in the form "epoch:major.minor.patchversion", where "epoch:" is
# an optional prefix used in order to have an increasing version number in case
# of marketing change.
#
# The result is cached in the acx_cv_c_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_CC_VERSION],
  [AC_REQUIRE([ACX_COMPILER_CC_VENDOR])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_version])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler version], [acx_cache_var],
     [AS_CASE([AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])],
        [intel],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^icc (ICC) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p'`]],
        [nag],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/^NAG Fortran Compiler Release \([0-9][0-9]*\.[0-9][0-9]*\).*]dnl
[Build \([0-9][0-9]*\)/\1.\2/p'`]],
        [portland],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V | dnl
[sed -n 's/pgcc \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\1.\2/p'`]],
        [cray],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*ersion \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p'`]],
        [gnu],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpversion 2>/dev/null`],
        [acx_cache_var=unknown])])
   m4_popdef([acx_cache_var])])
