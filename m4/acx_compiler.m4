# If you need to add support for a new vendor, add it to the language-specific
# version of the list _ACX_COMPILER_KNOWN_VENDORS. If the vendor sets
# additional preprocessor macros, put them in the list too, otherwise, you will
# have to implement additional detection algorithm (see the case of NAG in
# ACX_COMPILER_CC_VENDOR). All the language-specific version detection macros
# (ACX_COMPILER_FC_VERSION, ACX_COMPILER_CC_VERSION, etc.) expand the same
# macro _ACX_COMPILER_VERSION, which goes through <vendor>'s in the relevant
# version of _ACX_COMPILER_KNOWN_VENDORS and checks whether
# _ACX_COMPILER_VERSION_<VENDOR> or
# _ACX_COMPILER_VERSION_<VENDOR>(<current language>) is defined. If either of
# them (the language-specific version has higher priority) is defined,
# _ACX_COMPILER_VERSION will expand it. For example, Tiny CC (has vendor ID
# 'tcc') is only a C compiler and only _ACX_COMPILER_VERSION_TCC(C) is defined,
# which means that this option is added to the configure script only if the
# current language is set to C. At the same time, the macro for NAG compiler -
# _ACX_COMPILER_VERSION_NAG - does not have a language-specific version, which
# means, given that 'nag' is listed in _ACX_COMPILER_KNOWN_VENDORS(Fortran) and
# in _ACX_COMPILER_KNOWN_VENDORS(C), the version of NAG is checked the same way
# in the case of Fortran and C, respectively.
#
# If these macros look like an overkill to you, consider using
# ACX_COMPILER_*_SIMPLE macros (see below).

# ACX_COMPILER_FC_VENDOR()
# -----------------------------------------------------------------------------
# Detects the vendor of the Fortran compiler. The result is "unknown" or one of
# the vendor IDs (see _ACX_COMPILER_KNOWN_VENDORS(Fortran)).
#
# The result is cached in the acx_cv_fc_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_FC_VENDOR],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_PROVIDE_IFELSE([AC_FC_PP_SRCEXT], [],
     [m4_warn([syntax],
        [ACX_COMPILER_FC_VENDOR requires calling the Fortran compiler with ]dnl
[a preprocessor but no call to AC_FC_PP_SRCEXT is detected])])dnl
   AC_CACHE_CHECK([for Fortran compiler vendor], [acx_cv_fc_compiler_vendor],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^AMD clang version' >/dev/null 2>&1],
        [acx_cv_fc_compiler_vendor=amd],
        [_ACX_COMPILER_VENDOR])])])

# ACX_COMPILER_FC_VENDOR_SIMPLE()
# -----------------------------------------------------------------------------
# Detects the vendor of the Fortran compiler. The result is "intel", "nag",
# "portland", "cray", "nec", "gnu", "amd", "flang" or "unknown".
#
# This is a simplified version ACX_COMPILER_FC_VENDOR, which tries to detect
# the vendor based on the version output of the compiler, instead of checking
# whether vendor-specific macros are defined.
#
# The result is cached in the acx_cv_fc_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_FC_VENDOR_SIMPLE],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_provide([ACX_COMPILER_FC_VENDOR])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler vendor], [acx_cache_var],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^\(ifort\|ifx\) (IFORT)' >/dev/null 2>&1],
        [acx_cache_var=intel],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
grep '^NAG Fortran Compiler Release' >/dev/null 2>&1],
        [acx_cache_var=nag],
        [AS_VAR_GET([_AC_CC]) -V 2>/dev/null| dnl
grep '^Copyright.*\(The Portland Group\|NVIDIA CORPORATION\)' >/dev/null 2>&1],
        [acx_cache_var=portland],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | grep '^Cray Fortran' >/dev/null 2>&1],
        [acx_cache_var=cray],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^Copyright.*NEC Corporation' >/dev/null 2>&1],
        [acx_cache_var=nec],
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^GNU Fortran' >/dev/null 2>&1],
        [acx_cache_var=gnu],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^AMD clang version' >/dev/null 2>&1],
        [acx_cache_var=amd],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | grep '^f18 compiler' >/dev/null 2>&1],
        [acx_cache_var=flang],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^clang version' >/dev/null 2>&1],
        [acx_cache_var=flang],
        [acx_cache_var=unknown])
      rm -f a.out a.out.dSYM a.exe b.out])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_FC_VERSION()
# -----------------------------------------------------------------------------
# Detects the version of the C compiler. The result is either "unknown"
# or a string in the form "[epoch:]major[.minor[.patchversion]]", where
# "epoch:" is an optional prefix used in order to have an increasing version
# number in case of marketing change.
#
# The result is cached in the acx_cv_fc_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_FC_VERSION],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_REQUIRE([ACX_COMPILER_FC_VENDOR])_ACX_COMPILER_VERSION])

# ACX_COMPILER_FC_VERSION_SIMPLE()
# -----------------------------------------------------------------------------
# Detects the version of the C compiler. The result is either "unknown"
# or a string in the form "[epoch:]major[.minor[.patchversion]]", where
# "epoch:" is an optional prefix used in order to have an increasing version
# number in case of marketing change.
#
# This is a simplified version ACX_COMPILER_FC_VERSION, which tries to detect
# the version based on the version output of the compiler, instead of checking
# for the values of vendor-specific macros.
#
# The result is cached in the acx_cv_fc_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_FC_VERSION_SIMPLE],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_REQUIRE([ACX_COMPILER_FC_VENDOR])dnl
   m4_provide([ACX_COMPILER_FC_VERSION])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_version])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler version], [acx_cache_var],
     [AS_CASE([AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])],
        [intel],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^ifort (IFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
         AS_IF([test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^ifx (IFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
            AS_IF([test -n "$acx_cache_var"],
              [acx_cache_var="oneapi:${acx_cache_var}"])],
           [acx_cache_var="classic:${acx_cache_var}"])],
        [nag],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/^NAG Fortran Compiler Release \([0-9][0-9]*\.[0-9][0-9]*\).*]dnl
[Build \([0-9][0-9]*\)/\1.\2/p']`],
        [portland],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>/dev/null | dnl
[sed -n 's/\(pgfortran\|pgf90\) \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\2.\3/p']`
         AS_IF([test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>/dev/null | dnl
[sed -n 's/nvfortran \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\1.\2/p']`
            AS_IF([test -n "$acx_cache_var"],
              [acx_cache_var="nv:${acx_cache_var}"])],
           [acx_cache_var="pg:${acx_cache_var}"])],
        [cray],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*ersion \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`],
        [nec],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
[sed -n 's/^nfort (NFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`],
        [gnu],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpfullversion 2>/dev/null`
         AS_IF([test $? -ne 0 || test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpversion 2>/dev/null`])],
        [amd],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/.*AOCC_\([0-9][0-9]*\)[._]\([0-9][0-9]*\)[._]\([0-9][0-9]*\).*/\1.\2.\3/p']`],
        [flang],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*version \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\)/\1/p']`
         AS_IF([test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/.*version \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
            AS_IF([test -n "$acx_cache_var"],
              [acx_cache_var="classic:${acx_cache_var}"])],
           [acx_cache_var="f18:${acx_cache_var}"])],
        [acx_cache_var=unknown])
      rm -f a.out a.out.dSYM a.exe b.out
      AS_IF([test -z "$acx_cache_var"], [acx_cache_var=unknown])])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_CC_VENDOR()
# -----------------------------------------------------------------------------
# Detects the vendor of the C compiler. The result is "unknown" or one of
# the vendor IDs (see _ACX_COMPILER_KNOWN_VENDORS(C)).
#
# The result is cached in the acx_cv_c_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_CC_VENDOR],
  [AC_LANG_ASSERT([C])dnl
   AC_CACHE_CHECK([for C compiler vendor], [acx_cv_c_compiler_vendor],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
grep '^NAG Fortran Compiler Release' >/dev/null 2>&1],
        [acx_cv_c_compiler_vendor=nag],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^AMD clang version' >/dev/null 2>&1],
        [acx_cv_c_compiler_vendor=amd],
        [_ACX_COMPILER_VENDOR])])])

# ACX_COMPILER_CC_VENDOR_SIMPLE()
# -----------------------------------------------------------------------------
# Detects the vendor of the C compiler. The result is  "intel", "nag",
# "portland", "cray", "nec", "gnu", "amd", "clang" or "unknown".
#
# This is a simplified version ACX_COMPILER_CC_VENDOR, which tries to detect
# the vendor based on the version output of the compiler, instead of checking
# whether vendor-specific macros are defined.
#
# The result is cached in the acx_cv_c_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_CC_VENDOR_SIMPLE],
  [AC_LANG_ASSERT([C])dnl
   m4_provide([ACX_COMPILER_CC_VENDOR])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler vendor], [acx_cache_var],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^icc (ICC)' >/dev/null 2>&1],
        [acx_cache_var=intel],
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^Intel.*oneAPI.*Compiler' >/dev/null 2>&1],
        [acx_cache_var=intel],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
grep '^NAG Fortran Compiler Release' >/dev/null 2>&1],
        [acx_cache_var=nag],
        [AS_VAR_GET([_AC_CC]) -V 2>/dev/null| dnl
grep '^Copyright.*\(The Portland Group\|NVIDIA CORPORATION\)' >/dev/null 2>&1],
        [acx_cache_var=portland],
        [AS_VAR_GET([_AC_CC]) -V 2>&1 | grep '^Cray C' >/dev/null 2>&1],
        [acx_cache_var=cray],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^Cray clang' >/dev/null 2>&1],
        [acx_cache_var=cray],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^Copyright.*NEC Corporation' >/dev/null 2>&1],
        [acx_cache_var=nec],
        [AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
grep '^gcc' >/dev/null 2>&1],
        [acx_cache_var=gnu],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^AMD clang version' >/dev/null 2>&1],
        [acx_cache_var=amd],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^Apple \(LLVM\|clang\) version' >/dev/null 2>&1],
        [acx_cache_var=apple],
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^clang version' >/dev/null 2>&1],
        [acx_cache_var=clang],
        [acx_cache_var=unknown])
      rm -f a.out a.out.dSYM a.exe b.out])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_CC_VERSION()
# -----------------------------------------------------------------------------
# Detects the version of the C compiler. The result is either "unknown"
# or a string in the form "[epoch:]major[.minor[.patchversion]]", where
# "epoch:" is an optional prefix used in order to have an increasing version
# number in case of marketing change.
#
# The result is cached in the acx_cv_c_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_CC_VERSION],
  [AC_LANG_ASSERT([C])dnl
   AC_REQUIRE([ACX_COMPILER_CC_VENDOR])_ACX_COMPILER_VERSION])

# ACX_COMPILER_CC_VERSION_SIMPLE()
# -----------------------------------------------------------------------------
# Detects the version of the C compiler. The result is either "unknown"
# or a string in the form "[epoch:]major[.minor[.patchversion]]", where
# "epoch:" is an optional prefix used in order to have an increasing version
# number in case of marketing change.
#
# This is a simplified version ACX_COMPILER_CC_VERSION, which tries to detect
# the version based on the version output of the compiler, instead of checking
# for the values of vendor-specific macros.
#
# The result is cached in the acx_cv_c_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_CC_VERSION_SIMPLE],
  [AC_LANG_ASSERT([C])dnl
   AC_REQUIRE([ACX_COMPILER_CC_VENDOR])dnl
   m4_provide([ACX_COMPILER_CC_VERSION])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_version])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler version], [acx_cache_var],
     [AS_CASE([AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])],
        [intel],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^icc (ICC) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
         AS_IF([test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^Intel.*oneAPI.*Compiler \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
            AS_IF([test -n "$acx_cache_var"],
              [acx_cache_var="oneapi:${acx_cache_var}"])],
           [acx_cache_var="classic:${acx_cache_var}"])],
        [nag],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/^NAG Fortran Compiler Release \([0-9][0-9]*\.[0-9][0-9]*\).*]dnl
[Build \([0-9][0-9]*\)/\1.\2/p']`],
        [portland],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>/dev/null | dnl
[sed -n 's/pgcc \((.*) \)\?\([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\2.\3/p']`
         AS_IF([test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>/dev/null | dnl
[sed -n 's/nvc \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\1.\2/p']`
            AS_IF([test -n "$acx_cache_var"],
              [acx_cache_var="nv:${acx_cache_var}"])],
           [acx_cache_var="pg:${acx_cache_var}"])],
        [cray],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*ersion \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
         AS_IF([test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version | dnl
[sed -n 's/.*version \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
            AS_IF([test -n "$acx_cache_var"],
              [acx_cache_var="clang:${acx_cache_var}"])],
           [acx_cache_var="classic:${acx_cache_var}"])],
        [nec],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
[sed -n 's/^ncc (NCC) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`],
        [gnu],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpfullversion 2>/dev/null`
         AS_IF([test $? -ne 0 || test -z "$acx_cache_var"],
           [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpversion 2>/dev/null`])],
        [amd],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/.*AOCC_\([0-9][0-9]*\)[._]\([0-9][0-9]*\)[._]\([0-9][0-9]*\).*/\1.\2.\3/p']`],
        [apple],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -E -n 's/^Apple (LLVM|clang) version ([0-9]+\.[0-9]+\.[0-9]+).*/\2/p']`],
        [clang],
        [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/.*version \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`],
        [acx_cache_var=unknown])
      rm -f a.out a.out.dSYM a.exe b.out
      AS_IF([test -z "$acx_cache_var"], [acx_cache_var=unknown])])
   m4_popdef([acx_cache_var])])

# ACX_COMPILER_CXX_VENDOR()
# -----------------------------------------------------------------------------
# Detects the vendor of the C++ compiler. The result is "unknown" or one of
# the vendor IDs (see _ACX_COMPILER_KNOWN_VENDORS(C)).
#
# The result is cached in the acx_cv_cxx_compiler_vendor variable.
#
AC_DEFUN([ACX_COMPILER_CXX_VENDOR],
  [AC_LANG_ASSERT([C++])dnl
   AC_CACHE_CHECK([for C++ compiler vendor], [acx_cv_cxx_compiler_vendor],
     [AS_IF(
        [AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
grep '^AMD clang version' >/dev/null 2>&1],
        [acx_cv_cxx_compiler_vendor=amd],
        [_ACX_COMPILER_VENDOR])])])

# ACX_COMPILER_CXX_VERSION()
# -----------------------------------------------------------------------------
# Detects the version of the C compiler. The result is either "unknown"
# or a string in the form "[epoch:]major[.minor[.patchversion]]", where
# "epoch:" is an optional prefix used in order to have an increasing version
# number in case of marketing change.
#
# The result is cached in the acx_cv_cxx_compiler_version variable.
#
AC_DEFUN([ACX_COMPILER_CXX_VERSION],
  [AC_LANG_ASSERT([C++])dnl
   AC_REQUIRE([ACX_COMPILER_CXX_VENDOR])_ACX_COMPILER_VERSION])

# _ACX_COMPILER_KNOWN_VENDORS()
# -----------------------------------------------------------------------------
# Expands into a language-specific m4-quoted comma-separated list of pairs. The
# first value in each pair is the compiler vendor ID, the second value is a
# comma-separated list of the vendor-specific intrinsic preprocessor macros.
# By default, expands to m4_fatal with the message saying that _AC_LANG is not
# supported.
#
m4_define([_ACX_COMPILER_KNOWN_VENDORS],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_fatal([the list of ]_AC_LANG[ compiler vendors is undefined])])])

# _ACX_COMPILER_KNOWN_VENDORS(Fortran)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_COMPILER_KNOWN_VENDORS for Fortran language.
#
m4_define([_ACX_COMPILER_KNOWN_VENDORS(Fortran)],
[[[intel, [__INTEL_COMPILER,__INTEL_LLVM_COMPILER]],
  [cray, [_CRAYFTN]],
  [nec, [__NEC__]],
  [portland, [__PGI]],
  [nag, [NAGFOR]],
  [sun, [__SUNPRO_F95]],
  [gnu, [__GFORTRAN__]],
  [ibm, [__xlC__]],
  [amd],
  [flang, [__FLANG,__flang__]]]])

# _ACX_COMPILER_KNOWN_VENDORS(C)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_COMPILER_KNOWN_VENDORS for C language.
#
m4_define([_ACX_COMPILER_KNOWN_VENDORS(C)],
[[[intel, [__ICC,__ECC,__INTEL_COMPILER,__INTEL_LLVM_COMPILER]],
  [cray, [_CRAYC,__cray__]],
  [nec, [__NEC__]],
  [portland, [__PGI,__NVCOMPILER]],
  [amd],
  [ibm, [__xlc__,__xlC__,__IBMC__,__IBMCPP__]],
  [pathscale, [__PATHCC__,__PATHSCALE__]],
  [apple, [__apple_build_version__]],
  [clang, [__clang__]],
  [fujitsu, [__FUJITSU]],
  [sdcc, [SDCC,__SDCC]],
  [nag],
  [gnu, [__GNUC__]],
  [sun, [__SUNPRO_C,__SUNPRO_CC]],
  [hp, [__HP_cc,__HP_aCC]],
  [dec, [__DECC,__DECCXX,__DECC_VER,__DECCXX_VER]],
  [borland, [__BORLANDC__,__CODEGEARC__,__TURBOC__]],
  [comeau, [__COMO__]],
  [kai, [__KCC]],
  [lcc, [__LCC__]],
  [sgi, [__sgi,sgi]],
  [microsoft, [_MSC_VER]],
  [metrowerks, [__MWERKS__]],
  [watcom, [__WATCOMC__]],
  [tcc, [__TINYC__]]]])

# _ACX_COMPILER_KNOWN_VENDORS(C++)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_COMPILER_KNOWN_VENDORS for C++ language.
#
m4_copy([_ACX_COMPILER_KNOWN_VENDORS(C)], [_ACX_COMPILER_KNOWN_VENDORS(C++)])

# _ACX_COMPILER_VENDOR()
# -----------------------------------------------------------------------------
# Originally taken from Autoconf Archive where it is known as
# AX_COMPILER_VENDOR.
# -----------------------------------------------------------------------------
# Detects the vendor of the compiler based on its intrinsic preprocessor
# macro. The result is "unknown" or one of the vendor IDs
# (see _ACX_COMPILER_KNOWN_VENDORS).
#
# The result is stored in the acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor
# variable.
#
m4_define([_ACX_COMPILER_VENDOR],
  [acx_compiler_vendor_options=dnl
"m4_foreach([pair], _ACX_COMPILER_KNOWN_VENDORS,
    [m4_ifnblank(m4_quote(m4_shift(pair)), m4_n(m4_car(pair): m4_cdr(pair)))])dnl
unknown: CHOKEME"
      acx_success=no
      for acx_compiler_vendor_test in $acx_compiler_vendor_options; do
        AS_CASE([$acx_compiler_vendor_test],
          [*:], [acx_compiler_vendor_candidate=$acx_compiler_vendor_test
                 continue],
          [acx_compiler_vendor_macro_defs=dnl
"defined("`echo $acx_compiler_vendor_test | sed 's%,%) || defined(%g'`")"])
        AC_COMPILE_IFELSE(
          [AC_LANG_PROGRAM([],
[[#if !($acx_compiler_vendor_macro_defs)
      choke me
#endif]])],
          [acx_success=yes])
        test "x$acx_success" = xyes && break
      done
      AS_VAR_SET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor],
        [`echo $acx_compiler_vendor_candidate | cut -d: -f1`])])

# _ACX_COMPILER_VERSION()
# -----------------------------------------------------------------------------
# Originally taken from Autoconf Archive where it is known as
# AX_COMPILER_VERSION.
# -----------------------------------------------------------------------------
# Detects the version of the compiler. The result is either "unknown" or a
# string in the form "epoch:major.minor.patchversion", where "epoch:" is an
# optional prefix used in order to have an increasing version number in case of
# marketing change. The "patchversion" position might equal to "x" if the
# script fails to identify it correctly.
#
# The result is cached in the acx_cv_[]_AC_LANG_ABBREV[]_compiler_version
# variable.
#
m4_define([_ACX_COMPILER_VERSION],
  [m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_compiler_version])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler version], [acx_cache_var],
     [acx_cache_var=unknown
      case AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor]) in [#(]
        m4_foreach([pair], _ACX_COMPILER_KNOWN_VENDORS,
          [m4_ifdef([_ACX_COMPILER_VERSION_]m4_toupper(m4_car(pair))[(]_AC_LANG[)],
             [m4_car(pair)[)] m4_indir([_ACX_COMPILER_VERSION_]m4_toupper(m4_car(pair))[(]_AC_LANG[)]) ;; [#(]
        ],
             [m4_ifdef([_ACX_COMPILER_VERSION_]m4_toupper(m4_car(pair)),
                [m4_car(pair)[)] m4_indir([_ACX_COMPILER_VERSION_]m4_toupper(m4_car(pair))) ;; [#(]
        ])])])
        *[)] : acx_cache_var=unknown ;;
      esac])
   m4_popdef([acx_cache_var])])

# _ACX_COMPILER_VERSION_FROM_MACROS(MACRO_MAJOR,
#                                   MACRO_MINOR,
#                                   MACRO_PATCHLEVEL,
#                                   [PROLOGUE])
# -----------------------------------------------------------------------------
# Expands to a generic scripts that retrieves the C or C++ (or another
# language that supports AC_COMPUTE_INT) compiler version.
#
m4_define([_ACX_COMPILER_VERSION_FROM_MACROS],
  [acx_cache_var=unknown
   AC_COMPUTE_INT([acx_compiler_version_value],
     [$1], [$4], [acx_compiler_version_value=])
   AS_IF([test -n "$acx_compiler_version_value"],
     [acx_cache_var=$acx_compiler_version_value
      AC_COMPUTE_INT([acx_compiler_version_value],
        [$2], [$4], [acx_compiler_version_value=])
      AS_IF([test -n "$acx_compiler_version_value"],
        [AS_VAR_APPEND([acx_cache_var], [".$acx_compiler_version_value"])
         AC_COMPUTE_INT([acx_compiler_version_value],
           [$3], [$4], [acx_compiler_version_value=])
         AS_IF([test -n "$acx_compiler_version_value"],
           [AS_VAR_APPEND([acx_cache_var],
              [".$acx_compiler_version_value"])])])])])

# for GNU
m4_define([_ACX_COMPILER_VERSION_GNU(C)],
  [_ACX_COMPILER_VERSION_FROM_MACROS(
     [__GNUC__], [__GNUC_MINOR__], [__GNUC_PATCHLEVEL__])])
m4_copy([_ACX_COMPILER_VERSION_GNU(C)], [_ACX_COMPILER_VERSION_GNU(C++)])
m4_define([_ACX_COMPILER_VERSION_GNU(Fortran)],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpfullversion 2>/dev/null`
   AS_IF([test $? -ne 0 || test -z "$acx_cache_var"],
     [acx_cache_var=`AS_VAR_GET([_AC_CC]) -dumpversion 2>/dev/null`])
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for Intel
m4_define([_ACX_COMPILER_VERSION_INTEL(C)],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifdef __INTEL_LLVM_COMPILER
#else
      choke me
#endif]])],
     [acx_compiler_version_epoch='oneapi'
      _ACX_COMPILER_VERSION_FROM_MACROS(
        [__INTEL_LLVM_COMPILER/10000],
        [(__INTEL_LLVM_COMPILER%10000)/100],
        [__INTEL_LLVM_COMPILER%100])],
     [acx_compiler_version_epoch='classic'
      _ACX_COMPILER_VERSION_FROM_MACROS(
        [ACX_MACRO_MAJOR], [ACX_MACRO_MINOR], [ACX_MACRO_PATCHLEVEL], [[
#if __INTEL_COMPILER < 2021 || __INTEL_COMPILER == 202110 || __INTEL_COMPILER == 202111
# define ACX_MACRO_MAJOR __INTEL_COMPILER/100
# define ACX_MACRO_MINOR (__INTEL_COMPILER%100)/10
# ifdef __INTEL_COMPILER_UPDATE
#  define ACX_MACRO_PATCHLEVEL __INTEL_COMPILER_UPDATE
# else
#  define ACX_MACRO_PATCHLEVEL __INTEL_COMPILER%10
# endif
#else
# define ACX_MACRO_MAJOR __INTEL_COMPILER
# define ACX_MACRO_MINOR __INTEL_COMPILER_UPDATE
# define ACX_MACRO_PATCHLEVEL 0
#endif]])])
   AS_IF([test "x$acx_cache_var" != xunknown],
     [acx_cache_var="${acx_compiler_version_epoch}:${acx_cache_var}"])])
m4_copy([_ACX_COMPILER_VERSION_INTEL(C)], [_ACX_COMPILER_VERSION_INTEL(C++)])
m4_define([_ACX_COMPILER_VERSION_INTEL(Fortran)],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^ifort (IFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
   AS_IF([test -z "$acx_cache_var"],
     [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/^ifx (IFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
      AS_IF([test -n "$acx_cache_var"],
        [acx_cache_var="oneapi:${acx_cache_var}"])],
     [acx_cache_var="classic:${acx_cache_var}"])
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/^.*://' | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for NAG
m4_define([_ACX_COMPILER_VERSION_NAG],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/^NAG Fortran Compiler Release \([0-9][0-9]*\.[0-9][0-9]*\).*]dnl
[Build \([0-9][0-9]*\)/\1.\2/p']`
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for Cray
m4_define([_ACX_COMPILER_VERSION_CRAY(C)],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifdef __cray__
#else
      choke me
#endif]])],
     [acx_compiler_version_epoch='clang'
      _ACX_COMPILER_VERSION_FROM_MACROS(
        [__cray_major__], [__cray_minor__], [__cray_patchlevel__])],
     [acx_compiler_version_epoch='classic'
      _ACX_COMPILER_VERSION_FROM_MACROS(
        [_RELEASE_MAJOR], [_RELEASE_MINOR], [_RELEASE_PATCHLEVEL])])
   AS_IF([test "x$acx_cache_var" != xunknown],
     [acx_cache_var="${acx_compiler_version_epoch}:${acx_cache_var}"])])
m4_copy([_ACX_COMPILER_VERSION_CRAY(C)], [_ACX_COMPILER_VERSION_CRAY(C++)])
m4_define([_ACX_COMPILER_VERSION_CRAY(Fortran)],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*ersion \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for PGI
m4_define([_ACX_COMPILER_VERSION_PORTLAND(C)],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifdef __NVCOMPILER
#else
      choke me
#endif]])],
     [acx_compiler_version_epoch='nv'
      _ACX_COMPILER_VERSION_FROM_MACROS(
        [__NVCOMPILER_MAJOR__], [__NVCOMPILER_MINOR__],
        [__NVCOMPILER_PATCHLEVEL__])],
     [acx_compiler_version_epoch='pg'
      _ACX_COMPILER_VERSION_FROM_MACROS(
        [__PGIC__], [__PGIC_MINOR__], [__PGIC_PATCHLEVEL__])])
   AS_IF([test "x$acx_cache_var" != xunknown],
     [acx_cache_var="${acx_compiler_version_epoch}:${acx_cache_var}"])])
m4_copy([_ACX_COMPILER_VERSION_PORTLAND(C)],
  [_ACX_COMPILER_VERSION_PORTLAND(C++)])
m4_define([_ACX_COMPILER_VERSION_PORTLAND(Fortran)],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>/dev/null | dnl
[sed -n 's/\(pgfortran\|pgf90\) \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\2.\3/p']`
   AS_IF([test -z "$acx_cache_var"],
     [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>/dev/null | dnl
[sed -n 's/nvfortran \([0-9][0-9]*\.[0-9][0-9]*\)-\([0-9][0-9]*\).*/\1.\2/p']`
      AS_IF([test -n "$acx_cache_var"],
        [acx_cache_var="nv:${acx_cache_var}"])],
     [acx_cache_var="pg:${acx_cache_var}"])
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/^.*://' | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for LLVM
m4_define([_ACX_COMPILER_VERSION_CLANG(C)],
  [_ACX_COMPILER_VERSION_FROM_MACROS(
     [__clang_major__], [__clang_minor__], [__clang_patchlevel__])])
m4_copy([_ACX_COMPILER_VERSION_CLANG(C)], [_ACX_COMPILER_VERSION_CLANG(C++)])
m4_define([_ACX_COMPILER_VERSION_FLANG(Fortran)],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) -V 2>&1 | dnl
[sed -n 's/.*version \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\)/\1/p']`
   AS_IF([test $? -ne 0 || test -z "$acx_cache_var"],
     [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/.*version \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
      AS_IF([test $? -eq 0 && test -n "$acx_cache_var"],
        [acx_cache_var="classic:${acx_cache_var}"])],
     [acx_cache_var="f18:${acx_cache_var}"])
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/^.*://' | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for Apple Clang:
m4_copy([_ACX_COMPILER_VERSION_CLANG(C)], [_ACX_COMPILER_VERSION_APPLE(C)])
m4_copy([_ACX_COMPILER_VERSION_APPLE(C)], [_ACX_COMPILER_VERSION_APPLE(C++)])

# for Tiny CC
m4_define([_ACX_COMPILER_VERSION_TCC(C)],
  [_ACX_COMPILER_VERSION_FROM_MACROS(
     [__TINYC__/10000], [(__TINYC__%10000)/100], [__TINYC__%100])])

# for NEC
m4_define([_ACX_COMPILER_VERSION_NEC(C)],
  [_ACX_COMPILER_VERSION_FROM_MACROS(
     [__NEC_VERSION__/10000], [(__NEC_VERSION__%10000)/100],
     [__NEC_VERSION__%100])])
m4_copy([_ACX_COMPILER_VERSION_NEC(C)], [_ACX_COMPILER_VERSION_NEC(C++)])
m4_define([_ACX_COMPILER_VERSION_NEC(Fortran)],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>&1 | dnl
[sed -n 's/^nfort (NFORT) \([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p']`
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])

# for AMD (AOCC)
m4_define([_ACX_COMPILER_VERSION_AMD],
  [acx_cache_var=`AS_VAR_GET([_AC_CC]) --version 2>/dev/null | dnl
[sed -n 's/.*AOCC_\([0-9][0-9]*\)[._]\([0-9][0-9]*\)[._]\([0-9][0-9]*\).*/\1.\2.\3/p']`
   AS_IF([test dnl
"`echo $acx_cache_var | sed 's/@<:@0-9@:>@//g' 2>/dev/null`" != '..'],
     [acx_cache_var=unknown])])