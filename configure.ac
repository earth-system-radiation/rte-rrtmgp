AC_PREREQ([2.69])
AC_INIT([rte+rrtmgp], [1.7], [], [],
  [https://github.com/earth-system-radiation/rte-rrtmgp])
AC_CONFIG_MACRO_DIR([m4])
m4_pattern_forbid([^A[SC]X_])dnl
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([rte-frontend/mo_rte_sw.F90])

AM_INIT_AUTOMAKE([1.16 foreign subdir-objects])
AM_MAINTAINER_MODE([disable])
AM_SILENT_RULES([yes])
AM_MAKE_INCLUDE

AC_PROG_FC
AC_LANG([Fortran])
ACX_COMPILER_FC_VENDOR_SIMPLE
ACX_COMPILER_FC_VERSION_SIMPLE
m4_pushdef([AC_FC_SRCEXT], [ACX_FC_PP_SRCEXT($@)])dnl
AC_FC_SRCEXT([F90])
m4_popdef([AC_FC_SRCEXT])dnl

ACX_FC_LINE_LENGTH([unlimited])

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FC_MOD_IN], [$acx_cv_fc_module_in_flag])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FC_MOD_OUT], [$acx_cv_fc_module_out_flag])])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FC_MOD_FILE_UPPER], [$acx_cv_fc_module_naming_upper])
   AC_SUBST([FC_MOD_FILE_EXT], [$acx_cv_fc_module_naming_ext])])

ACX_FC_INCLUDE_FLAG(
  [AC_SUBST([FC_INC_FLAG], [$acx_cv_fc_ftn_include_flag])])
    ACX_FC_INCLUDE_ORDER(
      [FC_INC_ORDER=$acx_cv_fc_ftn_include_order],
      [FC_INC_ORDER=src,flg])
AC_SUBST([FC_INC_ORDER])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FC_INC_FLAG_PP_f90], [$acx_cv_fc_pp_include_flag])])
ACX_FC_INCLUDE_ORDER_PP(
  [FC_INC_ORDER_PP_f90=$acx_cv_fc_pp_include_order],
  [FC_INC_ORDER_PP_f90=inc,flg])
AC_SUBST([FC_INC_ORDER_PP_f90])

ACX_LANG_MACRO_FLAG([AC_SUBST([FC_PP_DEF], [$acx_cv_fc_macro_flag])])

ACX_FC_MODULE_CHECK([ISO_C_BINDING])

AC_ARG_ENABLE([rpaths],
  [AS_HELP_STRING([--enable-rpaths],
     [add directories specified with -L flags in LDFLAGS and LIBS to the
runtime library search paths (RPATH) @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rpaths=yes])],
  [enable_rpaths=yes])

acx_lib_search_flags='-L'
AC_SUBST([RPATH_FCFLAGS], [''])
AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([acx_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   for acx_L_arg in $acx_L_args; do
     ASX_VAR_APPEND_UNIQ([RPATH_FCFLAGS],
       ["$acx_cv_fc_rpath_flag$acx_L_arg"], [' '])
   done
   ACX_SHLIB_RPATH_FLAGS_CHECK([$RPATH_FCFLAGS])
   AS_VAR_APPEND([acx_lib_search_flags], [" '$acx_cv_fc_rpath_flag'"])])

AC_CACHE_CHECK([whether switching to C_BOOL type is required],
  [acx_cv_fc_use_cbool],
  [acx_cv_fc_use_cbool=unknown
   AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      subroutine conftest_foo(a) bind(C)
      use iso_c_binding
      implicit none
#ifdef RTE_USE_CBOOL
      integer, parameter :: wl = c_bool
#else
      integer, parameter :: wl = kind(.true.)
#endif
      logical(wl) :: a
      end subroutine]])])
   AC_COMPILE_IFELSE([],
     [acx_cv_fc_use_cbool=no],
     [acx_save_FCFLAGS=$FCFLAGS
      AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}RTE_USE_CBOOL"])
      AC_COMPILE_IFELSE([], [acx_cv_fc_use_cbool=yes])
      FCFLAGS=$acx_save_FCFLAGS])
   rm -f conftest.$ac_ext])
AM_CONDITIONAL([USE_CBOOL], [test "x$acx_cv_fc_use_cbool" = xyes])

AC_ARG_ENABLE([single-precision],
  [AC_HELP_STRING([--enable-single-precision],
     [enable single-precision floating-point model @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_single_precision=yes])],
  [enable_single_precision=no])

AS_VAR_IF([enable_single_precision], [no],
  [ACX_LANG_MACRO_CHECK_DEFINED([RTE_USE_SP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Single-precision floating-point model is enabled via dnl
FCFLAGS: either enable it via the configure option dnl
(--enable-single-precision) or remove the definition of the RTE_USE_SP dnl
preprocessor macro from FCFLAGS])])])
AM_CONDITIONAL([SP_ENABLED], [test "x$enable_single_precision" = xyes])

AC_SUBST([DEPGEN_FCFLAGS], [''])

AC_ARG_ENABLE([gpu],
[  --enable-gpu            enable GPU kernels. The value of the argument must
                          be one of the following @<:@default=no@:>@:
                            openacc|yes  enable OpenACC support
                            openmp       enable OpenMP support
                            no           disable GPU offload],
[AS_CASE(["$enableval"],
   [yes], [enable_gpu=openacc],
   [openacc|openmp|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-gpu='$enableval'; valid values are 'openacc', 'openmp', 'yes' dnl
(same as 'openacc'), 'no'])])],
[enable_gpu=no])

AS_CASE(["$enable_gpu"],
  [openacc],
  [ACX_LANG_OPENACC_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openacc_flag"])])
   ACX_LANG_MACRO_CHECK_VALUE([_OPENACC],
     [202111 202011 201911 201811 201711 201510 201308 201306 201111])
   AS_VAR_IF([acx_macro_value], [unknown],
     [AC_MSG_WARN([unable to detect supported OpenACC version (at least dnl
version 2.6 is required): the compilation of GPU kernels might fail])
      AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENACC"])],
     [AS_IF([test "$acx_macro_value" -lt 201711 >&AS_MESSAGE_LOG_FD 2>&1],
        [AC_MSG_WARN([Fortran compiler does not support OpenACC 2.6 dnl
(201711): the compilation of GPU kernels might fail])])
      AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FC_PP_DEF}_OPENACC=$acx_macro_value"])])
   ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenMP support by default: dnl
either enable GPU kernels with OpenMP support (--enable-gpu=openmp) or set dnl
FCFLAGS accordingly to disable OpenMP])])],
dnl
  [openmp],
  [ACX_LANG_OPENMP_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openmp_flag"])])
   ACX_LANG_MACRO_CHECK_VALUE([_OPENMP],
     [202111 202011 201811 201511 201307 201107 200805 200505 200011 199911 dnl
199710])
   AS_VAR_IF([acx_macro_value], [unknown],
     [AC_MSG_WARN([unable to detect supported OpenMP version (at least dnl
version 4.0 is required): the compilation of GPU kernels might fail])
      AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENMP"])],
     [AS_IF([test "$acx_macro_value" -lt 201711 >&AS_MESSAGE_LOG_FD 2>&1],
        [AC_MSG_WARN([Fortran compiler does not support OpenMP 4.0 (201307): dnl
the compilation of GPU kernels might fail])])
      AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FC_PP_DEF}_OPENMP=$acx_macro_value"])])
   ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenACC support by default: dnl
either enable GPU kernels with OpenACC support (--enable-gpu=openacc) or set dnl
FCFLAGS accordingly to disable OpenACC])])],
dnl
  [for acx_macro in _OPENACC _OPENMP; do
     ACX_LANG_MACRO_CHECK_DEFINED([$acx_macro])
     AS_VAR_IF([acx_macro_defined], [yes],
       [AC_MSG_WARN([GPU kernels are disabled but Fortran compiler defines dnl
$acx_macro macro: the compilation might be inconsistent])])
   done])
AM_CONDITIONAL([GPU_ENABLED], [test "x$enable_gpu" != xno])

AC_ARG_ENABLE([tests],
  [AC_HELP_STRING([--enable-tests], [enable tests @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_tests=yes])],
  [enable_tests=auto])

eval "set dummy $acx_lib_search_flags; shift"
ACX_LANG_PACKAGE_INIT([NetCDF], ["${FC_MOD_IN}"], ["$@"])

AS_IF([test x"$enable_tests" != xno],
  [acx_save_FCFLAGS=$FCFLAGS
   NETCDF_FCFLAGS=${NETCDF_FCFLAGS-$acx_fc_NetCDF_inc_search_args}
   FCFLAGS="$NETCDF_FCFLAGS $acx_save_FCFLAGS"
   ACX_FC_MODULE_CHECK([netcdf],
     [ACX_LANG_LIB_SEARCH([NETCDF_FCLIBS], [nf90_open],
        [ASX_PREPEND_LDFLAGS([$acx_fc_NetCDF_lib_search_args],
           [-lnetcdff], [-lnetcdff -lnetcdf])],
        [NETCDF_FCLIBS=$acx_cv_fc_lib_func_NF90_OPEN],
        [AS_VAR_IF([enable_tests], [yes],
           [AC_MSG_FAILURE([Fortran function NF90_OPEN is not available])],
           [enable_tests=no])],
        [AC_LANG_PROGRAM([],
[[      use netcdf, only: nf90_open, NF90_NOWRITE
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', NF90_NOWRITE, ncid)]])])],
     [AS_VAR_IF([enable_tests], [yes],
        [AC_MSG_FAILURE([Fortran module NETCDF is not available])],
        [enable_tests=no])])
   FCFLAGS=$acx_save_FCFLAGS])

AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python3 python],
  [acx_python_try="$acx_candidate $srcdir/mkhelper/depgen.py -h >/dev/null"
   _AC_DO_VAR([acx_python_try])],
  [PYTHON=$acx_cv_prog_PYTHON])

AS_IF([test x"$enable_tests" != xno],
  [AC_CACHE_CHECK([whether $PYTHON is compatible with the tests],
     [acx_cv_python_test_compatible],
     [acx_cv_python_test_compatible=no
      acx_python_try="$PYTHON -c dnl
'import sys; exit(sys.version_info@<:@0:2@:>@ < (3, 8))' >/dev/null"
      _AC_DO_VAR([acx_python_try])
      AS_IF([test $? -eq 0], [acx_cv_python_test_compatible=yes])])
   AS_VAR_IF([acx_cv_python_test_compatible], [no],
     [AS_VAR_IF([enable_tests], [yes],
        [AC_MSG_FAILURE([$PYTHON is not compatible with the tests])],
        [enable_tests=no])])])

AS_IF([test x"$enable_tests" != xno],
  [set 'numpy' 'import numpy' \
       'netCDF4 or h5netcdf+scipy' 'try:
    import netCDF4
except:
    import h5netcdf
    import scipy' \
       'xarray>=0.12.2' 'import xarray; exit(tuple(map(int, xarray.__version__.split("."))) < (0, 12, 2))' \
       'dask.array' 'import dask.array'
   while test $[]@%:@ -gt 1; do
     m4_pushdef([acx_cache_var], [acx_cv_has_pymodule_[]AS_TR_SH([$][1])])dnl
     AC_CACHE_CHECK([whether Python module $[]1 is available],
       [acx_cache_var],
       [AS_VAR_SET([acx_cache_var], [no])
        acx_pymodule_try="$PYTHON -c '$[]2' >/dev/null"
        _AC_DO_VAR([acx_pymodule_try])
        AS_IF([test $? -eq 0], [AS_VAR_SET([acx_cache_var], [yes])])])
     AS_VAR_IF([acx_cache_var], [no],
       [AS_VAR_IF([enable_tests], [yes],
          [AC_MSG_FAILURE([Python module $[]1 is not available: unable to dnl
run the tests])],
          [enable_tests=no])
        break])
     shift; shift
     m4_popdef([acx_cache_var])dnl
   done
   set dummy; shift])

AM_CONDITIONAL([TESTS_ENABLED], [test "x$enable_tests" != xno])
AC_ARG_VAR([LAUNCH_PREFIX], [execution prefix for the tests (e.g. ]dnl
[interactive job submission command)])
AC_SUBST([testcachedir], ['.testcache'])
AM_MISSING_PROG([CURL], [curl])

AC_PROG_RANLIB
AC_PROG_LN_S

AC_CONFIG_FILES([
   examples/all-sky/all-sky.test
   examples/rfmip-clear-sky/rfmip-clear-sky.test
   tests/check_equivalence.test
   tests/check_variants.test
   tests/test_zenith_angle_spherical_correction.test
], [chmod a+x "$ac_file"])

AC_CONFIG_FILES([Makefile])

AC_CONFIG_COMMANDS([fortran_dependencies],
  [touch Makefile.exe.deps Makefile.lib.deps; sleep 1; touch Makefile])

AC_OUTPUT
