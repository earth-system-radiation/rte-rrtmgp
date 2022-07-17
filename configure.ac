AC_PREREQ([2.69])
AC_INIT([rte+rrtmgp], [1.5], [], [],
  [https://github.com/earth-system-radiation/rte-rrtmgp])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([rte/mo_rte_sw.F90])

dnl Check Fortran compiler:
AC_PROG_FC
AC_LANG([Fortran])

ACX_FC_PP_SRCEXT([F90])

ACX_COMPILER_FC_VENDOR_SIMPLE
AC_SUBST([FC_VENDOR], [$acx_cv_fc_compiler_vendor])
ACX_COMPILER_FC_VERSION_SIMPLE
AC_SUBST([FC_VERSION], [$acx_cv_fc_compiler_version])

ACX_FC_LINE_LENGTH([unlimited])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])
   AC_SUBST([FCMODUC], ["$acx_cv_fc_module_naming_upper"])])

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FCMODOUT], ["$acx_cv_fc_module_out_flag"])])

ACX_FC_INCLUDE_FLAG([AC_SUBST([FCINCFLAG], ["$acx_cv_fc_ftn_include_flag"])])
ACX_FC_INCLUDE_ORDER(
  [AC_SUBST([FCINCORDER], ["$acx_cv_fc_ftn_include_order"])],
  [AC_SUBST([FCINCORDER], ['src,flg'])])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FCINCFLAG_PP], ["$acx_cv_fc_pp_include_flag"])])
ACX_FC_INCLUDE_ORDER_PP(
  [AC_SUBST([FCINCORDER_PP], ["$acx_cv_fc_pp_include_order"])],
  [AC_SUBST([FCINCORDER_PP], ['inc,flg'])])

ACX_LANG_MACRO_FLAG([AC_SUBST([FCDEF_PP], ["$acx_cv_fc_macro_flag"])])

ACX_FC_MODULE_CHECK([ISO_C_BINDING])
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
      AS_VAR_APPEND([FCFLAGS], [" ${FCDEF_PP}RTE_USE_CBOOL"])
      AC_COMPILE_IFELSE([], [acx_cv_fc_use_cbool=yes])
      FCFLAGS=$acx_save_FCFLAGS])
   rm -f conftest.$ac_ext])
AS_VAR_IF([acx_cv_fc_use_cbool], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FCDEF_PP}RTE_USE_CBOOL"])])

dnl Additional macro definitions to aid the dependency generator,
dnl e.g. compiler intrinsic macros:
AC_SUBST([DEPGEN_FCFLAGS], [''])

dnl Enable/disable silent rules:
AC_ARG_ENABLE([silent-rules],
  [AC_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_silent_rules=yes])],
  [enable_silent_rules=yes])

AC_SUBST([DEFAULT_VERBOSITY], [0])
AS_VAR_IF([enable_silent_rules], [no], [DEFAULT_VERBOSITY=1])

dnl Generate RPATH flags if needed:
AC_ARG_ENABLE([rpaths],
  [AC_HELP_STRING([--enable-rpaths],
     [add directories specified with -L flags in LDFLAGS and LIBS to the
runtime library search paths (RPATH) @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rpaths=yes])],
  [enable_rpaths=yes])

AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([acx_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   for acx_L_arg in $acx_L_args; do
     ASX_VAR_APPEND_UNIQ([LDFLAGS],
       ["$acx_cv_fc_rpath_flag$acx_L_arg"], [' '])
   done
   ACX_SHLIB_RPATH_FLAGS_CHECK])

AC_ARG_ENABLE([gpu],
[  --enable-gpu            enable GPU kernels. The value of the argument must
                          be one of the following @<:@default=no@:>@:
                            acc|yes  enable OpenACC support
                            mp       enable OpenMP support
                            no       disable GPU offload],
[AS_CASE(["$enableval"],
   [yes], [enable_gpu=acc],
   [acc|mp|no], [],
   [AC_MSG_ERROR([unexpected value for the argument
--enable-gpu='$enableval'; valid values are 'acc', 'mp', dnl
'yes' (same as 'acc'), 'no'])])],
[enable_gpu=no])

AS_CASE(["$enable_gpu"],
  [acc],
  [ACX_LANG_OPENACC_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openacc_flag"])])
   acx_cv_fc_macro__OPENACC_defined=${acx_cv_fc_macro__OPENACC_defined-yes}
   ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenMP support by default: dnl
either enable GPU kernels with OpenMP support (--enable-gpu=mp) or set dnl
FCFLAGS accordingly to disable OpenMP])])],
dnl
  [mp],
  [ACX_LANG_OPENMP_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openmp_flag"])])
   acx_cv_fc_macro__OPENMP_defined=${acx_cv_fc_macro__OPENMP_defined-yes}
   ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenACC support by default: dnl
either enable GPU kernels with OpenACC support (--enable-gpu=acc) or set dnl
FCFLAGS accordingly to disable OpenACC])])],
dnl
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])])

AS_VAR_IF([acx_cv_fc_macro__OPENACC_defined], [yes],
  [ACX_LANG_MACRO_CHECK_VALUE([_OPENACC], [201711 201510 201308 201306 201111])
   AS_VAR_IF([acx_macro_value], [unknown],
     [AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FCDEF_PP}_OPENACC"])],
     [AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FCDEF_PP}_OPENACC=$acx_macro_value"])])
   AS_VAR_IF([enable_gpu], [no],
     [AC_MSG_WARN([GPU kernels are disabled but Fortran compiler defines dnl
_OPENACC macro: the compilation might be inconsistent])],
     [AS_IF([test x"$acx_cv_fc_macro__OPENACC_value" = xunknown],
        [AC_MSG_WARN([unable to detect supported OpenACC version (at least dnl
version 2.6 is required): the compilation of GPU kernels might fail])],
     [test "$acx_cv_fc_macro__OPENACC_value" -lt 201711 dnl
>&AS_MESSAGE_LOG_FD 2>&1],
     [AC_MSG_WARN([Fortran compiler does not support OpenACC 2.6 (201711): dnl
the compilation of GPU kernels might fail])])])])

AS_VAR_IF([acx_cv_fc_macro__OPENMP_defined], [yes],
  [ACX_LANG_MACRO_CHECK_VALUE([_OPENMP],
     [201511 201307 201107 200805 200505 200011 199911 199710])
   AS_VAR_IF([acx_macro_value], [unknown],
     [AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FCDEF_PP}_OPENMP"])],
     [AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FCDEF_PP}_OPENMP=$acx_macro_value"])])
   AS_VAR_IF([enable_gpu], [no],
     [AC_MSG_WARN([GPU kernels are disabled but Fortran compiler defines dnl
_OPENMP macro: the compilation might be inconsistent])],
     [AS_IF([test x"$acx_cv_fc_macro__OPENMP_value" = xunknown],
        [AC_MSG_WARN([unable to detect supported OpenMP version (at least dnl
version 4.0 is required): the compilation of GPU kernels might fail])],
     [test "$acx_cv_fc_macro__OPENMP_value" -lt 201307 dnl
>&AS_MESSAGE_LOG_FD 2>&1],
     [AC_MSG_WARN([Fortran compiler does not support OpenMP 4.0 (201307): dnl
the compilation of GPU kernels might fail])])])])

AS_VAR_IF([enable_gpu], [no],
  [AC_SUBST([GPU_ENABLED], ['#'])
   AC_SUBST([GPU_DISABLED], [''])],
  [AC_SUBST([GPU_ENABLED], [''])
   AC_SUBST([GPU_DISABLED], ['#'])])

AC_ARG_ENABLE([examples],
  [AC_HELP_STRING([--enable-examples],
     [enable examples @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_examples=yes])],
  [enable_examples=no])

AC_ARG_ENABLE([tests],
  [AC_HELP_STRING([--enable-tests],
     [enable tests @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_tests=yes])],
  [enable_tests=auto])

dnl Check NetCDF-Fortran:
netcdf_fortran_requested=no
AC_ARG_WITH([netcdf-fortran],
  [AC_HELP_STRING([--with-netcdf-fortran=<yes|no|directory>],
     [location of NetCDF Fortran library (lib and include subdirs)])],
     [AS_IF([test x"$withval" != xno], [netcdf_fortran_requested=yes])],
  [with_netcdf_fortran=yes])
AC_ARG_VAR([NETCDF_FCFLAGS],
  [exact Fortran compiler flags enabling NetCDF Fortran])
AC_ARG_VAR([NETCDF_FCLIBS],
  [exact linker flags enabling NetCDF Fortran when linking with Fortran ]dnl
[compiler])

AS_IF([test x"$with_netcdf_fortran" != xno],
  [acx_inc_search_args=
   acx_lib_search_args=
   AS_IF([test x"$with_netcdf_fortran" != xyes],
     [acx_inc_search_args="${FCMODINC}$with_netcdf_fortran/include"
      acx_lib_search_args="-L$with_netcdf_fortran/lib"
      AS_VAR_IF([enable_rpaths], [yes],
        [AS_VAR_APPEND([acx_lib_search_args],
           [" ${acx_cv_fc_rpath_flag}$with_netcdf_fortran/lib"])])])
   acx_save_FCFLAGS=$FCFLAGS
   NETCDF_FCFLAGS=${NETCDF_FCFLAGS-$acx_inc_search_args}
   FCFLAGS="$NETCDF_FCFLAGS $acx_save_FCFLAGS"
   ACX_FC_MODULE_CHECK([netcdf],
     [ACX_LANG_LIB_SEARCH([NETCDF_FCLIBS], [nf90_open],
        [ASX_PREPEND_LDFLAGS([$acx_lib_search_args],
           [-lnetcdff], [-lnetcdff -lnetcdf])],
        [NETCDF_FCLIBS=$acx_cv_fc_lib_func_NF90_OPEN],
        [with_netcdf_fortran=no],
        [AC_LANG_PROGRAM([],
[[      use netcdf, only: nf90_open, NF90_NOWRITE
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', NF90_NOWRITE, ncid)]])])],
     [with_netcdf_fortran=no])
   FCFLAGS=$acx_save_FCFLAGS])

AS_VAR_IF([with_netcdf_fortran], [no],
  [AS_VAR_IF([netcdf_fortran_requested], [yes],
     [AC_MSG_FAILURE([NetCDF Fortran is not available])])
   AS_VAR_IF([enable_examples], [yes],
     [AC_MSG_FAILURE([NetCDF Fortran is not available: unable to build the dnl
examples])])
   AS_IF(
     [test x"$enable_tests" = xyes],
     [AC_MSG_FAILURE([NetCDF Fortran is not available: unable to build the dnl
tests])],
     [test x"$enable_tests" = xauto],
     [AC_MSG_WARN([NetCDF Fortran is not available: unable to build the tests])
      enable_tests=no])
   NETCDF_FCFLAGS=;NETCDF_FCLIBS=])

dnl Check GPTL-Fortran:
gptl_fortran_requested=no
AC_ARG_WITH([gptl-fortran],
  [AC_HELP_STRING([--with-gptl-fortran=<yes|no|directory>],
     [location of GPTL Fortran library (lib and include subdirs)])],
  [AS_IF([test x"$withval" != xno], [gptl_fortran_requested=yes])],
  [with_gptl_fortran=no])
AC_ARG_VAR([GPTL_FCFLAGS],
  [exact Fortran compiler flags enabling GPTL Fortran])
AC_ARG_VAR([GPTL_FCLIBS],
  [exact linker flags enabling GPTL Fortran when linking with Fortran ]dnl
[compiler])

AS_IF([test x"$with_gptl_fortran" != xno],
  [acx_inc_search_args=
   acx_lib_search_args=
   AS_IF([test x"$with_gptl_fortran" != xyes],
     [acx_inc_search_args="${FCMODINC}$with_gptl_fortran/include"
      acx_lib_search_args="-L$with_gptl_fortran/lib"
      AS_VAR_IF([enable_rpaths], [yes],
        [AS_VAR_APPEND([acx_lib_search_args],
           [" ${acx_cv_fc_rpath_flag}$with_gptl_fortran/lib"])])])
   acx_save_FCFLAGS=$FCFLAGS
   GPTL_FCFLAGS=${GPTL_FCFLAGS-$acx_inc_search_args}
   FCFLAGS="$GPTL_FCFLAGS $acx_save_FCFLAGS"
   ACX_FC_MODULE_CHECK([gptl],
     [ACX_LANG_LIB_SEARCH([GPTL_FCLIBS], [gptlstart],
        [ASX_PREPEND_LDFLAGS([$acx_lib_search_args],
           [-lgptlf], [-lgptlf -lgptl])],
        [GPTL_FCLIBS=$acx_cv_fc_lib_func_GPTLSTART],
        [with_gptl_fortran=no],
        [AC_LANG_PROGRAM([],
[[      use gptl, only: gptlstart
      implicit none
      integer :: status
      status = gptlstart('conftest')]])])],
     [with_gptl_fortran=no])
   FCFLAGS=$acx_save_FCFLAGS])

AS_VAR_IF([with_gptl_fortran], [no],
  [AS_VAR_IF([gptl_fortran_requested], [yes],
     [AC_MSG_FAILURE([GPTL Fortran is not available])])
   GPTL_FCFLAGS=;GPTL_FCLIBS=],
  [AS_VAR_APPEND([FCFLAGS], [" ${FCDEF_PP}USE_TIMING"])])

dnl Check AR and ARFLAGS:
AC_ARG_VAR([AR], [archiver command])
AC_ARG_VAR([ARFLAGS], [archiver flags])
ARFLAGS=${ARFLAGS-cr}
ACX_PROG_SEARCH([AR], [ar],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
     [acx_ar_try="$acx_candidate ${ARFLAGS} libconftest.a dnl
conftest.$ac_objext >&AS_MESSAGE_LOG_FD"
      _AC_DO_VAR([acx_ar_try])
      rm -f libconftest.a])
   test "$ac_status" -eq 0],
  [AR=$acx_cv_prog_AR])

dnl Check PYTHON:
AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python3 python],
  [acx_python_try="$acx_candidate $srcdir/mkhelper/depgen.py -h dnl
>&AS_MESSAGE_LOG_FD"
   _AC_DO_VAR([acx_python_try])],
  [PYTHON=$acx_cv_prog_PYTHON])

dnl We need Python 3 for the tests:
AS_IF([test x"$enable_tests" != xno],
  [AC_CACHE_CHECK([whether $PYTHON is compatible with the tests],
     [acx_cv_python_test_compatible],
     [acx_cv_python_test_compatible=no
      acx_python_try="$PYTHON -c dnl
'import sys; import sys; exit(sys.version_info@<:@0:2@:>@ < (3, 8))' dnl
>&AS_MESSAGE_LOG_FD"
      _AC_DO_VAR([acx_python_try])
      AS_IF([test $? -eq 0], [acx_cv_python_test_compatible=yes])])
   AS_VAR_IF([acx_cv_python_test_compatible], [no],
     [AS_VAR_IF([enable_tests], [auto],
        [AC_MSG_WARN([$PYTHON is not compatible with the tests])
         enable_tests=no],
        [AC_MSG_FAILURE([$PYTHON is not compatible with the tests])])])])

dnl We need some additional Python modules for the tests:
AS_IF([test x"$enable_tests" != xno],
  [set 'netCDF4'           'import netCDF4' \
       'xarray (v0.12.2+)' 'import xarray; exit(tuple(map(int, xarray.__version__.split("."))) < (0, 12, 2))' \
       'dask.array'        'import dask.array'
   while test $[]@%:@ -gt 1; do
     m4_pushdef([acx_cache_var], [acx_cv_has_pymodule_[]AS_TR_SH([$][1])])dnl
     AC_CACHE_CHECK([whether Python module '$[]1' is available],
       [acx_cache_var],
       [AS_VAR_SET([acx_cache_var], [no])
        acx_pymodule_try="$PYTHON -c '$[]2' >&AS_MESSAGE_LOG_FD"
        _AC_DO_VAR([acx_pymodule_try])
        AS_IF([test $? -eq 0], [AS_VAR_SET([acx_cache_var], [yes])])])
     AS_VAR_IF([acx_cache_var], [no],
       [AS_VAR_IF([enable_tests], [auto],
          [AC_MSG_WARN([Python module '$[]1' is not available: dnl
unable to run the tests])
           enable_tests=no],
          [AC_MSG_FAILURE([Python module '$[]1' is not available: dnl
unable to run the tests])])
        break])
     shift; shift
     m4_popdef([acx_cache_var])dnl
   done
   set dummy; shift])

AS_VAR_IF([enable_examples], [no],
  [AC_SUBST([EXAMPLES_ENABLED], ['#'])
   AC_SUBST([EXAMPLES_DISABLED], [''])],
  [AC_SUBST([EXAMPLES_ENABLED], [''])
   AC_SUBST([EXAMPLES_DISABLED], ['#'])])

AS_VAR_IF([enable_tests], [no],
  [AC_SUBST([TESTS_ENABLED], ['#'])
   AC_SUBST([TESTS_DISABLED], [''])],
  [AC_SUBST([TESTS_ENABLED], [''])
   AC_SUBST([TESTS_DISABLED], ['#'])])
AC_SUBST([testcachedir], ['.testcache'])

AC_ARG_VAR([LAUNCH_PREFIX], [execution prefix for the tests ]dnl
[(e.g. interactive job submission command)])

AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_INSTALL

AC_CONFIG_FILES([
   examples/all-sky/all-sky.test
   examples/rfmip-clear-sky/rfmip-clear-sky.test
   tests/clear_sky_regression.test
   tests/test_zenith_angle_spherical_correction.test
], [chmod a+x "$ac_file"])

AC_CONFIG_FILES([Makefile])

dnl An additional step to generate dummy dependencies:
AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch Makefile])

AC_OUTPUT
