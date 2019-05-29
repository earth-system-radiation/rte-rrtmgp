AC_PREREQ([2.69])
AC_INIT([RTE+RRTMGP], [1.0.0], [], [],
  [https://github.com/RobertPincus/rte-rrtmgp])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([rte/mo_rte_sw.F90])

AM_INIT_AUTOMAKE([1.16.1 foreign subdir-objects])
dnl Enable silent rules by default with a workaround for the case when the
dnl version constraint above is relaxed:
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE

AC_ARG_ENABLE([openacc],
  [AC_HELP_STRING([--enable-openacc],
     [enable OpenACC kernels @<:@default=no@:>@])], [],
  [enable_openacc=no])
AM_CONDITIONAL([ENABLE_OPENACC], [test x"$enable_openacc" = xyes])

AC_PROG_FC
AC_LANG([Fortran])
AC_FC_PP_SRCEXT([F90])
ACX_COMPILER_FC_VENDOR

AC_FC_LINE_LENGTH([unlimited])
ACX_FC_MODULE_NAMING([AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])])
AM_CONDITIONAL([FCMODUC], [test x"$acx_cv_fc_module_naming_upper" = xyes])
ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FCMODOUT], ["$acx_cv_fc_module_out_flag"])])
ACX_LANG_MACRO_FLAG([], [:])

AS_VAR_IF([enable_openacc], [yes],
  [ACX_LANG_OPENACC_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openacc_flag"])])
   ACX_LANG_MACRO_CHECK_VALUE([_OPENACC], [201711 201510 201308 201306 201111])
   AS_IF(
     [test x"$acx_macro_value" = xunknown],
     [AC_MSG_WARN([unable to detect supported OpenACC version (at least dnl
version 2.6 is required): the compilation of the OpenACC kernels might fail])],
     [test "$acx_macro_value" -lt 201711 >&AS_MESSAGE_LOG_FD 2>&1],
     [AC_MSG_ERROR([Fortran compiler does not support OpenACC 2.6 (201711): dnl
disable OpenACC support (--disable-openacc)])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenACC support by default: dnl
set FCFLAGS accordingly to disable it])])])

ACX_SHLIB_FC_RPATH_FLAG
ACX_LANG_PACKAGE_INIT([NetCDF], ["$FCMODINC"], [-L "$acx_cv_fc_rpath_flag"])
ACX_LANG_PACKAGE_INIT([GPTL], ["$FCMODINC"], [-L "$acx_cv_fc_rpath_flag"])

dnl Check for NetCDF
acx_have_netcdf=no
acx_save_FCFLAGS=$FCFLAGS
NETCDF_FCFLAGS=${NETCDF_FCFLAGS-$acx_fc_NetCDF_inc_search_args}
FCFLAGS="$NETCDF_FCFLAGS $acx_save_FCFLAGS"
ACX_FC_MODULE_CHECK([netcdf],
  [ACX_LANG_LIB_SEARCH([NETCDF_FCLIBS], [nf90_open],
     [ASX_PREPEND_LDFLAGS([$acx_fc_NetCDF_lib_search_args],
        [-lnetcdff], [-lnetcdff -lnetcdf])],
     [NETCDF_FCLIBS=$acx_cv_fc_lib_func_NF90_OPEN
      acx_have_netcdf=yes], [:],
     [AC_LANG_PROGRAM([],
[[      use netcdf, only: nf90_open, NF90_NOWRITE
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', NF90_NOWRITE, ncid)]])])], [:])
FCFLAGS=$acx_save_FCFLAGS
AS_VAR_IF([acx_have_netcdf], [no],
  [AC_MSG_WARN([NetCDF is not available: the test suite is disabled])
   NETCDF_FCFLAGS= ; NETCDF_FCLIBS=])
AM_CONDITIONAL([HAVE_NETCDF], [test x"$acx_have_netcdf" = xyes])

dnl Check for GPTL
acx_have_gptl=no
acx_save_FCFLAGS=$FCFLAGS
GPTL_FCFLAGS=${GPTL_FCFLAGS-$acx_fc_GPTL_inc_search_args}
FCFLAGS="$GPTL_FCFLAGS $acx_save_FCFLAGS"
ACX_FC_MODULE_CHECK([gptl],
  [ACX_LANG_LIB_SEARCH([GPTL_FCLIBS], [gptlstart],
     [ASX_PREPEND_LDFLAGS([$acx_fc_GPTL_lib_search_args],
        [-lgptlf], [-lgptlf -lgptl])],
     [GPTL_FCLIBS=$acx_cv_fc_lib_func_GPTLSTART
      acx_have_gptl=yes], [:],
     [AC_LANG_PROGRAM([],
[[      use gptl, only: gptlstart
      implicit none
      integer :: status
      status = gptlstart('conftest')]])])], [:])
FCFLAGS=$acx_save_FCFLAGS
AS_IF([test x"$acx_have_gptl" = xyes && dnl
test x"acx_cv_fc_macro_flag" != xunknown],
  [AS_VAR_APPEND([FCFLAGS], [" ${acx_cv_fc_macro_flag}USE_TIMING"])],
  [GPTL_FCFLAGS= ; GPTL_FCLIBS=])

AC_ARG_VAR([AR], [archiver command])
AC_CHECK_TOOL([AR], [ar], [:])
AS_VAR_IF([AR], [:], [AC_MSG_ERROR([the archiver command not found])])
AC_ARG_VAR([ARFLAGS], [archiver flags])
ARFLAGS=${ARFLAGS-cr}
AC_ARG_VAR([RANLIB], [archive indexer command])
AC_PROG_RANLIB

dnl Check for PYTHON
AC_ARG_VAR([PYTHON],
  [Python interpreter command (used for testing)])dnl
ACX_PROG_SEARCH([PYTHON], [python3 python],
  [# 'generate-output-file-templates.py' needs Python 3, once the script is
   # compatible with Python 2.7 or later, replace the following line with
   # acx_python_try="$acx_candidate -c 'import sys; exit(not(sys.version_info[0] == 3 or (sys.version_info[0] == 2 and sys.version_info[1] > 6)))' >&AS_MESSAGE_LOG_FD"
   acx_python_try="$acx_candidate -c 'import sys; exit(sys.version_info@<:@0@:>@ < 3)' >&AS_MESSAGE_LOG_FD"
   _AC_DO_VAR([acx_python_try])],
  [PYTHON=$acx_candidate],
  [PYTHON=false
   AC_MSG_WARN(
     [unable to find Python interpreter: the test suite is disabled])])

acx_have_python_mods=yes
AS_IF([test x"$PYTHON" != xfalse],
  [for acx_python_mod in 'netCDF4' 'xarray' 'dask.base'; do
     AC_MSG_CHECKING([whether Python module '$acx_python_mod' is available])
     acx_python_mod_try="$PYTHON -c 'import $acx_python_mod' >&AS_MESSAGE_LOG_FD"
     _AC_DO_VAR([acx_python_mod_try])
     AS_IF([test $? -eq 0],
       [AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no])
        acx_have_python_mods=no
        AC_MSG_WARN([Python module '$acx_python_mod' is not available: the dnl
test suite is disabled])
        break])
   done])
AM_CONDITIONAL([HAVE_PYTHON_MODS], [test x"$acx_have_python_mods" = xyes])

AC_CONFIG_FILES([
    examples/rfmip-clear-sky.test
], [chmod a+x "$ac_file"])

AC_CONFIG_FILES([
    Makefile
    examples/Makefile
    rrtmgp/Makefile
    rte/Makefile
])

AC_OUTPUT
