# ACX_FC_MODULE_IN_FLAG([ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is known as
# AC_FC_MODULE_FLAG.
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to specify module search paths.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_fc_module_in_flag variable, which may
# contain a significant trailing whitespace.
#
# Known flags:
# gfortran: -Idir, -I dir (-M dir, -Mdir (deprecated),
#                          -Jdir for writing)
# g95: -I dir (-fmod=dir for writing)
# SUN: -Mdir, -M dir (-moddir=dir for writing;
#                     -Idir for includes is also searched)
# HP: -Idir, -I dir (+moddir=dir for writing)
# IBM: -Idir (-qmoddir=dir for writing)
# Intel: -Idir -I dir (-mod dir for writing)
# Absoft: -pdir
# Lahey: -mod dir
# Cray: -module dir, -p dir (-J dir for writing)
#       -e m is needed to enable writing .mod files at all
# Compaq: -Idir
# NAGWare: -I dir
# PathScale: -I dir  (but -module dir is looked at first)
# Portland: -module dir (first -module also names dir for writing)
# Fujitsu: -Am -Idir (-Mdir for writing is searched first, then '.',
#                     then -I)
#                    (-Am indicates how module information is saved)
#
AC_DEFUN([ACX_FC_MODULE_IN_FLAG],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_CACHE_CHECK([for Fortran compiler flag needed to specify search paths dnl
for module files],
     [acx_cv_fc_module_in_flag],
     [acx_cv_fc_module_in_flag=unknown
      AS_MKDIR_P([conftest.dir])
      cd conftest.dir
      AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[[      module conftest_module
      implicit none
      public
      contains
      subroutine conftest_routine
      end subroutine
      end module]])],
        [cd ..
         acx_save_FCFLAGS=$FCFLAGS
         AC_LANG_CONFTEST([AC_LANG_PROGRAM([],
[[      use conftest_module, only : conftest_routine
      implicit none
      call conftest_routine()]])])
         for acx_flag in -M -I '-I ' '-M ' -p '-mod ' '-module ' '-Am -I'; do
           FCFLAGS="$acx_save_FCFLAGS ${acx_flag}conftest.dir dnl
dnl Add the flag twice to prevent matching an output flag.
${acx_flag}conftest.dir"
           AC_COMPILE_IFELSE([], [acx_cv_fc_module_in_flag=$acx_flag])
           test "x$acx_cv_fc_module_in_flag" != xunknown && break
         done
         rm -f conftest.$ac_ext
         FCFLAGS=$acx_save_FCFLAGS])
      rm -rf conftest.dir])
   AS_VAR_IF([acx_cv_fc_module_in_flag], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify search paths for module files])])], [$1])])

# ACX_FC_MODULE_OUT_FLAG([ACTION-IF-SUCCESS],
#                        [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is known as
# AC_FC_MODULE_OUTPUT_FLAG.
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to specify module output path.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_fc_module_out_flag variable, which may
# contain a significant trailing whitespace.
#
# See ACX_FC_MODULE_IN_FLAG for the known flags.
#
AC_DEFUN([ACX_FC_MODULE_OUT_FLAG],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_CACHE_CHECK([for Fortran compiler flag needed to specify output path dnl
for module files],
     [acx_cv_fc_module_out_flag],
     [acx_cv_fc_module_out_flag=unknown
      AS_MKDIR_P([conftest.dir/sub])
      cd conftest.dir
      acx_save_FCFLAGS=$FCFLAGS
      AC_LANG_CONFTEST([AC_LANG_PROGRAM([],
[[      use conftest_module, only : conftest_routine
      implicit none
      call conftest_routine()]])])
      mv conftest.$ac_ext sub/conftest.$ac_ext
      AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      module conftest_module
      implicit none
      public
      contains
      subroutine conftest_routine
      end subroutine
      end module]])])
      for acx_flag in -J '-J ' -fmod= -moddir= +moddir= -qmoddir= '-mdir ' dnl
'-mod ' '-module ' -M '-Am -M' '-e m -J '; do
        FCFLAGS="${acx_flag}sub $acx_save_FCFLAGS"
        AC_COMPILE_IFELSE([],
          [cd sub
           AC_COMPILE_IFELSE([], [acx_cv_fc_module_out_flag=$acx_flag])
           cd ..])
        test "x$acx_cv_fc_module_out_flag" != xunknown && break
      done
      FCFLAGS=$acx_save_FCFLAGS
      cd ..
      rm -rf conftest.dir])
   AS_VAR_IF([acx_cv_fc_module_out_flag], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify output path for module files])])], [$1])])

# ACX_FC_MODULE_NAMING([ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is known as
# AC_FC_MODULE_EXTENSION.
# -----------------------------------------------------------------------------
# Finds the Fortran compiler module file naming template.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_module_naming_upper and
# acx_cv_fc_module_naming_ext variables. If output module files have uppercase
# names, acx_cv_fc_module_naming_upper is "yes", and "no" otherwise. The
# acx_cv_fc_module_naming_ext variable stores the file extension without the
# leading dot. Either of the variables can have value "unknown". The result is
# successful only if both variables are detected.
#
AC_DEFUN([ACX_FC_MODULE_NAMING],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_MSG_CHECKING([for Fortran compiler module file naming template])
   AS_IF([AS_VAR_TEST_SET([acx_cv_fc_module_naming_upper]) && dnl
AS_VAR_TEST_SET([acx_cv_fc_module_naming_ext])],
     [AS_ECHO_N(["(cached) "]) >&AS_MESSAGE_FD],
     [AS_MKDIR_P([conftest.dir])
      cd conftest.dir
      AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[[      module conftest_module
      implicit none
      public
      contains
      subroutine conftest_routine
      end subroutine
      end module]])],
        [AS_CASE(["$acx_cv_fc_module_naming_upper"],
           [yes],
           [acx_tmp='CONFTEST_MODULE.*'
            acx_cv_fc_module_naming_ext=unknown],
           [no],
           [acx_tmp='conftest_module.*'
            acx_cv_fc_module_naming_ext=unknown],
           [AS_VAR_SET_IF([acx_cv_fc_module_naming_ext],
              [acx_tmp="CONFTEST_MODULE.$acx_cv_fc_module_naming_ext dnl
conftest_module.$acx_cv_fc_module_naming_ext"
               acx_cv_fc_module_naming_upper=unknown],
              [acx_tmp='CONFTEST_MODULE.* conftest_module.*'
               acx_cv_fc_module_naming_upper=unknown
               acx_cv_fc_module_naming_ext=unknown])])
         acx_tmp=`ls $acx_tmp 2>/dev/null`
         AS_IF([test 1 -eq `AS_ECHO(["$acx_tmp"]) | wc -l` 2>/dev/null],
           [AS_CASE(["$acx_tmp"],
              [CONFTEST_MODULE.*],
              [acx_cv_fc_module_naming_upper=yes
               acx_cv_fc_module_naming_ext=`echo $acx_tmp | dnl
sed -n 's,CONFTEST_MODULE\.,,p'`],
              [conftest_module.*],
              [acx_cv_fc_module_naming_upper=no
               acx_cv_fc_module_naming_ext=`echo $acx_tmp | dnl
sed -n 's,conftest_module\.,,p'`])])])
      cd ..
      rm -rf conftest.dir])
   AS_IF([test x"$acx_cv_fc_module_naming_upper" = xunknown || dnl
test x"$acx_cv_fc_module_naming_ext" = xunknown],
     [AC_MSG_RESULT([unknown])
      m4_default([$2], [AC_MSG_FAILURE([unable to detect Fortran compiler dnl
module file naming template])])],
     [AS_VAR_IF([acx_cv_fc_module_naming_upper], [yes],
        [AC_MSG_RESULT([NAME.$acx_cv_fc_module_naming_ext])],
        [AC_MSG_RESULT([name.$acx_cv_fc_module_naming_ext])])
      $1])])

# ACX_FC_MODULE_CHECK(MODULE-NAME,
#                     [ACTION-IF-SUCCESS],
#                     [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the Fortran module MODULE-NAME is available. The result is
# either "yes" or "no".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_module_[]AS_TR_CPP(MODULE-NAME)
# variable.
#
AC_DEFUN([ACX_FC_MODULE_CHECK],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_module_[]AS_TR_CPP([$1])])dnl
   AC_CACHE_CHECK([for Fortran module AS_TR_CPP([$1])], [acx_cache_var],
     [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[      use $1]])],
        [AS_VAR_SET([acx_cache_var], [yes])],
        [AS_VAR_SET([acx_cache_var], [no])])])
   AS_VAR_IF([acx_cache_var], [yes], [$2], [m4_default([$3],
     [AC_MSG_FAILURE([Fortran module 'AS_TR_CPP([$1])' is not available])])])
   m4_popdef([acx_cache_var])])

# ACX_FC_MODULE_PROC_CHECK(MODULE-NAME,
#                          PROCEDURE-NAME,
#                          [CALL-CODE = "      CALL PROCEDURE-NAME()"],
#                          [ACTION-IF-SUCCESS],
#                          [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the Fortran module procedure PROCEDURE-NAME from the module
# MODULE-NAME is available. The check is performed by linking a program that
# uses the module MODULE-NAME as "USE MODULE-NAME, ONLY : PROCEDURE-NAME"
# followed by the "IMPLICIT NONE" statement and the CALL-CODE (defaults to
# calling the PROCEDURE-NAME without parameters, which means that if
# PROCEDURE-NAME is a function or a subroutine with parameters, the CALL-CODE
# must be provided). The result is either "yes" or "no".
#
# If successful, runs ACTION-IF-SUCCESS (defaults to nothing), otherwise runs
# ACTION-IF-FAILURE (defaults to failing with an error message).
#
# The result is cached in the
# acx_cv_fc_module_proc_[]AS_TR_CPP(MODULE-NAME)_[]AS_TR_CPP(PROCEDURE-NAME)
# variable.
#
AC_DEFUN([ACX_FC_MODULE_PROC_CHECK],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var],
     [acx_cv_fc_module_proc_[]AS_TR_CPP([$1])_[]AS_TR_CPP([$2])])dnl
   AC_CACHE_CHECK([for Fortran procedure AS_TR_CPP([$2]) from module dnl
AS_TR_CPP([$1])],
     [acx_cache_var],
     [AC_LINK_IFELSE([AC_LANG_PROGRAM([],[[      use $1, only : $2
      implicit none]
m4_default([$3], [[      call $2 ()]])])],
        [AS_VAR_SET([acx_cache_var], [yes])],
        [AS_VAR_SET([acx_cache_var], [no])])])
   AS_VAR_IF([acx_cache_var], [yes], [$4], [m4_default([$5],
     [AC_MSG_FAILURE([Fortran module procedure 'AS_TR_CPP([$2])' from dnl
module 'AS_TR_CPP([$1])' is not available])])])
   m4_popdef([acx_cache_var])])
