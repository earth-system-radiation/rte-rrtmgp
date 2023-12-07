# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# ACX_FC_MODULE_IN_FLAG([ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to specify module search paths.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_fc_module_in_flag variable, which may
# contain a significant trailing whitespace.
#
# The implementation patches the standard Autoconf macro AC_FC_MODULE_FLAG to
# reduce the number of LANG switches and to avoid false negative results with
# the GFortran '-fmodule-private' flag.
#
AC_DEFUN([ACX_FC_MODULE_IN_FLAG],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([ac_cv_fc_module_flag], [acx_cv_fc_module_in_flag])dnl
   m4_pushdef([AC_CACHE_CHECK],
     m4_bpatsubst(m4_dquote(m4_defn([AC_CACHE_CHECK])),
       [\$][1],
       [for Fortran compiler flag needed to specify search paths for module dnl
files]))dnl
   m4_pushdef([AC_SUBST], [dn][l ])dnl
   m4_pushdef([AC_CONFIG_COMMANDS_PRE], [dn][l ])dnl
   m4_pushdef([acx_orig_macro],
     m4_bpatsubsts(m4_dquote(m4_defn([AC_FC_MODULE_FLAG])),
       [^      module conftest_module], [\&
      implicit none
      public],
       [^      use conftest_module], [\&, only : conftest_routine
      implicit none],
       [AC_LANG_P\(OP\|USH\)(\[?Fortran\]?)], [dn][l ],
       [FC_MODINC=.*], [dn][l ],
       [^ *#], [dn][l ]))dnl
   acx_orig_macro([:], [:])dnl
   m4_popdef([acx_orig_macro])dnl
   m4_popdef([AC_SUBST])dnl
   m4_popdef([AC_CONFIG_COMMANDS_PRE])dnl
   m4_popdef([AC_CACHE_CHECK])dnl
   m4_popdef([ac_cv_fc_module_flag])dnl
   AS_VAR_IF([acx_cv_fc_module_in_flag], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify search paths for module files])])], [$1])])

# ACX_FC_MODULE_OUT_FLAG([ACTION-IF-SUCCESS],
#                        [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to specify module output path.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_fc_module_out_flag variable, which may
# contain a significant trailing whitespace.
#
# The implementation patches the standard Autoconf macro
# AC_FC_MODULE_OUTPUT_FLAG to reduce the number of LANG switches and to avoid
# false negative results with the GFortran '-fmodule-private' flag.
#
AC_DEFUN([ACX_FC_MODULE_OUT_FLAG],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([ac_cv_fc_module_output_flag], [acx_cv_fc_module_out_flag])dnl
   m4_pushdef([AC_CACHE_CHECK],
     m4_bpatsubst(m4_dquote(m4_defn([AC_CACHE_CHECK])),
       [\$][1],
       [for Fortran compiler flag needed to specify output path for module dnl
files]))dnl
   m4_pushdef([AC_SUBST], [dn][l ])dnl
   m4_pushdef([AC_CONFIG_COMMANDS_PRE], [dn][l ])dnl
   m4_pushdef([acx_orig_macro],
     m4_bpatsubsts(m4_dquote(m4_defn([AC_FC_MODULE_OUTPUT_FLAG])),
       [^      module conftest_module], [\&
      implicit none
      public],
       [^      use conftest_module], [\&, only : conftest_routine
      implicit none],
       [AC_LANG_P\(OP\|USH\)(\[?Fortran\]?)], [dn][l ],
       [FC_MODOUT=.*], [dn][l ],
       [^ *#], [dn][l ]))dnl
   m4_version_prereq([2.70], [],
     [m4_define([acx_orig_macro],
        m4_bpatsubsts(m4_dquote(m4_defn([acx_orig_macro])),
          ['-mod '], ['-mdir ' \&],))])dnl
   acx_orig_macro([:], [:])dnl
   m4_popdef([acx_orig_macro])dnl
   m4_popdef([AC_SUBST])dnl
   m4_popdef([AC_CONFIG_COMMANDS_PRE])dnl
   m4_popdef([AC_CACHE_CHECK])dnl
   m4_popdef([ac_cv_fc_module_output_flag])dnl
   AS_VAR_IF([acx_cv_fc_module_out_flag], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify output path for module files])])], [$1])])

# ACX_FC_MODULE_NAMING([ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE])
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
   AS_IF([test "x$acx_cv_fc_module_naming_upper" = xunknown || dnl
test "x$acx_cv_fc_module_naming_ext" = xunknown],
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
