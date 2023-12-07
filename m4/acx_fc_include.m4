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

# ACX_FC_INCLUDE_FLAG([ACTION-IF-SUCCESS],
#                     [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to specify search paths for the Fortran
# "INCLUDE" statement. The result is either "unknown" or the actual compiler
# flag, which may contain a significant trailing whitespace.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_ftn_include_flag variable.
#
AC_DEFUN([ACX_FC_INCLUDE_FLAG],
  [_ACX_FC_INCLUDE_FLAG
   AS_VAR_IF([acx_cv_fc_ftn_include_flag], [unknown],
     [m4_default([$2],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify search paths for _ACX_FC_INCLUDE_DESC([ftn])])])], [$1])])

# ACX_FC_INCLUDE_FLAG_PP([ACTION-IF-SUCCESS],
#                        [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to specify search paths for the quoted form of
# the preprocessor "#include" directive. The result is either "unknown" or the
# actual compiler flag, which may contain a significant trailing whitespace.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_pp_include_flag variable.
#
AC_DEFUN([ACX_FC_INCLUDE_FLAG_PP],
  [_ACX_FC_INCLUDE_FLAG_PP
   AS_VAR_IF([acx_cv_fc_pp_include_flag], [unknown],
     [m4_default([$2],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify search paths for _ACX_FC_INCLUDE_DESC([pp])])])], [$1])])

# ACX_FC_INCLUDE_FLAG_PP_SYS([ACTION-IF-SUCCESS],
#                            [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to specify search paths for the angle-bracket
# form of the preprocessor "#include" directive. The result is either "unknown"
# or the actual compiler flag, which may contain a significant trailing
# whitespace.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_pp_sys_include_flag variable.
#
AC_DEFUN([ACX_FC_INCLUDE_FLAG_PP_SYS],
  [_ACX_FC_INCLUDE_FLAG_PP_SYS
   AS_VAR_IF([acx_cv_fc_pp_sys_include_flag], [unknown],
     [m4_default([$2],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
specify search paths for _ACX_FC_INCLUDE_DESC([pp_sys])])])], [$1])])

# ACX_FC_INCLUDE_ORDER([ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the search path order for the Fortran "INCLUDE" statement. The result
# is either "unknown" (e.g. in the case of cross-compilation) or a
# comma-separated list of identifiers that denote directories, in which the
# compiler searches for an included file:
#   "cwd": current working directory;
#   "flg": directories specified with the search path flags;
#   "src": directory containing the compiled source file;
#   "inc": directory containing the file with the the include statement or
#          directive.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_ftn_include_order variable.
#
AC_DEFUN([ACX_FC_INCLUDE_ORDER],
  [AC_REQUIRE([_ACX_FC_INCLUDE_FLAG])_ACX_FC_INCLUDE_ORDER([ftn],$@)])

# ACX_FC_INCLUDE_ORDER_PP([ACTION-IF-SUCCESS],
#                         [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the search path order for the quoted form of the preprocessor
# "#include" directive. See ACX_FC_INCLUDE_ORDER for the description of the
# result.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_pp_include_order variable.
#
AC_DEFUN([ACX_FC_INCLUDE_ORDER_PP],
  [AC_REQUIRE([_ACX_FC_INCLUDE_FLAG_PP])_ACX_FC_INCLUDE_ORDER([pp],$@)])

# ACX_FC_INCLUDE_ORDER_PP_SYS([ACTION-IF-SUCCESS],
#                             [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the search path order for the angle-bracket form of the preprocessor
# "#include" directive. See ACX_FC_INCLUDE_ORDER for the description of the
# result.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_pp_sys_include_order variable.
#
AC_DEFUN([ACX_FC_INCLUDE_ORDER_PP_SYS],
  [AC_REQUIRE([_ACX_FC_INCLUDE_FLAG_PP_SYS])dnl
_ACX_FC_INCLUDE_ORDER([pp_sys],$@)])

# ACX_FC_INCLUDE_CHECK(HEADER-FILE,
#                      [ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the header HEADER-FILE included with the Fortran "INCLUDE"
# statement is available. The result is either "yes" or "no".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the
# acx_cv_fc_ftn_header_[]AS_TR_SH(HEADER-FILE) variable.
#
AC_DEFUN([ACX_FC_INCLUDE_CHECK], [_ACX_FC_INCLUDE_CHECK([ftn],$@)])

# _ACX_FC_INCLUDE_LINE(HEADER-TYPE, HEADER-FILE)
# -----------------------------------------------------------------------------
# Expands into a line of code with the include statement or directive for the
# HEADER-FILE. The HEADER-TYPE defines, which actual statement or directive is
# used, and can be one of the following:
#   "ftn"    for the Fortran "INCLUDE" statement;
#   "pp"     for the quoted form of the preprocessor "#include" directive;
#   "pp_sys" for the angle-bracket form of the preprocessor "#include"
#            directive.
#
m4_define([_ACX_FC_INCLUDE_LINE],
  [m4_case([$1],
     [ftn], [m4_n([[      include "$2"]])],
     [pp], [m4_n([[@%:@include "$2"]])],
     [pp_sys], [m4_n([[@%:@include <$2>]])],
     [m4_fatal([unexpected header type: '$1'])])])

# _ACX_FC_INCLUDE_DESC(HEADER-TYPE)
# -----------------------------------------------------------------------------
# Expands into a shell string with the description of the HEADER-TYPE (see
# ACX_FC_INCLUDE_LINE).
#
m4_define([_ACX_FC_INCLUDE_DESC],
  [m4_case([$1],
     [ftn], [[the \"INCLUDE\" statement]],
     [pp], [[the quoted form of the \"#include\" directive]],
     [pp_sys], [[the angle-bracket form of the \"#include\" directive]],
     [m4_fatal([unexpected header type: '$1'])])])

# _ACX_FC_INCLUDE_KNOWN_FLAGS(HEADER-TYPE)
# -----------------------------------------------------------------------------
# Expands into a space-separated list of known flags needed to specify search
# paths for the HEADER-TYPE (see _ACX_FC_INCLUDE_LINE).
#
m4_define([_ACX_FC_INCLUDE_KNOWN_FLAGS],
  [m4_bmatch([$1],
     [ftn], [[-I '-I ']],
     [[pp|pp_sys]], [[-I '-I ' '-WF,-I' '-Wp,-I']],
     [m4_fatal([unexpected header type: '$1'])])])

# _ACX_FC_INCLUDE_FLAG()
# -----------------------------------------------------------------------------
# A parameterless alias for __ACX_FC_INCLUDE_FLAG([ftn]) to be used as an
# argument for AC_REQUIRE.
#
AC_DEFUN([_ACX_FC_INCLUDE_FLAG], [__ACX_FC_INCLUDE_FLAG([ftn])])

# _ACX_FC_INCLUDE_FLAG_PP()
# -----------------------------------------------------------------------------
# A parameterless alias for __ACX_FC_INCLUDE_FLAG([pp]) to be used as an
# argument for AC_REQUIRE.
#
AC_DEFUN([_ACX_FC_INCLUDE_FLAG_PP], [__ACX_FC_INCLUDE_FLAG([pp])])

# _ACX_FC_INCLUDE_FLAG_PP_SYS()
# -----------------------------------------------------------------------------
# A parameterless alias for __ACX_FC_INCLUDE_FLAG([pp_sys]) to be used as an
# argument for AC_REQUIRE.
#
AC_DEFUN([_ACX_FC_INCLUDE_FLAG_PP_SYS], [__ACX_FC_INCLUDE_FLAG([pp_sys])])

# __ACX_FC_INCLUDE_FLAG(HEADER-TYPE)
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to specify search paths for the HEADER-TYPE
# (see _ACX_FC_INCLUDE_LINE).
#
# The result is cached in the acx_cv_fc_[]HEADER-TYPE[]_include_flag variable.
#
# See _ACX_LANG_KNOWN_INC_FLAGS for the known flags.
#
m4_define([__ACX_FC_INCLUDE_FLAG],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_[]$1[]_include_flag])dnl
   AC_CACHE_CHECK([for Fortran compiler flag needed to specify search dnl
paths for _ACX_FC_INCLUDE_DESC([$1])], [acx_cache_var],
     [acx_cache_var=unknown
      AS_MKDIR_P([conftest.dir])
      AC_LANG_CONFTEST([AC_LANG_PROGRAM])
      mv conftest.$ac_ext conftest.dir/conftest.inc
      AC_LANG_CONFTEST([AC_LANG_SOURCE(
        [_ACX_FC_INCLUDE_LINE([$1], [conftest.inc])])])
      acx_save_FCFLAGS=$FCFLAGS
      for acx_flag in _ACX_FC_INCLUDE_KNOWN_FLAGS([$1]); do
        FCFLAGS="$acx_save_FCFLAGS ${acx_flag}conftest.dir"
        AC_LINK_IFELSE([], [AS_VAR_COPY([acx_cache_var], [acx_flag])])
        test "x$acx_cache_var" != xunknown && break
      done
      FCFLAGS=$acx_save_FCFLAGS
      rm -f conftest.$ac_ext
      rm -rf conftest.dir])
   m4_popdef([acx_cache_var])])

# _ACX_FC_INCLUDE_ORDER(HEADER-TYPE,
#                       [ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the search path order for the for the HEADER-TYPE (see
# _ACX_FC_INCLUDE_LINE). See ACX_FC_INCLUDE_ORDER for the description of the
# result.
#
# It is implied that __ACX_FC_INCLUDE_FLAG(HEADER-TYPE) has already been
# called.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_[]HEADER-TYPE[]_include_order
# variable.
#
m4_define([_ACX_FC_INCLUDE_ORDER],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_[]$1[]_include_order])dnl
   AC_CACHE_CHECK([for Fortran compiler search path order for dnl
_ACX_FC_INCLUDE_DESC([$1])], [acx_cache_var],
     [acx_cache_var=
      AS_VAR_IF([cross_compiling], [no],
        [AS_MKDIR_P([conftest.dir/src/inc])
         AS_MKDIR_P([conftest.dir/build])
         AS_MKDIR_P([conftest.dir/src/inc2])
         AC_LANG_CONFTEST([AC_LANG_PROGRAM([],
           [_ACX_FC_INCLUDE_LINE([$1], [conftest.inc])])])
dnl Copy the file to the build dir to keep _AC_MSG_LOG_CONFTEST happy.
dnl This copy does not get compiled.
         cp conftest.$ac_ext conftest.dir/build/conftest.$ac_ext
dnl This instance of the file will be compiled.
         mv conftest.$ac_ext conftest.dir/src/conftest.$ac_ext
         AC_LANG_CONFTEST([AC_LANG_SOURCE(
           [_ACX_FC_INCLUDE_LINE([$1], [conftest.write])])])
         mv conftest.$ac_ext conftest.dir/src/inc2/conftest.inc
         set "src" "/src/" "flg" "/src/inc/" "inc" "/src/inc2/" "cwd" "/build/"
         while test $[]@%:@ != 0; do
           AC_LANG_CONFTEST([AC_LANG_SOURCE([[      write(*,"(a)") "${1}"]])])
           shift; mv conftest.$ac_ext conftest.dir${1}conftest.write; shift
         done
         cd conftest.dir/build
         acx_save_FCFLAGS=$FCFLAGS
         FCFLAGS="$FCFLAGS ${acx_cv_fc_[]$1[]_include_flag}../src/inc dnl
${acx_cv_fc_[]$1[]_include_flag}../src/inc2"
         acx_save_ac_link=$ac_link
         ac_link=`AS_ECHO(["$ac_link"]) | sed 's%conftest\.\$ac_ext%../src/&%'`
         while :; do
           AC_LINK_IFELSE([],
             [acx_exec_result=`./conftest$ac_exeext`
              AS_IF([test $? -eq 0],
                [AS_CASE([$acx_exec_result],
                   [src], [rm -f ../src/conftest.write],
                   [inc], [rm -f ../src/inc2/conftest.write],
                   [cwd], [rm -f ./conftest.write],
                   [flg], [rm -f ../src/inc/conftest.write dnl
../src/inc2/conftest.write],
                   [break])
                 AS_IF([test -z "$acx_cache_var"],
                   [acx_cache_var=$acx_exec_result],
                   [AS_VAR_APPEND([acx_cache_var], [",$acx_exec_result"])])
                 rm -f conftest$ac_exeext],
                [break])],
             [break])
         done
         ac_link=$acx_save_ac_link
         FCFLAGS=$acx_save_FCFLAGS
         cd ../..
         rm -rf conftest.dir])
      AS_IF([test -z "$acx_cache_var"], [acx_cache_var=unknown])])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$3],
     [AC_MSG_FAILURE([unable to detect Fortran compiler search path dnl
order for _ACX_FC_INCLUDE_DESC([$1])])])], [$2])
   m4_popdef([acx_cache_var])])

# _ACX_FC_INCLUDE_CHECK(HEADER-TYPE,
#                       HEADER-FILE,
#                       [ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks the availability of the header file HEADER-FILE of the type
# HEADER-TYPE (see _ACX_FC_INCLUDE_LINE).
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the
# acx_cv_fc_[]HEADER-TYPE[]_header_[]AS_TR_SH(HEADER-FILE)
# variable.
#
m4_define([_ACX_FC_INCLUDE_CHECK],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_[]$1[]_header_[]AS_TR_SH([$2])])dnl
   AC_CACHE_CHECK([for $2], [acx_cache_var],
     [AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([],[_ACX_FC_INCLUDE_LINE([$1], [$2])])],
        [AS_VAR_SET([acx_cache_var], [yes])],
        [AS_VAR_SET([acx_cache_var], [no])])])
   AS_VAR_IF([acx_cache_var], [yes], [$3], [m4_default([$4],
     [AC_MSG_FAILURE([Fortran header file '$2' included with dnl
_ACX_FC_INCLUDE_DESC([$1]) is not available])])])
    m4_popdef([acx_cache_var])])
