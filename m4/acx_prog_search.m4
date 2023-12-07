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

# ACX_PROG_SEARCH(VARIABLE,
#                 [CANDIDATES],
#                 [CHECK-SCRIPT = 'eval $acx_candidate'],
#                 [ACTION-IF-SUCCESS],
#                 [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Searches for the program (command) that results into a zero exit status of
# the CHECK-SCRIPT (defaults to running the candidate command). CHECK-SCRIPT
# can get the tested command from the shell variable $acx_candidate. If the
# shell variable VARIABLE is set, checks whether the value it stores passes the
# test. If VARIABLE is not set, iterates over the values of the blank-separated
# list CANDIDATES and stops when the first valid command is found. The value of
# VARIABLE is never set or changed.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# A positive result of this test is cached in the
# acx_cv_prog_[]AS_TR_SH(VARIABLE) variable.
#
AC_DEFUN([ACX_PROG_SEARCH],
  [AS_VAR_PUSHDEF([acx_cache_var], [acx_cv_prog_$1])dnl
   AS_LITERAL_IF([$1],
     [AC_MSG_CHECKING([for m4_tolower([$1])])],
     [acx_tmp=`AS_ECHO(["$1"]) | tr 'm4_cr_LETTERS' 'm4_cr_letters'`
      AC_MSG_CHECKING([for $acx_tmp])])
   AC_CACHE_VAL([acx_cache_var],
     [AS_VAR_SET_IF([$1], [set dummy "AS_VAR_GET([$1])"], [set dummy $2])
      shift
      for acx_candidate in "$[@]"; do
        m4_default([$3],
          [AC_TRY_COMMAND([$acx_candidate >&AS_MESSAGE_LOG_FD])])
        AS_IF([test $? -eq 0],
          [AS_VAR_SET([acx_cache_var], [$acx_candidate])
           break])
      done])
   AS_VAR_SET_IF([acx_cache_var],
     [AC_MSG_RESULT([AS_VAR_GET(acx_cache_var)])
      $4],
     [AC_MSG_RESULT([unknown])
      m4_default([$5],
        [AS_LITERAL_IF([$1],
           [AC_MSG_FAILURE([unable to find m4_tolower([$1])])],
           [acx_tmp=`AS_ECHO(["$1"]) | tr 'm4_cr_LETTERS' 'm4_cr_letters'`
            AC_MSG_FAILURE([unable to find $acx_tmp])])])])
   AS_VAR_POPDEF([acx_cache_var])])

# ACX_PROG_SEARCH_ABSPATH(PROG-TO-CHECK-FOR,
#                         [ACTION-IF-SUCCESS],
#                         [ACTION-IF-FAILURE = FAILURE],
#                         [PATH = $PATH])
# -----------------------------------------------------------------------------
# Searches for the absolute path to the PROG-TO-CHECK-FOR executable. If
# PROG-TO-CHECK-FOR contains slashes (i.e. specified as either absolute or
# relative path), the macro checks whether it is an executable and returns an
# absolute path to it (not necessarily in the canonical form). If the specified
# path is not a path to executable, the result is "unknown". If
# PROG-TO-CHECK-FOR does not contain slashes, the macro tries to find the
# executable in the list of directories stored in the PATH (defaults to the
# value of the $PATH shell variable). The value of the variable is interpreted
# as a list separated with the value of the $PATH_SEPARATOR shell variable,
# which is set by the configure script during the initialization (the usual
# value is ':'). If the path that contains PROG-TO-CHECK-FOR is a relative one,
# it will be converted to the absolute one. If PROG-TO-CHECK-FOR contains
# arguments, they will be preserved in the result.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is stored in the acx_prog_search_abspath shell variable.
#
AC_DEFUN([ACX_PROG_SEARCH_ABSPATH],
  [acx_prog_search_abspath=unknown
   set dummy $1; shift; acx_prog_exec=$[1]; shift; acx_prog_args="$[@]"
   AC_MSG_CHECKING([for the absolute path to $acx_prog_exec])
   AS_CASE([$acx_prog_exec],
     [*[[\\/]]*],
     [AS_IF([AS_EXECUTABLE_P([$acx_prog_exec])],
        [acx_prog_search_abspath=$acx_prog_exec])],
     [_AS_PATH_WALK([$4],
        [AS_IF([AS_EXECUTABLE_P(["$as_dir/$acx_prog_exec"])],
           [acx_prog_search_abspath="$as_dir/$acx_prog_exec"; break])])])
dnl If acx_prog_search_abspath is not "unknown", it is a path to an executable
dnl (without arguments).
   AS_CASE([$acx_prog_search_abspath],
     [unknown], [],
     [[[\\/]]* | ?:[[\\/]]*], [],
     [asx_dir=`echo "$acx_prog_search_abspath" | dnl
sed 's%/@<:@^/@:>@*$%%' 2>/dev/null`
      asx_file=`echo "$acx_prog_search_abspath" | sed 's%.*/%%' 2>/dev/null`
      asx_dir=`cd "$asx_dir" >/dev/null 2>&1 && pwd 2>/dev/null`
dnl Set the result to unknown until we make sure that we can provide a correct
dnl one.
      acx_prog_search_abspath=unknown
      AS_CASE([$asx_dir],
        [[[\\/]]* | ?:[[\\/]]*],
        [AS_IF([AS_EXECUTABLE_P(["$asx_dir/$asx_file"])],
           [acx_prog_search_abspath="$asx_dir/$asx_file"])])])
   AC_MSG_RESULT([$acx_prog_search_abspath])
   AS_VAR_IF([acx_prog_search_abspath], [unknown],
     [m4_default([$3],
        [AC_MSG_FAILURE(
           [unable to find the absolute path to $acx_prog_exec])])],
     [AS_IF([test -n "$acx_prog_args"],
        [AS_VAR_APPEND([acx_prog_search_abspath], [" $acx_prog_args"])])
      $2])])
