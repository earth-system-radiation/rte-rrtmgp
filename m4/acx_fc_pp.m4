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

# ACX_FC_PP_SRCEXT(EXTENSION,
#                  [ACTION-IF-SUCCESS],
#                  [ACTION-IF-FAILURE = FAILURE])
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
# The implementation patches the standard Autoconf macro AC_FC_PP_SRCEXT to
# reduce the number of LANG switches and to support additional known compiler
# flags:
# Cray: -e T (must precede -e Z, which triggers generation of unwanted *.i
#             flags and crashes old versions of the compiler at the linking
#             stage)
#
AC_DEFUN([ACX_FC_PP_SRCEXT],
  [AC_LANG_ASSERT([Fortran])dnl
   acx_ext_save=$ac_ext
   m4_pushdef([acx_cache_var], [acx_cv_fc_pp_srcext_$1])dnl
   m4_pushdef([ac_cv_fc_pp_srcext_$1], [acx_cache_var])dnl
   m4_pushdef([AC_CACHE_CHECK],
     m4_bpatsubst(m4_dquote(m4_defn([AC_CACHE_CHECK])),
       [\$][1],
       [for Fortran compiler flag needed to compile preprocessed .$1 files]))dnl
   m4_pushdef([AC_FC_PP_SRCEXT],
     m4_bpatsubsts(m4_dquote(m4_defn([AC_FC_PP_SRCEXT])),
       ["-e Z"], ["-e T" \&],
       [AC_LANG_P\(OP\|USH\)(\[?Fortran\]?)], [dn][l ]))dnl
   AC_FC_PP_SRCEXT([$1], [], [:])
   m4_popdef([AC_FC_PP_SRCEXT])dnl
   m4_popdef([AC_CACHE_CHECK])dnl
   m4_popdef([ac_cv_fc_pp_srcext_$1])dnl
   AS_VAR_IF([acx_cache_var], [unknown],
     [ac_ext=$acx_ext_save
      m4_default([$3],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
compile preprocessed .$1 files])])],
     [ac_ext=$ac_fc_srcext
      $2])
   m4_popdef([acx_cache_var])])
