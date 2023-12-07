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

# ACX_SHLIB_FC_RPATH_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the Fortran compiler flag needed to add a directory to the
# runtime library search path.
#
# The result is cached in the acx_cv_fc_rpath_flag variable.
#
AC_DEFUN([ACX_SHLIB_FC_RPATH_FLAG],
  [AC_REQUIRE([ACX_COMPILER_FC_VENDOR])_ACX_SHLIB_RPATH_FLAG])

# ACX_SHLIB_CC_RPATH_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the C compiler flag needed to add a directory to the
# runtime library search path.
#
# The result is cached in the acx_cv_c_rpath_flag variable.
#
AC_DEFUN([ACX_SHLIB_CC_RPATH_FLAG],
  [AC_REQUIRE([ACX_COMPILER_CC_VENDOR])_ACX_SHLIB_RPATH_FLAG])

# ACX_SHLIB_CXX_RPATH_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the C++ compiler flag needed to add a directory to the
# runtime library search path.
#
# The result is cached in the acx_cv_cxx_rpath_flag variable.
#
AC_DEFUN([ACX_SHLIB_CXX_RPATH_FLAG],
  [AC_REQUIRE([ACX_COMPILER_CXX_VENDOR])_ACX_SHLIB_RPATH_FLAG])

# ACX_SHLIB_RPATH_FLAGS_CHECK([RPATH-FLAGS],
#                             [ACTION-IF-SUCCESS],
#                             [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Expands to a shell script that checks whether the current compiler accepts
# the automatically generated RPATH flags RPATH-FLAGS by trying to link a dummy
# program with LDFLAGS set to "RPATH-FLAGS $LDFLAGS".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
AC_DEFUN([ACX_SHLIB_RPATH_FLAGS_CHECK],
  [acx_shlib_rpath_flags_check_result=no
   AC_MSG_CHECKING([whether _AC_LANG compiler accepts the automatically dnl
generated RPATH flags])
   m4_ifnblank([$1],
     [acx_save_LDFLAGS=$LDFLAGS
      LDFLAGS="$1 $LDFLAGS"])
   AC_LINK_IFELSE([AC_LANG_PROGRAM],
     [acx_shlib_rpath_flags_check_result=yes])
   m4_ifnblank([$1], [LDFLAGS=$acx_save_LDFLAGS])
   AC_MSG_RESULT([$acx_shlib_rpath_flags_check_result])
   AS_VAR_IF([acx_shlib_rpath_flags_check_result], [yes], [$2],
     [m4_default([$3], [AC_MSG_FAILURE([_AC_LANG compiler does not accept dnl
the automatically generated RPATH flags[]m4_ifnblank([$1],[ '$1'])])])])])

# ACX_SHLIB_FC_PIC_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the Fortran compiler flag needed to generate the position
# independent code (PIC).
#
# The result is cached in the acx_cv_fc_pic_flag variable.
#
AC_DEFUN([ACX_SHLIB_FC_PIC_FLAG],
  [AC_REQUIRE([ACX_COMPILER_FC_VENDOR])_ACX_SHLIB_PIC_FLAG])

# ACX_SHLIB_CC_PIC_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the C compiler flag needed to generate the position
# independent code (PIC).
#
# The result is cached in the acx_cv_c_pic_flag variable.
#
AC_DEFUN([ACX_SHLIB_CC_PIC_FLAG],
  [AC_REQUIRE([ACX_COMPILER_CC_VENDOR])_ACX_SHLIB_PIC_FLAG([-DPIC])])

# ACX_SHLIB_CXX_PIC_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the C++ compiler flag needed to generate the position
# independent code (PIC).
#
# The result is cached in the acx_cv_cxx_pic_flag variable.
#
AC_DEFUN([ACX_SHLIB_CXX_PIC_FLAG],
  [AC_REQUIRE([ACX_COMPILER_CXX_VENDOR])_ACX_SHLIB_PIC_FLAG([-DPIC])])

# ACX_SHLIB_PIC_FLAGS_CHECK([PIC-FLAGS],
#                           [ACTION-IF-SUCCESS],
#                           [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Expands to a shell script that checks whether the current compiler accepts
# the automatically generated PIC flags PIC-FLAGS by trying to link a dummy
# program with the compiler-specific flags appended with PIC-FLAGS.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
AC_DEFUN([ACX_SHLIB_PIC_FLAGS_CHECK],
  [acx_shlib_pic_flags_check_result=no
   AC_MSG_CHECKING([whether _AC_LANG compiler accepts the automatically dnl
generated PIC flags])
   m4_ifnblank([$1],
     [acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      AS_VAR_APPEND([_AC_LANG_PREFIX[]FLAGS], [" $1"])])
   AC_LINK_IFELSE([AC_LANG_PROGRAM],
     [acx_shlib_pic_flags_check_result=yes])
   m4_ifnblank([$1],
     [_AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS])
   AC_MSG_RESULT([$acx_shlib_pic_flags_check_result])
   AS_VAR_IF([acx_shlib_pic_flags_check_result], [yes], [$2],
     [m4_default([$3], [AC_MSG_FAILURE([_AC_LANG compiler does not accept dnl
the automatically generated PIC flags[]m4_ifnblank([$1],[ '$1'])])])])])

# ACX_SHLIB_PATH_VAR()
# -----------------------------------------------------------------------------
# Sets the result to the name of the environment variable specifying the search
# paths for shared libraries.
#
# The result is cached in the acx_cv_shlib_path_var variable.
#
AC_DEFUN([ACX_SHLIB_PATH_VAR],
  [AC_REQUIRE([AC_CANONICAL_HOST])dnl
   AC_CACHE_CHECK([for the name of the environment variable specifying the dnl
search paths for shared libraries], [acx_cv_shlib_path_var],
     [AS_CASE([$host_os],
        [darwin*], [acx_cv_shlib_path_var=DYLD_LIBRARY_PATH],
        [acx_cv_shlib_path_var=LD_LIBRARY_PATH])])])

# _ACX_SHLIB_RPATH_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the compiler flag needed to add a directory to the runtime
# library search path (requires calling _ACX_COMPILER_VENDOR first).
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_rpath_flag variable.
#
m4_define([_ACX_SHLIB_RPATH_FLAG],
  [m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_rpath_flag])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler flag needed to add a directory to dnl
the runtime library search path], [acx_cache_var],
     [AS_CASE([AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])],
        [nag], [acx_cache_var="-Wl,-Wl,,-rpath -Wl,-Wl,,"],
        [acx_cache_var="-Wl,-rpath -Wl,"])])
   m4_popdef([acx_cache_var])])

# _ACX_SHLIB_PIC_FLAG([COMMON-EXTRA-FLAG])
# -----------------------------------------------------------------------------
# Sets the result to the compiler flag needed to generate the position
# independent code (PIC) (requires calling _ACX_COMPILER_VENDOR first). When
# provided, COMMON-EXTRA-FLAG is appended to the result.
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_pic_flag variable.
#
m4_define([_ACX_SHLIB_PIC_FLAG],
  [m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_pic_flag])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler flag needed to produce PIC],
     [acx_cache_var],
     [AS_CASE([AS_VAR_GET([acx_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])],
        [nag], [acx_cache_var='-PIC'],
        [portland], [acx_cache_var='-fpic'],
        [sun], [acx_cache_var='-KPIC'],
        [ibm], [acx_cache_var='-qpic'],
        [acx_cache_var='-fPIC'])
      m4_ifnblank([$1], [AS_VAR_APPEND([acx_cache_var], [" $1"])])])
   m4_popdef([acx_cache_var])])
