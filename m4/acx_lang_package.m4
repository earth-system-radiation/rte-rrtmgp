# Copyright (c) 2010-2024, DKRZ, MPI-M
#
# Authors:
#   Thomas Jahns <jahns@dkrz.de>
#   Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
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

# ACX_LANG_PACKAGE_INIT(PACKAGE-NAME,
#                       [INC-SEARCH-FLAGS],
#                       [LIB-SEARCH-FLAGS],
#                       [INC-SEARCH-SUFFIX = /include],
#                       [LIB-SEARCH-SUFFIX = /lib])
# -----------------------------------------------------------------------------
# Sets command-line arguments of the configure script that allows for setting
# search paths for the PACKAGE-NAME. By default, sets only the
# "--with-package-name-root" argument.
#
# If the argument INC-SEARCH-FLAGS is not blank, adds a command-line argument
# "--with-package-name-include", declares a precious variable
# AS_TR_CPP(PACKAGE-NAME)_AC_LANG_PREFIX[]FLAGS (e.g. PACKAGE_NAME_FCFLAGS),
# and sets a shell variable
# acx_[]_AC_LANG_ABBREV[]_[]AS_TR_SH(PACKAGE-NAME)_inc_search_args
# (e.g. acx_fc_Package_Name_inc_search_args). The latter variable is set to
# a string containing each flag from the space-separated list INC-SEARCH-FLAGS
# appended either with the value of the command-line argument
# "--with-package-name-include" (if given) or with the concatenation of the
# value of the command-line argument "--with-package-name-root" (if given) and
# INC-SEARCH-SUFFIX (defaults to /include). If neither of the two mentioned
# command-line arguments is given the variable
# acx_[]_AC_LANG_ABBREV[]_[]AS_TR_SH(PACKAGE-NAME)_inc_search_args is empty.
#
# If the argument LIB-SEARCH-FLAGS is given, adds a command-line argument
# "--with-package-name-lib", declares a precious variable
# AS_TR_CPP(PACKAGE-NAME)_AC_LANG_PREFIX[]LIBS (e.g. PACKAGE_NAME_FCLIBS), and
# sets a shell variable
# acx_[]_AC_LANG_ABBREV[]_[]AS_TR_SH(PACKAGE-NAME)_lib_search_args (e.g.
# acx_fc_Package_Name_lib_search_args). The value for the latter variable is
# set according to the same rules as for the variable related to the include
# search path (LIB-SEARCH-SUFFIX defaults to /lib).
#
AC_DEFUN([ACX_LANG_PACKAGE_INIT],
  [m4_pushdef([acx_package_ROOT], [AS_TR_CPP([$1]_ROOT)])dnl
   m4_pushdef([acx_package_with_root],
     [with_[]AS_TR_SH([ASX_TR_ARG([$1])])_root])dnl
   AC_ARG_WITH(ASX_TR_ARG([$1])[-root],
     [AS_HELP_STRING([--with-ASX_TR_ARG([$1])-root=[]acx_package_ROOT],
        [root search path for $1 headers and libraries])])
   m4_ifnblank([$2],
     [m4_pushdef([acx_package_with_include],
        [with_[]AS_TR_SH([ASX_TR_ARG([$1])])_include])dnl
      AC_ARG_WITH(ASX_TR_ARG([$1])[-include],
        [AS_HELP_STRING([--with-ASX_TR_ARG([$1])-include=DIR],
           [search path for $1 headers @<:@]acx_package_ROOT[]dnl
m4_default([$4], [/include])[@:>@])], [],
           [AS_VAR_SET_IF([acx_package_with_root],
              [acx_package_with_include="${acx_package_with_root}dnl
m4_default([$4], [/include])"])])
      AC_ARG_VAR(AS_TR_CPP([$1])_[]_AC_LANG_PREFIX[FLAGS],
        [exact ]_AC_LANG[ compiler flags enabling $1])
      AS_VAR_SET_IF([acx_package_with_include],
        [for acx_flag in $2; do
           ASX_VAR_APPEND_UNIQ(
             [acx_[]_AC_LANG_ABBREV[]_[]AS_TR_SH([$1])_inc_search_args],
             ["${acx_flag}${acx_package_with_include}"], [" "])
         done])
      m4_popdef([acx_package_with_include])])
   m4_ifnblank([$3],
     [m4_pushdef([acx_package_with_lib],
        [with_[]AS_TR_SH([ASX_TR_ARG([$1])])_lib])dnl
      AC_ARG_WITH(ASX_TR_ARG([$1])[-lib],
        [AS_HELP_STRING([--with-ASX_TR_ARG([$1])-lib=DIR],
           [search path for $1 libraries @<:@]acx_package_ROOT[]dnl
m4_default([$5], [/lib])[@:>@])],
           [],
           [AS_VAR_SET_IF([acx_package_with_root],
              [acx_package_with_lib="${acx_package_with_root}dnl
m4_default([$5], [/lib])"])])
      AC_ARG_VAR(AS_TR_CPP([$1])_[]_AC_LANG_PREFIX[LIBS],
        [exact linker flags enabling $1 when linking with ]_AC_LANG[ compiler])
      AS_VAR_SET_IF([acx_package_with_lib],
        [for acx_flag in $3; do
           ASX_VAR_APPEND_UNIQ(
             [acx_[]_AC_LANG_ABBREV[]_[]AS_TR_SH([$1])_lib_search_args],
             ["${acx_flag}${acx_package_with_lib}"], [" "])
         done])
      m4_popdef([acx_package_with_lib])])
   m4_popdef([acx_package_ROOT])dnl
   m4_popdef([acx_package_with_root])])
