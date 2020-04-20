# ACX_LANG_LIB_CHECK(FUNC-NAME,
#                    [ACTION-IF-SUCCESS],
#                    [ACTION-IF-FAILURE = FAILURE],
#                    [CHECK-PROGRAM = AC_LANG_CALL([],FUNC-NAME)])
# -----------------------------------------------------------------------------
# Checks whether the function FUNC-NAME is available for the current language.
# The check is performed by linking a program CHECK-PROGRAM (defaults to
# AC_LANG_CALL([], FUNC-NAME)). The result is either "yes" or "no".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_func_[]AS_TR_SH(FUNC-NAME) variable if the current
# language is case-sensitive (e.g. C), or in the
# acx_cv_[]_AC_LANG_ABBREV[]_func_[]AS_TR_CPP(FUNC-NAME) variable if the
# current language is case-insensitive (e.g. Fortran).
#
AC_DEFUN([ACX_LANG_LIB_CHECK],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_indir([$0()], $@)])])

# ACX_LANG_LIB_SEARCH(VARIABLE,
#                     FUNC-NAME,
#                     [CANDIDATES],
#                     [ACTION-IF-SUCCESS],
#                     [ACTION-IF-FAILURE = FAILURE],
#                     [CHECK-PROGRAM = AC_LANG_CALL([],FUNC-NAME)])
# -----------------------------------------------------------------------------
# Searches for a set of linker flags enabling the function FUNC-NAME for the
# current language. If the shell variable VARIABLE is set, checks whether the
# linker flags it stores enable FUNC-NAME. If VARIABLE is not set, checks
# whether additional flags are needed at all to enable FUNC-NAME. If the latter
# is the case, iterates over the values of the blank-separated list CANDIDATES
# and stops when the first value corresponding to the valid set of flags
# enabling FUNC-NAME is found. The checks are performed by trying to compile
# and link the program CHECK-PROGRAM (defaults to AC_LANG_CALL for FUNC-NAME).
# The result of the macro is either an empty string (i.e. no additional flags
# are needed), or the first successful element of the list CANDIDATES. The
# value of VARIABLE is never set or changed.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# A positive result of this test is cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_lib_func_[]AS_TR_SH(FUNC-NAME) variable if the
# current language is case-sensitive (e.g. C), or in the
# acx_cv_[]_AC_LANG_ABBREV[]_lib_func_[]AS_TR_CPP(FUNC-NAME) variable if the
# current language is case-insensitive (e.g. Fortran).
#
AC_DEFUN([ACX_LANG_LIB_SEARCH],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_indir([$0()], $@)])])

# ACX_LANG_LIB_CHECK()(FUNC-NAME,
#                      [ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE],
#                      [CHECK-PROGRAM = AC_LANG_CALL([],FUNC-NAME)])
# -----------------------------------------------------------------------------
# Default implementation of ACX_LANG_LIB_CHECK for case-sensitive languages.
#
# The result is cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_func_[]AS_TR_SH(FUNC-NAME) variable.
#
m4_define([ACX_LANG_LIB_CHECK()],
  [m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_func_[]AS_TR_SH([$1])])dnl
   AC_CACHE_CHECK([for _AC_LANG function $1],
     [acx_cache_var],
     [AC_LINK_IFELSE([m4_default([$4], [AC_LANG_CALL([], [$1])])],
        [AS_VAR_SET([acx_cache_var], [yes])],
        [AS_VAR_SET([acx_cache_var], [no])])])
   AS_IF([test x"AS_VAR_GET(acx_cache_var)" = xyes], [$2],
     [m4_default([$3],
        [AC_MSG_FAILURE([_AC_LANG function $1 is not available])])])
   m4_popdef([acx_cache_var])])

# ACX_LANG_LIB_CHECK(Fortran)(FUNC-NAME,
#                             [ACTION-IF-SUCCESS],
#                             [ACTION-IF-FAILURE = FAILURE],
#                             [CHECK-PROGRAM = AC_LANG_CALL([],FUNC-NAME)])
# -----------------------------------------------------------------------------
# Implementation of ACX_LANG_LIB_CHECK for Fortran language. Accounts for the
# case-insensitivity of the language.
#
# The result is cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_func_[]AS_TR_CPP(FUNC-NAME) variable.
#
m4_define([ACX_LANG_LIB_CHECK(Fortran)],
  [AS_VAR_SET([acx_tmp], [AS_TR_CPP([$1])])
   m4_indir([ACX_LANG_LIB_CHECK()], [$acx_tmp], m4_shift($@))])

# ACX_LANG_LIB_SEARCH()(VARIABLE,
#                       FUNC-NAME,
#                       [CANDIDATES],
#                       [ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE],
#                       [CHECK-PROGRAM = AC_LANG_CALL([],FUNC-NAME)])
# -----------------------------------------------------------------------------
# Default implementation of ACX_LANG_LIB_SEARCH for case-sensitive languages.
#
# The result is cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_lib_func_[]AS_TR_SH(FUNC-NAME) variable.
#
m4_define([ACX_LANG_LIB_SEARCH()],
  [m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_lib_func_[]AS_TR_SH([$2])])dnl
   AC_MSG_CHECKING([for linker flags enabling _AC_LANG function $2])
   AC_CACHE_VAL([acx_cache_var],
     [AC_LANG_CONFTEST([m4_default([$6], [AC_LANG_CALL([], [$2])])])
      acx_save_LIBS=$LIBS
      AS_VAR_SET_IF([$1], [set dummy "AS_VAR_GET([$1])"], [set dummy '' $3])
      shift
      for acx_libs in "$[@]"; do
        LIBS="$acx_libs $acx_save_LIBS"
        AC_LINK_IFELSE([], [AS_VAR_COPY([acx_cache_var], [acx_libs])])
        AS_VAR_SET_IF([acx_cache_var], [break])
      done
      rm -f conftest.$ac_ext
      LIBS=$acx_save_LIBS])
   AS_VAR_SET_IF([acx_cache_var],
     [AS_IF([test -n "AS_VAR_GET(acx_cache_var)"],
        [AC_MSG_RESULT([AS_VAR_GET(acx_cache_var)])],
        [AC_MSG_RESULT([none needed])])
      $4],
     [AC_MSG_RESULT([unknown])
      m4_default([$5], [AC_MSG_FAILURE(
        [unable to find linker flags enabling _AC_LANG function $2])])])
   m4_popdef([acx_cache_var])])

# ACX_LANG_LIB_SEARCH(Fortran)(VARIABLE,
#                              FUNC-NAME,
#                              [CANDIDATES],
#                              [ACTION-IF-SUCCESS],
#                              [ACTION-IF-FAILURE = FAILURE],
#                              [CHECK-PROGRAM = AC_LANG_CALL([],FUNC-NAME)])
# -----------------------------------------------------------------------------
# Implementation of ACX_LANG_LIB_SEARCH for Fortran language. Accounts for the
# case-insensitivity of the language.
#
# The result is cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_lib_func_[]AS_TR_CPP(FUNC-NAME) variable.
#
m4_define([ACX_LANG_LIB_SEARCH(Fortran)],
  [AS_VAR_SET([acx_fc_lib_search], [AS_TR_CPP([$2])])
   m4_indir([ACX_LANG_LIB_SEARCH()],
     [$1], [$acx_fc_lib_search], m4_shift2($@))])
