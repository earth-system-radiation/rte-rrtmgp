# ACX_LANG_MACRO_FLAG([ACTION-IF-SUCCESS],
#                     [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is known as
# AC_FC_PP_DEFINE.
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to specify a preprocessor macro definition.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_macro_flag variable.
#
# See _ACX_LANG_KNOWN_MACRO_FLAGS for the known flags.
#
AC_DEFUN([ACX_LANG_MACRO_FLAG],
  [m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_macro_flag])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler flag needed to define a dnl
preprocessor macro],
     [acx_cache_var],
     [acx_cache_var=unknown
      AC_LANG_CONFTEST(
        [AC_LANG_PROGRAM([],
[[#ifndef CONFTEST_ONE
      choke me
#endif
#if CONFTEST_TWO != 42
      choke me
#endif]])])
      acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      for acx_lang_macro_flag in _ACX_LANG_KNOWN_MACRO_FLAGS; do
        _AC_LANG_PREFIX[]FLAGS="${acx_save_[]_AC_LANG_PREFIX[]FLAGS} dnl
${acx_lang_macro_flag}CONFTEST_ONE ${acx_lang_macro_flag}CONFTEST_TWO=42"
        AC_COMPILE_IFELSE([],
          [acx_cache_var=$acx_lang_macro_flag
           break])
      done
      rm -f conftest.$ac_ext
      _AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect _AC_LANG compiler flag needed to dnl
define a preprocessor macro])])], [$1])
   m4_popdef([acx_cache_var])])

# ACX_LANG_MACRO_CHECK_DEFINED(MACRO-NAME)
# -----------------------------------------------------------------------------
# Checks whether the preprocessor macro MACRO-NAME is defined. The result is
# either "yes", "no" or "unsupported" (if the current language does not
# support preprocessor directives).
#
# The result is stored in the acx_macro_defined variable and cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_macro_[]AS_TR_SH(MACRO-NAME)_defined variable.
#
AC_DEFUN([ACX_LANG_MACRO_CHECK_DEFINED],
  [m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_macro_[]AS_TR_SH([$1])_defined])dnl
   AC_CACHE_CHECK([whether the _AC_LANG preprocessor macro $1 is defined],
     [acx_cache_var],
     [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifdef $1
#else
      choke me
#endif]])],
        [AS_VAR_SET([acx_cache_var], [yes])],
        [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[#ifndef $1
#else
      choke me
#endif]])],
           [AS_VAR_SET([acx_cache_var], [no])],
           [AS_VAR_SET([acx_cache_var], [unsupported])])])])
   AS_VAR_COPY([acx_macro_defined], [acx_cache_var])
   m4_popdef([acx_cache_var])])

# ACX_LANG_MACRO_CHECK_VALUE(MACRO-NAME,
#                            [KNOWN-INTEGER-VALUES])
# -----------------------------------------------------------------------------
# Detects the value of the preprocessor macro MACRO-NAME. First, tries to link
# and to run a program that prints the value of MACRO-NAME. If that is
# successful, returns the output of the program. Otherwise (e.g. in the case of
# cross-compilation), goes through the optionally provided space-separated list
# of integers KNOWN-INTEGER-VALUES and checks whether MACRO-NAME expands to one
# of them. The result is either "unknown" or the actual value of the macro.
#
# The result is stored in the acx_macro_value variable and cached in the
# acx_cv_[]_AC_LANG_ABBREV[]_macro_[]AS_TR_SH(MACRO-NAME)_value variable.
#
AC_DEFUN([ACX_LANG_MACRO_CHECK_VALUE],
  [m4_pushdef([acx_cache_var],
     [acx_cv_[]_AC_LANG_ABBREV[]_macro_[]AS_TR_SH([$1])_value])dnl
   AC_CACHE_CHECK([for the value of the _AC_LANG preprocessor macro $1],
     [acx_cache_var],
     [AS_VAR_SET([acx_cache_var], [unknown])
      AS_VAR_IF([cross_compiling], [no],
        [AC_LINK_IFELSE([_ACX_LANG_MACRO_PRINT_PROGRAM([$1])],
           [acx_exec_result=`./conftest$ac_exeext 2>&AS_MESSAGE_LOG_FD`
            AS_IF([test $? -eq 0],
              [AS_VAR_COPY([acx_cache_var], [acx_exec_result])])])])
      m4_ifnblank([$2],
        [AS_VAR_IF([acx_cache_var], [unknown],
           [set dummy $2; shift
            while test $[]@%:@ != 0; do
              AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
[[#if $1 == _CONFTEST_UNDEFINED_OR_EMPTY || $1 != $][1
      choke me
#endif]])],
                [AS_VAR_COPY([acx_cache_var], [1])
                 set dummy; shift],
                [shift])
            done])])])
   AS_VAR_COPY([acx_macro_value], [acx_cache_var])
   m4_popdef([acx_cache_var])])

# _ACX_LANG_KNOWN_MACRO_FLAGS()
# -----------------------------------------------------------------------------
# Expands into a language-specific space-separated list of known flags needed
# to specify a preprocessor macro definition. By default, expands to m4_fatal
# with the message saying that _AC_LANG is not supported.
#
m4_define([_ACX_LANG_KNOWN_MACRO_FLAGS],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_fatal([the list of known ]_AC_LANG[ compiler flags needed to ]dnl
[specify a preprocessor macro definition is undefined])])])

# _ACX_LANG_KNOWN_MACRO_FLAGS(C)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_KNOWN_MACRO_FLAGS for C language.
#
m4_define([_ACX_LANG_KNOWN_MACRO_FLAGS(C)], [-D])

# _ACX_LANG_KNOWN_MACRO_FLAGS(Fortran)(HEADER-TYPE)
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_KNOWN_MACRO_FLAGS for Fortran language.
#
# Known flags:
# IBM: -WF,-D
# Lahey/Fujitsu: -Wp,-D     older versions???
# f2c: -D or -Wc,-D
# others: -D
#
m4_define([_ACX_LANG_KNOWN_MACRO_FLAGS(Fortran)], [-D -WF,-D -Wp,-D -Wc,-D])

# _ACX_LANG_MACRO_PRINT_PROGRAM(MACRO-NAME)
# -----------------------------------------------------------------------------
# Expands into the source code of a program in the current language that prints
# the value of the preprocessor macro MACRO-NAME. The program fails if
# MACRO-NAME is not defined. By default, expands to m4_fatal with the message
# saying that _AC_LANG is not supported.
#
m4_define([_ACX_LANG_MACRO_PRINT_PROGRAM],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_fatal([the macro print program is not defined for ]dnl
_AC_LANG[ language])])])

# _ACX_LANG_MACRO_PRINT_PROGRAM(C)(MACRO-NAME)
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_MACRO_PRINT_PROGRAM for C language.
#
m4_define([_ACX_LANG_MACRO_PRINT_PROGRAM(C)],
  [AC_LANG_PROGRAM([[#include <stdio.h>]],
[[#ifndef $1
choke me
#else
#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
printf("%s\n", STRINGIFY($1));
#endif]])])

# _ACX_LANG_MACRO_PRINT_PROGRAM(C)(MACRO-NAME)
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_MACRO_PRINT_PROGRAM for C++ language.
#
m4_copy([_ACX_LANG_MACRO_PRINT_PROGRAM(C)],
  [_ACX_LANG_MACRO_PRINT_PROGRAM(C++)])

# _ACX_LANG_MACRO_PRINT_PROGRAM(Fortran)(MACRO-NAME)
# -----------------------------------------------------------------------------
# Implementation of _ACX_LANG_MACRO_PRINT_PROGRAM for Fortran language. The
# program compilation succeeds only if MACRO-NAME expands either to an integer,
# to an empty string or to a quoted string. If MACRO-NAME expands to a quoted
# string, the output of the program is always quoted with the double quotation
# marks (""), even if the actual value of MACRO-NAME is a string quoted with
# the single quotation marks ('').
#
m4_define([_ACX_LANG_MACRO_PRINT_PROGRAM(Fortran)],
  [AC_LANG_SOURCE(
[[#ifndef $1
      choke me
#else
      subroutine p_str(s)
      implicit none
      character(len=*) :: s
      write(*, "(a)") """"//s//""""
      end subroutine

      subroutine p_int(i)
      implicit none
      integer :: i
      write(*, "(i0)") i
      end subroutine

      subroutine p_none()
      write(*, "(a)") ""
      end subroutine

      program main
      implicit none
      interface p_val
      subroutine p_str(s)
      implicit none
      character(len=*) :: s
      end subroutine
      subroutine p_int(i)
      implicit none
      integer :: i
      end subroutine
      subroutine p_none()
      end subroutine
      end interface
      call p_val($1)
      end
#endif]])])
