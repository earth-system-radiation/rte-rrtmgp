# ACX_FC_LINE_LENGTH([LENGTH = 132],
#                    [ACTION-IF-SUCCESS = APPEND-FCFLAGS],
#                    [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Originally taken from the master branch of Autoconf where it is known as
# AC_FC_LINE_LENGTH.
# -----------------------------------------------------------------------------
# Finds the Fortran compiler flag needed to accept lines of length LENGTH,
# where LENGTH may be 80, 132 (default), or 'unlimited' for longer lines. The
# result is either "unknown", or the actual compiler flag required to to accept
# lines of length LENGTH.
#
# If successful, runs ACTION-IF-SUCCESS (defaults to appending the result to
# FCFLAGS), otherwise runs ACTION-IF-FAILURE (defaults to failing with an error
# message).
#
# The flag is cached in the acx_cv_fc_line_length_[]LENGTH.
#
# Known flags:
# -f{free,fixed}-line-length-N with N 72, 80, 132, or 0 or none for none.
# gfortran: -ffree-line-length-none
# g95: -ffree-line-length-huge (also -ffixed-line-length-N as above)
# IBM: -qfixed=132 80 72
# Cray: -Mextend
# Intel: -132 -80 -72 (needs to come before -extend_source because ifort
#                      accepts that as well with an optional parameter and
#                      does not fail but only warns about unknown arguments)
# SGI: -extend_source
# Absoft: -W, -WNN (132, 80, 72)
# HP: +es, +extend_source (254 in either form, default is 72 fixed, 132 free)
# Lahey/Fujitsu: -w, (-)-wide (255 cols in fixed form)
# SUN: -e (132 characters)
# NAGWare: -132
# f2c: -72, -f, -Wf,-f: f2c (a weak form of "free-form" and long lines)
# Open Watcom: /XLine
#
AC_DEFUN([ACX_FC_LINE_LENGTH],
  [AC_LANG_ASSERT([Fortran])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_line_length_$1])dnl
   AC_MSG_CHECKING([for fortran flag needed to accept dnl
m4_default([$1], [132]) column source lines])
   AC_CACHE_VAL([acx_cache_var],
     [acx_cache_var=unknown
      m4_case(m4_default([$1], [132]),
        [unlimited],
        [acx_flag_suffix=0
         AC_LANG_CONFTEST([AC_LANG_SOURCE([[
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)
        implicit none
        integer :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19
      end subroutine]])])],
        [132],
        [acx_flag_suffix=132
         AC_LANG_CONFTEST([AC_LANG_SOURCE([[
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
        implicit none
        integer :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10
      end subroutine]])])],
        [80],
        [acx_flag_suffix=80
         AC_LANG_CONFTEST([AC_LANG_SOURCE([[
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
        implicit none
        integer :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9
      end subroutine]])])],
      [m4_fatal([Invalid LENGTH argument for ACX_FC_LINE_LENGTH: '$1'])])
      acx_save_FCFLAGS=$FCFLAGS
      for acx_flag in '' \
               -ffree-line-length-none -ffixed-line-length-none \
               -ffree-line-length-huge \
               -ffree-line-length-$acx_flag_suffix \
               -ffixed-line-length-$acx_flag_suffix \
               -qfixed=$acx_flag_suffix -Mextend \
               -$acx_flag_suffix -extend_source \
               -W$acx_flag_suffix -W +extend_source +es -wide --wide -w -e \
               -f -Wf,-f -xline; do
        FCFLAGS="$acx_save_FCFLAGS $acx_flag"
        AC_COMPILE_IFELSE([],[acx_cache_var=$acx_flag])
        test "x$acx_cache_var" != xunknown && break
      done
      rm -f conftest.$ac_ext
      FCFLAGS=$acx_save_FCFLAGS])
   AS_IF([test -n "$acx_cache_var"],
     [AC_MSG_RESULT([$acx_cache_var])],
     [AC_MSG_RESULT([none needed])])
   AS_VAR_IF([acx_cache_var], [unknown],
     [m4_default([$3],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
accept m4_default([$1], [132]) column source lines])])],
     [m4_default([$2],
        [AS_IF([test -n "$acx_cache_var"],
           [AS_VAR_APPEND([FCFLAGS], [" $acx_cache_var"])])])])
   m4_popdef([acx_cache_var])])
