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
   AS_VAR_SET_IF([acx_cache_var],
     [AS_ECHO_N(["(cached) "]) >&AS_MESSAGE_FD],
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
