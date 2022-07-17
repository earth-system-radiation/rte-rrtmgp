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
# The implementation is based on the standard Autoconf macro AC_FC_PP_SRCEXT,
# which is monkey-patched to reduce the number of LANG switches and to support
# additional known compiler flags:
# Cray: -e T (must precede -e Z, which triggers generation of unwanted *.i flags
#             and crashes old versions of the compiler at the linking stage)
#
AC_DEFUN([ACX_FC_PP_SRCEXT],
  [AC_LANG_ASSERT([Fortran])dnl
   acx_ext_save=$ac_ext
   m4_pushdef([acx_cache_var], [acx_cv_fc_pp_srcext_$1])dnl
   m4_pushdef([ac_cv_fc_pp_srcext_$1], [acx_cache_var])dnl
   m4_pushdef([AC_FC_PP_SRCEXT],
     m4_bpatsubsts(m4_dquote(m4_defn([AC_FC_PP_SRCEXT])),
       ["-e Z"], ["-e T" "-e Z"],
       [AC_LANG_P\(OP\|USH\)(\[?Fortran\]?)]))dnl
   AC_FC_PP_SRCEXT([$1], [], [:])
   m4_popdef([AC_FC_PP_SRCEXT])dnl
   m4_popdef([ac_cv_fc_pp_srcext_$1])dnl
   AS_VAR_IF([acx_cache_var], [unknown],
     [ac_ext=$acx_ext_save
      m4_default([$3],
        [AC_MSG_FAILURE([unable to detect Fortran compiler flag needed to dnl
compile preprocessed .$1 files])])],
     [ac_ext=$ac_fc_srcext
      $2])
   m4_popdef([acx_cache_var])])
