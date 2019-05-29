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

# ACX_SHLIB_RPATH_FLAGS_CHECK(RPATH_FLAGS)
#                             [ACTION-IF-SUCCESS],
#                             [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Expands to a shell script that checks whether the current compiler accepts
# the automatically generated RPATH flags RPATH_FLAGS by trying to link a dummy
# program with LDFLAGS set to "RPATH_FLAGS $LDFLAGS".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
AC_DEFUN([ACX_SHLIB_RPATH_FLAGS_CHECK],
  [acx_shlib_rpath_flags_check_result=no
   AC_MSG_CHECKING([whether _AC_LANG compiler accepts the automatically dnl
generated RPATH flags])
   acx_save_LDFLAGS=$LDFLAGS
   LDFLAGS="$1 $LDFLAGS"
   AC_LINK_IFELSE([AC_LANG_PROGRAM],
     [acx_shlib_rpath_flags_check_result=yes])
   LDFLAGS=$acx_save_LDFLAGS
   AC_MSG_RESULT([$acx_shlib_rpath_flags_check_result])
   AS_VAR_IF([acx_shlib_rpath_flags_check_result], [yes], [$2],
     [m4_default([$3], [AC_MSG_FAILURE([_AC_LANG compiler does not accept dnl
the automatically generated RPATH flags '$1'])])])])

# ACX_SHLIB_PATH_VAR()
# -----------------------------------------------------------------------------
# Originally taken from Libtools where it is part of libtool.m4.
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
        [aix3*], [acx_cv_shlib_path_var=LIBPATH],
        [aix[[4-9]]*],
           [AS_VAR_IF([host_cpu],
              [ia64], [acx_cv_shlib_path_var=LD_LIBRARY_PATH],
              [acx_cv_shlib_path_var=LIBPATH])],
        [beos* | haiku*], [acx_cv_shlib_path_var=LIBRARY_PATH],
        [cygwin* | mingw* | pw32* | cegcc*], [acx_cv_shlib_path_var=PATH],
        [darwin* | rhapsody*], [acx_cv_shlib_path_var=DYLD_LIBRARY_PATH],
        [hpux9* | hpux10* | hpux11*],
           [AS_CASE([$host_cpu],
              [ia64* | hppa*64*], [acx_cv_shlib_path_var=LD_LIBRARY_PATH],
              [acx_cv_shlib_path_var=SHLIB_PATH])],
        [irix6*],
           [AS_CASE([$LD],
              [*-n32|*"-n32 "|*-melf32bmipn32|*"-melf32bmipn32 "],
                 [acx_cv_shlib_path_var=LD_LIBRARYN32_PATH],
              [*-64|*"-64 "|*-melf64bmip|*"-melf64bmip "],
                 [acx_cv_shlib_path_var=LD_LIBRARY64_PATH],
                 [acx_cv_shlib_path_var=LD_LIBRARY_PATH])],
        [os2*], [acx_cv_shlib_path_var=BEGINLIBPATH],
        [acx_cv_shlib_path_var=LD_LIBRARY_PATH])])])

# _ACX_SHLIB_RPATH_FLAG()
# -----------------------------------------------------------------------------
# Sets the result to the compiler    flag needed to add a directory to the runtime
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
