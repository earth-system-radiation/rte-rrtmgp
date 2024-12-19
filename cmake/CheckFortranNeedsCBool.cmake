# ~~~
# check_fortran_needs_c_bool(<variable>)
# ~~~
# Checks whether the Fortran compiler requires the c_bool kind of the logical
# type. Sets the cache <variable> to the result of the check.
#
function(check_fortran_needs_cbool var)
  if(DEFINED ${var})
    return()
  endif()

  if(NOT CMAKE_REQUIRED_QUIET)
    message(CHECK_START "Checking if the Fortran compiler requires C_BOOL")
  endif()

  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

  set(check_source_code
      "
      subroutine conftest_foo(a) bind(C)
      use iso_c_binding
      implicit none
#ifdef RTE_USE_CBOOL
      integer, parameter :: wl = c_bool
#else
      integer, parameter :: wl = kind(.true.)
#endif
      logical(wl) :: a
      end subroutine
"
  )

  set(check_source_file
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90"
  )

  unset(check_result)
  file(WRITE "${check_source_file}" "${check_source_code}")
  try_compile(try_result "${CMAKE_BINARY_DIR}" "${check_source_file}")
  if(try_result)
    set(check_result FALSE)
  else()
    file(WRITE "${check_source_file}" "${check_source_code}")
    try_compile(
      try_result "${CMAKE_BINARY_DIR}"
      "${check_source_file}"
      COMPILE_DEFINITIONS "-DRTE_USE_CBOOL"
    )
    if(try_result)
      set(check_result TRUE)
    endif()
  endif()

  if(NOT CMAKE_REQUIRED_QUIET)
    if(NOT DEFINED check_result)
      message(CHECK_FAIL "unknown (assuming no)")
    elseif(check_result)
      message(CHECK_PASS "yes")
    else()
      message(CHECK_PASS "no")
    endif()
  endif()

  if(NOT DEFINED check_result)
    set(check_result FALSE)
  endif()

  set(${var}
      "${check_result}"
      CACHE BOOL "Whether the Fortran compiler requires CBOOL type"
  )
endfunction()
