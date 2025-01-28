# ~~~
# check_compiler_is_nec(<variable> <language>)
# ~~~
# Checks whether the <language> compiler in use is actually the NEC one (CMake
# does not support NEC compiler yet and usually mistakes it for GNU). Sets the
# internal cache <variable> to the result of the check. The requested <language>
# must be enabled. Supported languages are C, CXX and Fortran.
#
function(check_compiler_is_nec var lang)
  if(DEFINED ${var})
    return()
  endif()

  get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  if(NOT ${lang} IN_LIST languages)
    message(FATAL_ERROR "Language ${lang} is not enabled")
  endif()

  if(${lang} STREQUAL "C")
    include(CheckSymbolExists)
    set(CMAKE_REQUIRED_QUIET 1)
    check_symbol_exists(__NEC__ "" ${var})
  elseif(${lang} STREQUAL "CXX")
    include(CheckCXXSymbolExists)
    set(CMAKE_REQUIRED_QUIET 1)
    check_cxx_symbol_exists(__NEC__ "" ${var})
  elseif(${lang} STREQUAL "Fortran")
    set(check_source_file
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90"
    )
    file(
      WRITE "${check_source_file}"
      "      program main
      implicit none
#ifdef __NEC__
      integer a
#else
      choke me
#endif
#ifndef __NEC__
      choke me
#else
      integer b
#endif
      a = 4
      b = 2
      end
"
    )
    try_compile(${var} "${PROJECT_BINARY_DIR}" "${check_source_file}")
  else()
    message(FATAL_ERROR "Language ${lang} is not supported")
  endif()

  if(${var})
    set(${var} TRUE)
  else()
    set(${var} FALSE)
  endif()

  set(${var}
      "${${var}}"
      CACHE INTERNAL "Whether the ${lang} compiler is NEC"
  )
endfunction()
