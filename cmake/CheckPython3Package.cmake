# ~~~
# check_python3_package(<package>
#                       [CODE <code>])
# ~~~
# Checks whether the Python3 <package> is available by running the <code> with
# ${Python3_EXECUTABLE} (defaults to 'import <package>'). Fails the
# configuration if the result is negative.
#
function(check_python3_package package)
  cmake_parse_arguments(PARSE_ARGV 1 arg "" "CODE" "")

  if(DEFINED ${var})
    return()
  endif()

  if(NOT CMAKE_REQUIRED_QUIET)
    message(
      CHECK_START "Checking if the Python3 package ${package} is available"
    )
  endif()

  if(NOT arg_CODE)
    set(arg_CODE "import ${package}")
  endif()

  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "${arg_CODE}"
    RESULT_VARIABLE exit_status
    OUTPUT_QUIET ERROR_QUIET
  )

  if(NOT CMAKE_REQUIRED_QUIET)
    if(NOT exit_status)
      message(CHECK_PASS "yes")
    else()
      message(CHECK_FAIL "no")
    endif()
  endif()

  if(exit_status)
    message(FATAL_ERROR "Required Python3 package ${package} is not available")
  endif()
endfunction()
