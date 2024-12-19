# ~~~
# Build a list of files to compile based on input file names list
# For each file in the FILENAMES input array, we check if
# the path ${CMAKE_CURRENT_SOURCE_DIR}/${KERNEL_MODE}/{$FILENAME}
# is a valid path and if it is, we append the file to the output list
# If not, we append the default file instead.
# ~~~
function(find_files_in_folders OUTPUT_VAR FILENAMES KERNEL_MODE)
  # Variable to hold the output list
  if(KERNEL_MODE STREQUAL "accel")
    set(KERNEL_SUB_FOLDER "accel")
  elseif(KERNEL_MODE STREQUAL "extern")
    set(KERNEL_SUB_FOLDER "api")
  elseif(KERNEL_MODE STREQUAL "default")
    set(KERNEL_SUB_FOLDER ".")
  else()
    message(FATAL_ERROR "Invalid KERNEL_MODE ${KERNEL_MODE}")
  endif()

  set(RESULT_LIST) # Initialize an empty list to store results

  foreach(FILENAME IN LISTS FILENAMES)
    set(FULL_PATH
        "${CMAKE_CURRENT_SOURCE_DIR}/${KERNEL_SUB_FOLDER}/${FILENAME}"
    )
    if(EXISTS "${FULL_PATH}")
      list(APPEND RESULT_LIST "${FULL_PATH}")
      message(STATUS "Using ${FULL_PATH}")
      continue()
    endif()
    set(FULL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
    list(APPEND RESULT_LIST "${FULL_PATH}")
    message(STATUS "Using ${FULL_PATH}")
  endforeach()

  # Return the result list to the caller
  set(${OUTPUT_VAR}
      "${RESULT_LIST}"
      PARENT_SCOPE
  )
endfunction()
