# Verify that the rte-examples data files are present in ${RTE_EXAMPLES_DATA}.
#
# Run as: cmake -D RTE_EXAMPLES_DATA=<dir> -P CheckRTEExamplesData.cmake
#
# The data is fetched at build time by the rte-examples-data ExternalProject
# target (part of the default ALL target). This script is used as a test fixture
# so that tests requiring the data fail with an actionable message if the
# project was not built before running ctest.

set(required_files
    # cmake-format: sort
    ckdmip-states.nc
    rce-states.nc
    rfmip-states.nc
)

set(missing_files "")
foreach(data_file IN LISTS required_files)
  if(NOT EXISTS "${RTE_EXAMPLES_DATA}/${data_file}")
    list(APPEND missing_files "${RTE_EXAMPLES_DATA}/${data_file}")
  endif()
endforeach()

if(missing_files)
  string(REPLACE ";" "\n  " missing_pretty "${missing_files}")
  message(
    FATAL_ERROR
      "rte-examples data files are missing:\n  ${missing_pretty}\n"
      "Build the project first (e.g. 'cmake --build <build-dir>') so the "
      "rte-examples-data target fetches the data before running ctest."
  )
endif()
