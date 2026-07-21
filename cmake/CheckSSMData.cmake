# Verify that the SSM data files are present in ${SSM_DATA}.
#
# Run as: cmake -D SSM_DATA=<dir> -P CheckSSMData.cmake
#
# The data is fetched at build time by the ssm-data ExternalProject target
# (part of the default ALL target). This script is used as a test fixture so
# that tests requiring the data fail with an actionable message if the project
# was not built before running ctest.

set(_required_files # cmake-format: sort
    ssm-lw-ckdmip.nc
    ssm-lw-rce.nc
    ssm-lw-rfmip.nc
    ssm-sw-ckdmip.nc
    ssm-sw-rce.nc
    ssm-sw-rfmip.nc
)

set(_missing "")
foreach(_file IN LISTS _required_files)
  if(NOT EXISTS "${SSM_DATA}/${_file}")
    list(APPEND _missing "${SSM_DATA}/${_file}")
  endif()
endforeach()

if(_missing)
  string(REPLACE ";" "\n  " _missing_pretty "${_missing}")
  message(
    FATAL_ERROR
    "SSM data files are missing:\n  ${_missing_pretty}\n"
    "Build the project first (e.g. 'cmake --build <build-dir>') so the "
    "ssm-data target fetches the data before running ctest."
  )
endif()
