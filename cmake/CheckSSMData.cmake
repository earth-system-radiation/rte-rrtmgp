# Verify that the SSM data files are present in ${SSM_DATA}.
#
# Run as: cmake -D SSM_DATA=<dir> -P CheckSSMData.cmake
#
# The data is fetched at build time by the ssm-data ExternalProject target (part
# of the default ALL target). This script is used as a test fixture so that
# tests requiring the data fail with an actionable message if the project was
# not built before running ctest.

set(required_files
    # cmake-format: sort
    ssm-lw-ckdmip.nc
    ssm-lw-rce.nc
    ssm-lw-rfmip.nc
    ssm-sw-ckdmip.nc
    ssm-sw-rce.nc
    ssm-sw-rfmip.nc
)

set(missing_files "")
foreach(data_file IN LISTS required_files)
  if(NOT EXISTS "${SSM_DATA}/${data_file}")
    list(APPEND missing_files "${SSM_DATA}/${data_file}")
  endif()
endforeach()

if(missing_files)
  string(REPLACE ";" "\n  " missing_pretty "${missing_files}")
  message(
    FATAL_ERROR
      "SSM data files are missing:\n  ${missing_pretty}\n"
      "Build the project first (e.g. 'cmake --build <build-dir>') so the "
      "ssm-data target fetches the data before running ctest."
  )
endif()
