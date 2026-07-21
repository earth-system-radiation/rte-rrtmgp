# Verify that the RRTMGP data files are present in ${RRTMGP_DATA}.
#
# Run as: cmake -D RRTMGP_DATA=<dir> -P CheckRRTMGPData.cmake
#
# The data is fetched at build time by the rrtmgp-data ExternalProject target
# (part of the default ALL target). This script is used as a test fixture so
# that tests requiring the data fail with an actionable message if the project
# was not built before running ctest.

set(_required_files
    # cmake-format: sort
    rrtmgp-aerosols-merra-lw.nc
    rrtmgp-aerosols-merra-sw.nc
    rrtmgp-clouds-lw-bnd.nc
    rrtmgp-clouds-sw-bnd.nc
    rrtmgp-gas-lw-g128.nc
    rrtmgp-gas-lw-g256.nc
    rrtmgp-gas-sw-g112.nc
    rrtmgp-gas-sw-g224.nc
)

set(_missing "")
foreach(_file IN LISTS _required_files)
  if(NOT EXISTS "${RRTMGP_DATA}/${_file}")
    list(APPEND _missing "${RRTMGP_DATA}/${_file}")
  endif()
endforeach()

if(_missing)
  string(REPLACE ";" "\n  " _missing_pretty "${_missing}")
  message(
    FATAL_ERROR
      "RRTMGP data files are missing:\n  ${_missing_pretty}\n"
      "Build the project first (e.g. 'cmake --build <build-dir>') so the "
      "rrtmgp-data target fetches the data before running ctest."
  )
endif()
