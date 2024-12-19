find_library(
  NetCDF_Fortran_LIBRARY
  NAMES netcdff
  DOC "NetCDF-Fortran library"
)
mark_as_advanced(NetCDF_Fortran_LIBRARY)

find_path(
  NetCDF_Fortran_INCLUDE_DIR
  NAMES netcdf.mod NETCDF.mod
  DOC "NetCDF_Fortran include directory"
)
mark_as_advanced(NetCDF_Fortran_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  NetCDF_Fortran
  REQUIRED_VARS NetCDF_Fortran_LIBRARY NetCDF_Fortran_INCLUDE_DIR
)
