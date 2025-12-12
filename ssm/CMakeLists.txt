add_library(
  ssm # cmake-format: sort
          mo_optics_ssm.F90
          mo_optics_ssm_kernels.F90
)

if(BUILD_SHARED_LIBS)
  set_target_properties(ssm PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(
  ssm PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY}
)
