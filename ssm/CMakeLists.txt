add_library(
  ssm # cmake-format: sort
  mo_optics_ssm.F90 mo_optics_ssm_kernels.F90
)

add_library(rte-rrtmgp::ssm ALIAS ssm)

if(BUILD_SHARED_LIBS)
  set_target_properties(ssm PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(
  ssm
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)

target_link_libraries(ssm PUBLIC rte)

set_target_properties(
  ssm
  PROPERTIES EXPORT_NAME rte-rrtmgp::ssm
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

install(TARGETS ssm EXPORT rte-rrtmgp-targets)
