add_library(rrtmgpkernels OBJECT)

if(KERNEL_MODE STREQUAL "extern")
  target_sources(
    rrtmgpkernels
    PRIVATE # cmake-format: sort
            api/mo_cloud_optics_rrtmgp_kernels.F90
            api/mo_gas_optics_rrtmgp_kernels.F90
            api/rrtmgp_kernels.h
  )
else()
  target_sources(
    rrtmgpkernels
    PRIVATE # cmake-format: sort
            mo_cloud_optics_rrtmgp_kernels.F90
  )
  if(KERNEL_MODE STREQUAL "accel")
    target_sources(
      rrtmgpkernels
      PRIVATE # cmake-format: sort
              accel/mo_gas_optics_rrtmgp_kernels.F90
    )
  else()
    target_sources(
      rrtmgpkernels
      PRIVATE # cmake-format: sort
              mo_gas_optics_rrtmgp_kernels.F90
    )
  endif()
endif()

if(BUILD_SHARED_LIBS)
  set_target_properties(rrtmgpkernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(
  rrtmgpkernels PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(rrtmgpkernels PRIVATE rte)
