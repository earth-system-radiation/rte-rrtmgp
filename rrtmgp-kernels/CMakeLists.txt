add_library(rrtmgpkernels STATIC)

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

target_link_libraries(rrtmgpkernels PRIVATE rtekernels)
