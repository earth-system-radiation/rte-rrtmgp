add_library(rrtmgpkernels OBJECT)

if(KERNEL_MODE STREQUAL "extern")
  target_sources(
    rrtmgpkernels
    PRIVATE # cmake-format: sort
            api/mo_cloud_optics_rrtmgp_kernels.F90
            api/mo_gas_optics_rrtmgp_kernels.F90
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

if(BUILD_C_HEADERS)
  get_target_property(RRTMGPKERNELS_SOURCES rrtmgpkernels SOURCES)

  if(RRTMGPKERNELS_SOURCES)
    set(RRTMGP_HEADER ${CMAKE_BINARY_DIR}/headers/rrtmgp_kernels.h)
    # Convert sources to absolute paths
    set(ABSOLUTE_RRTMGPKERNELS_SOURCES "")

    foreach(src ${RRTMGPKERNELS_SOURCES})
      get_filename_component(
        ABS_PATH ${src} ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR}
      )
      list(APPEND ABSOLUTE_RRTMGPKERNELS_SOURCES ${ABS_PATH})
    endforeach()

    add_custom_command(
      OUTPUT ${RRTMGP_HEADER}
      COMMAND
        ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cbind_generator.py -if
        ${ABSOLUTE_RRTMGPKERNELS_SOURCES} -of ${RRTMGP_HEADER}
      COMMENT "Generating RRTMGPKernels' C bindings header"
      DEPENDS ${ABSOLUTE_RRTMGPKERNELS_SOURCES}
    )

    add_custom_target(
      rrtmgpkernels_header
      DEPENDS ${RRTMGP_HEADER}
      COMMENT
        "Ensure rrtmgpkernels_header is built before rrtmgpkernels can use its output"
    )

    if(KERNEL_MODE STREQUAL "extern")
      # Ensure that rrtmgpkernels depends on rrtmgpkernels_header
      add_dependencies(rrtmgpkernels rrtmgpkernels_header)

      target_sources(
        rrtmgpkernels
        PRIVATE # cmake-format: sort
                ${RRTMGP_HEADER}
      )

      install(FILES ${RRTMGP_HEADER} TYPE INCLUDE)
    endif()
  else()
    message(WARNING "No sources found for rrtmgpkernels!")
  endif()
endif()
