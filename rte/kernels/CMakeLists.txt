add_library(rtekernels OBJECT mo_rte_kind.F90)

if(KERNEL_MODE STREQUAL "extern")
  target_sources(
    rtekernels
    PRIVATE # cmake-format: sort
            api/mo_fluxes_broadband_kernels.F90
            api/mo_optical_props_kernels.F90
            api/mo_rte_solver_kernels.F90
            api/mo_rte_util_array.F90
  )
else()
  target_sources(
    rtekernels
    PRIVATE # cmake-format: sort
            mo_fluxes_broadband_kernels.F90
            mo_rte_util_array.F90
  )
  if(KERNEL_MODE STREQUAL "accel")
    target_sources(
      rtekernels
      PRIVATE # cmake-format: sort
              accel/mo_optical_props_kernels.F90
              accel/mo_rte_solver_kernels.F90
    )
  else()
    target_sources(
      rtekernels
      PRIVATE # cmake-format: sort
              mo_optical_props_kernels.F90
              mo_rte_solver_kernels.F90
    )
  endif()
endif()

if(BUILD_SHARED_LIBS)
  set_target_properties(rtekernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_compile_definitions(
  rtekernels PRIVATE $<$<BOOL:${RTE_ENABLE_SP}>:RTE_USE_SP>
)

target_include_directories(rtekernels PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

if(BUILD_C_HEADERS)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/api/rte_types.h.in
    ${CMAKE_BINARY_DIR}/headers/rte_types.h
    @ONLY
  )

  get_target_property(RTEKERNELS_SOURCES rtekernels SOURCES)

  if(RTEKERNELS_SOURCES)
    set(RTE_HEADER ${CMAKE_BINARY_DIR}/headers/rte_kernels.h)
    # Convert sources to absolute paths
    set(ABSOLUTE_RTEKERNELS_SOURCES "")

    foreach(src ${RTEKERNELS_SOURCES})
      get_filename_component(
        ABS_PATH ${src} ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR}
      )
      list(APPEND ABSOLUTE_RTEKERNELS_SOURCES ${ABS_PATH})
    endforeach()

    add_custom_command(
      OUTPUT ${RTE_HEADER}
      COMMAND
        ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cbind_generator.py -if
        ${ABSOLUTE_RTEKERNELS_SOURCES} -of ${RTE_HEADER}
      COMMENT "Generating RTEKernels' C bindings header"
      DEPENDS ${ABSOLUTE_RTEKERNELS_SOURCES}
    )

    add_custom_target(
      rtekernels_header
      DEPENDS ${RTE_HEADER}
      COMMENT
        "Ensure rtekernels_header is built before rtekernels can use its output"
    )
    # Ensure that rtekernels depends on rtekernels_header
    add_dependencies(rtekernels rtekernels_header)

    if(KERNEL_MODE STREQUAL "extern")
      target_sources(
        rtekernels
        PRIVATE # cmake-format: sort
                ${CMAKE_BINARY_DIR}/headers/rte_types.h
                ${RTE_HEADER}
      )
    endif()

    install(FILES ${RTE_HEADER} TYPE INCLUDE)
  else()
    message(WARNING "No sources found for rtekernels!")
  endif()

  install(FILES ${CMAKE_BINARY_DIR}/headers/rte_types.h TYPE INCLUDE)
endif()
