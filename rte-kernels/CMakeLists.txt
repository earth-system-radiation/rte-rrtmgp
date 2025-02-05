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
  rtekernels
  PRIVATE $<$<BOOL:${RTE_ENABLE_SP}>:RTE_USE_SP> rtekernels
  PRIVATE $<$<BOOL:${RTE_BOOL}>:RTE_BOOL> rtekernels
  PRIVATE $<$<BOOL:${RTE_PRECISION}>:RTE_PRECISION> rtekernels
  PRIVATE $<$<BOOL:${RTE_EPSILON}>:RTE_EPSILON>
)

target_include_directories(rtekernels PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

if(BUILD_C_HEADERS)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/api/rte_types.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/api/rte_types.h
    @ONLY
  )

  get_target_property(RTEKERNELS_SOURCES rtekernels SOURCES)

  if(RTEKERNELS_SOURCES)
    # Convert sources to absolute paths
    set(ABSOLUTE_RTEKERNELS_SOURCES "")

    foreach(SRC ${RTEKERNELS_SOURCES})
      get_filename_component(
        ABS_PATH ${SRC} ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR}
      )
      list(APPEND ABSOLUTE_RTEKERNELS_SOURCES ${ABS_PATH})
    endforeach()

    add_custom_target(
      rtekernels_header ALL
      COMMAND
        ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cbind_generator.py -if
        ${ABSOLUTE_RTEKERNELS_SOURCES} -of
        ${CMAKE_CURRENT_SOURCE_DIR}/api/rte_kernels.h
    )

    if(KERNEL_MODE STREQUAL "extern")
      target_sources(
        rtekernels
        PRIVATE # cmake-format: sort
                api/rte_kernels.h
                api/rte_types.h
      )
    endif()
  else()
    message(WARNING "No sources found for rtekernels!")
  endif()
endif()
