set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

set(extensions_source_dir ${PROJECT_SOURCE_DIR}/extensions)

add_library(
  test_utils STATIC # cmake-format: sort
  ${extensions_source_dir}/mo_compute_bc.F90
  ${extensions_source_dir}/mo_heating_rates.F90
  ${extensions_source_dir}/mo_rrtmgp_clr_all_sky.F90
  ${extensions_source_dir}/mo_zenith_angle_spherical_correction.F90
  ${extensions_source_dir}/solar_variability/mo_solar_variability.F90
  mo_comparisons.F90
  mo_rcemip_profiles.F90
)

target_include_directories(
  test_utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)
target_link_libraries(test_utils PUBLIC rte-rrtmgp::rte rte-rrtmgp::rrtmgp)

foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        check_equivalence
        check_variants
        rte_lw_solver_unit_tests
        rte_optic_prop_unit_tests
        rte_sw_solver_unit_tests
        test_zenith_angle_spherical_correction
)
  add_executable(${test_executable} ${test_executable}.F90)
  target_link_libraries(
    ${test_executable} PUBLIC test_utils rte-rrtmgp::rte rte-rrtmgp::rrtmgp
  )
endforeach()

#
# Adding tests
#
foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        rte_lw_solver_unit_tests
        rte_optic_prop_unit_tests
        rte_sw_solver_unit_tests
)
  add_test(NAME ${test_executable} COMMAND ${test_executable})
endforeach()

foreach(g_value IN ITEMS 256 128)
  add_test(
    NAME check_equivalence_lw_g${g_value}
    COMMAND
      check_equivalence
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
      ${RRTMGP_DATA}/rrtmgp-gas-lw-g${g_value}.nc
  )
  set_tests_properties(
    check_equivalence_lw_g${g_value}
    PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data
  )
endforeach()

foreach(g_value IN ITEMS 224 112)
  add_test(
    NAME check_equivalence_sw_g${g_value}
    COMMAND
      check_equivalence
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
      ${RRTMGP_DATA}/rrtmgp-gas-sw-g${g_value}.nc
  )
  set_tests_properties(
    check_equivalence_sw_g${g_value}
    PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data
  )
endforeach()

add_test(
  NAME check_variants_lw
  COMMAND
    check_variants
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    lw_flux_variants.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)
set_tests_properties(
  check_variants_lw PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data
)

add_test(
  NAME check_variants_sw
  COMMAND
    check_variants
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    sw_flux_variants.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)
set_tests_properties(
  check_variants_sw PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data
)

add_custom_target(
  validation-plots
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validation-plots.py
    --state_file
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    --lw_vars_file lw_flux_variants.nc --sw_vars_file sw_flux_variants.nc
  COMMENT "Generating validation plots"
)
