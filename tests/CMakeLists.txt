foreach(TEST rte_optic_prop_unit_tests rte_lw_solver_unit_tests rte_sw_solver_unit_tests)
 	add_executable(${TEST})
	target_sources(${TEST} PRIVATE mo_testing_utils.F90 ${TEST}.F90)
	target_link_libraries(${TEST} PRIVATE rtekernels rtef)
endforeach() 


#[[

For the rest of the tests we need to be able to look for 
	mo_heating_rates mo_compute_bc mo_rrtmgp_clr_all_sky
	mo_simple_netcdf mo_rfmip_io
	mo_testing_io mo_testing_utils
	mo_cloud_optics_rrtmgp mo_load_cloud_coefficients mo_garand_atmos_io
	mo_solar_variability
in 
	examples examples/rfmip-clear-sk examples/all-sky
    extensions extensions/solar_variability

foreach(TEST check_equivalence check_variants)
endforeach() 
]]