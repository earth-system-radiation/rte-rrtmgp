cmake_minimum_required(VERSION 3.0)
project(standalone)

set(CPP_SRC 
    ../mo_load_coefficients.cpp
    ../../extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
    mo_garand_atmos_io.cpp 
    mo_load_cloud_coefficients.cpp
    rrtmgp_allsky.cpp
)
add_executable(allsky ${CPP_SRC})

target_link_libraries(allsky yakl rrtmgp ${NCFLAGS})
target_compile_features(allsky PUBLIC cxx_std_14)
target_include_directories(allsky PUBLIC ../)
target_include_directories(allsky PUBLIC ../../)
target_include_directories(allsky PUBLIC ../../extensions/cloud_optics)
target_include_directories(allsky PUBLIC ../../extensions/fluxes_byband)
target_include_directories(allsky PUBLIC ../../rte)
target_include_directories(allsky PUBLIC ../../rte/kernels)
target_include_directories(allsky PUBLIC ../../rrtmgp)
target_include_directories(allsky PUBLIC ../../rrtmgp/kernels)
target_include_directories(allsky PUBLIC ${NCINCLUDE})
process_cxx_source_files(${CPP_SRC})
