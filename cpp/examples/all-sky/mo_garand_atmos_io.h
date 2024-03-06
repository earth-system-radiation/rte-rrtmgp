
#pragma once
#include "rrtmgp_const.h"
#include "YAKL.h"
#include "YAKL_netcdf.h"
#include "mo_gas_concentrations.h"

#ifdef RRTMGP_ENABLE_YAKL
void read_atmos(std::string input_file, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, real2d &col_dry, int ncol);


void write_sw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, real2d const &flux_dir, int ncol);


void write_lw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, int ncol);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
void read_atmos(std::string input_file, real2dk &p_lay, real2dk &t_lay, real2dk &p_lev, real2dk &t_lev,
                GasConcsK &gas_concs, real2dk &col_dry, int ncol);


void write_sw_fluxes(std::string fileName, real2dk const &flux_up, real2dk const &flux_dn, real2dk const &flux_dir, int ncol);


void write_lw_fluxes(std::string fileName, real2dk const &flux_up, real2dk const &flux_dn, int ncol);
#endif
