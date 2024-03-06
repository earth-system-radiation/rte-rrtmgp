
#pragma once

#include "rrtmgp_const.h"
#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include <string>

#ifdef RRTMGP_ENABLE_YAKL
void load_and_init(GasOpticsRRTMGP &kdist, std::string filename, GasConcs const &available_gases);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
void load_and_init(GasOpticsRRTMGPK &kdist, std::string filename, GasConcsK const &available_gases);
#endif
