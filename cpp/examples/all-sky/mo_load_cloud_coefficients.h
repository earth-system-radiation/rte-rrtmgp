
#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props.h"
#include "mo_cloud_optics.h"

#ifdef RRTMGP_ENABLE_YAKL
void load_cld_lutcoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);

void load_cld_padecoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
void load_cld_lutcoeff(CloudOpticsK &cloud_spec, std::string cld_coeff_file);

void load_cld_padecoeff(CloudOpticsK &cloud_spec, std::string cld_coeff_file);
#endif
