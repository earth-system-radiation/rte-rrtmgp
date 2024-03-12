#pragma once
#include "rrtmgp_const.h"

#ifdef RRTMGP_ENABLE_YAKL
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &bnd_lims, real3d const &spectral_flux, real3d &byband_flux);
void net_byband(int ncol, int nlev, int nbnd, real3d const &bnd_flux_dn, real3d const &bnd_flux_up, real3d &bnd_flux_net);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2dk const &bnd_lims, real3dk const &spectral_flux, real3dk &byband_flux);
void net_byband(int ncol, int nlev, int nbnd, real3dk const &bnd_flux_dn, real3dk const &bnd_flux_up, real3dk &bnd_flux_net);
#endif
