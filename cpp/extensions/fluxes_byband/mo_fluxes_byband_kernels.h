#pragma once
#include "const.h"
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &bnd_lims, real3d const &spectral_flux, real3d &byband_flux);
void net_byband(int ncol, int nlev, int nbnd, real3d const &bnd_flux_dn, real3d const &bnd_flux_up, real3d &bnd_flux_net);
