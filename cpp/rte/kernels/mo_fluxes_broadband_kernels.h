
#pragma once

#include "rrtmgp_const.h"


// Spectral reduction over all points
void sum_broadband(int ncol, int nlev, int ngpt, real3d const &spectral_flux, real2d const &broadband_flux);


// Net flux: Spectral reduction over all points
void net_broadband(int ncol, int nlev, int ngpt, real3d const &spectral_flux_dn, real3d const &spectral_flux_up, real2d const &broadband_flux_net);


// Net flux when bradband flux up and down are already available
void net_broadband(int ncol, int nlev, real2d const &flux_dn, real2d const &flux_up, real2d const &broadband_flux_net);


