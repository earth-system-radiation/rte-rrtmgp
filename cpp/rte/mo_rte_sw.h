
#pragma once

#include "const.h"
#include "mo_optical_props.h"
#include "mo_source_functions.h"
#include "mo_fluxes.h"
#include "expand_and_transpose.h"
#include "mo_rte_solver_kernels.h"

// This code is part of Radiative Transfer for Energetics (RTE)
//
// Contacts: Robert Pincus and Eli Mlawer
// email:  rrtmgp@aer.com
//
// Copyright 2015-2018,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
// -------------------------------------------------------------------------------------------------
//
//  Contains a single routine to compute direct and diffuse fluxes of solar radiation given
//    atmospheric optical properties on a spectral grid
//    information about vertical ordering
//    boundary conditions
//      solar zenith angle, spectrally-resolved incident colimated flux, surface albedos for direct and diffuse radiation
//    optionally, a boundary condition for incident diffuse radiation
//
// It is the user's responsibility to ensure that boundary conditions (incident fluxes, surface albedos) are on the same
//   spectral grid as the optical properties.
//
// Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
//   whatever summary the user needs.
//
// The routine does error checking and choses which lower-level kernel to invoke based on
//   what kinds of optical properties are supplied
//
// -------------------------------------------------------------------------------------------------


void rte_sw(OpticalProps2str const &atmos, bool top_at_1, real1d const &mu0, real2d const &inc_flux,
            real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, FluxesBroadband &fluxes, real2d const &inc_flux_dif=real2d());



