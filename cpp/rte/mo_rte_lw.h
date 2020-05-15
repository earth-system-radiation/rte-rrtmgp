
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
//    atmospheric optical properties, spectrally-resolved
//    information about vertical ordering
//    internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
//    boundary conditions: surface emissivity defined per band
//    optionally, a boundary condition for incident diffuse radiation
//    optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected
//
// If optical properties are supplied via class ty_optical_props_1scl (absorption optical thickenss only)
//    then an emission/absorption solver is called
//    If optical properties are supplied via class ty_optical_props_2str fluxes are computed via
//    two-stream calculations and adding.
//
// It is the user's responsibility to ensure that emissivity is on the same
//   spectral grid as the optical properties.
//
// Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
//   whatever summary the user needs.
//
// The routine does error checking and choses which lower-level kernel to invoke based on
//   what kinds of optical properties are supplied
//
// -------------------------------------------------------------------------------------------------


// Interface using only optical properties and source functions as inputs; fluxes as outputs.
void rte_lw(int max_gauss_pts, real2d const &gauss_Ds, real2d const &gauss_wts, OpticalProps1scl const &optical_props, bool top_at_1, SourceFuncLW const &sources, real2d const &sfc_emis,
            FluxesBroadband &fluxes, real2d const &inc_flux=real2d(), int n_gauss_angles=-1);


