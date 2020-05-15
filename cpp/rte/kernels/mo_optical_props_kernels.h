
#pragma once

#include "const.h"



// increment 2stream by 2stream
void inc_2stream_by_2stream_bybnd(int ncol, int nlay, int ngpt,
                                  real3d       &tau1, real3d       &ssa1, real3d       &g1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2,
                                  int nbnd, int2d const &gpt_lims);



// Incrementing when the second set of optical properties is defined at lower spectral resolution
//   (e.g. by band instead of by gpoint)
void inc_1scalar_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2,
                                  int nbnd, int2d const &gpt_lims);



// Delta-scale
//   f = g*g
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, real3d &tau, real3d &ssa, real3d &g);


// Delta-scaling, provided only for two-stream properties at present
// -------------------------------------------------------------------------------------------------
// Delta-scale
//   user-provided value of f (forward scattering)
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, real3d &tau, real3d &ssa, real3d &g, real3d const &f);



// Addition of optical properties: the first set are incremented by the second set.
//
//   There are three possible representations of optical properties (scalar = optical depth only;
//   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
//   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
//   each choice of representation on the left hand side times three representations of the
//   optical properties to be added.
//
//   There are two sets of these nine routines. In the first the two sets of optical
//   properties are defined at the same spectral resolution. There is also a set of routines
//   to add properties defined at lower spectral resolution to a set defined at higher spectral
//   resolution (adding properties defined by band to those defined by g-point)
void increment_1scalar_by_1scalar(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2);



// increment 1scalar by 2stream
void increment_1scalar_by_2stream(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2, real3d const &ssa2);



// increment 1scalar by nstream
void increment_1scalar_by_nstream(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2, real3d const &ssa2);



// increment 2stream by 1scalar
void increment_2stream_by_1scalar(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d const &tau2);



// increment 2stream by 2stream
void increment_2stream_by_2stream(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d &g1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2);


// increment 2stream by nstream
void increment_2stream_by_nstream(int ncol, int nlay, int ngpt, int nmom2, real3d &tau1, real3d &ssa1,
                                  real3d &g1, real3d const &tau2, real3d const &ssa2, real4d const &p2);



// increment nstream by 1scalar
void increment_nstream_by_1scalar(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d const &tau2);



// increment nstream by 2stream
void increment_nstream_by_2stream(int ncol, int nlay, int ngpt, int nmom1, real3d &tau1, real3d &ssa1,
                                  real4d &p1, real3d const &tau2, real3d const &ssa2, real3d const &g2);



// increment nstream by nstream
void increment_nstream_by_nstream(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real3d &tau1,
                                  real3d &ssa1, real4d &p1, real3d const &tau2, real3d const &ssa2, real4d const &p2);


// increment 1scalar by 2stream
void inc_1scalar_by_2stream_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2,
                                  real3d const &ssa2, int nbnd, int2d const &gpt_lims);



// increment 1scalar by nstream
void inc_1scalar_by_nstream_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2,
                                  real3d const &ssa2, int nbnd, int2d const &gpt_lims);



// increment 2stream by 1scalar
void inc_2stream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1,
                                  real3d const &tau2, int nbnd, int2d const &gpt_lims);


// increment 2stream by nstream
void inc_2stream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom2, real3d &tau1, real3d &ssa1,
                                  real3d &g1, real3d const &tau2, real3d const &ssa2, real4d const &p2,
                                  int nbnd, int2d const &gpt_lims);



// increment nstream by 1scalar
void inc_nstream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d const &tau2,
                                  int nbnd, int2d const &gpt_lims);



// increment nstream by 2stream
void inc_nstream_by_2stream_bybnd(int ncol, int nlay, int ngpt, int nmom1, real3d &tau1, real3d &ssa1, real4d &p1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2, int nbnd,
                                  int2d const &gpt_lims);



// increment nstream by nstream
void inc_nstream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real3d &tau1, real3d &ssa1, real4d &p1,
                                  real3d const &tau2, real3d const &ssa2, real4d const &p2, int nbnd, int2d const &gpt_lims);



// Subsetting, meaning extracting some portion of the 3D domain
void extract_subset_dim1_3d(int ncol, int nlay, int ngpt, real3d const &array_in, int colS, int colE, real3d &array_out);



void extract_subset_dim2_4d(int nmom, int ncol, int nlay, int ngpt, real3d const &array_in, int colS, int colE, real3d &array_out);



// Extract the absorption optical thickness which requires mulitplying by 1 - ssa
void extract_subset_absorption_tau(int ncol, int nlay, int ngpt, real3d const &tau_in, real3d const &ssa_in, int colS, int colE,
                                   real3d &tau_out);




