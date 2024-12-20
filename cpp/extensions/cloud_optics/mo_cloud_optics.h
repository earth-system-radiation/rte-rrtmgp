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
// Provides cloud optical properties as a function of effective radius for the RRTMGP bands
//   Based on Mie calculations for liquid
//     and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
//   Can use either look-up tables or Pade approximates according to which data has been loaded
//   Mike Iacono (AER) is the original author
//
// The class can be used as-is but is also intended as an example of how to extend the RTE framework
// -------------------------------------------------------------------------------------------------

#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props.h"

#ifdef RRTMGP_ENABLE_YAKL
class CloudOptics : public OpticalProps {
public:
  // Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  int icergh;  // (1 = none, 2 = medium, 3 = high)
  // Lookup table information
  // Upper and lower limits of the tables
  real radliq_lwr;
  real radliq_upr;
  real radice_lwr;
  real radice_upr;
  // How many steps in the table? (for convenience)
  int liq_nsteps;
  int ice_nsteps;
  // How big is each step in the table?
  real liq_step_size;
  real ice_step_size;
  // The tables themselves.
  real2d lut_extliq; // (nsize_liq, nbnd)
  real2d lut_ssaliq; // (nsize_liq, nbnd)
  real2d lut_asyliq; // (nsize_liq, nbnd)
  real3d lut_extice; // (nsize_ice, nbnd, nrghice)
  real3d lut_ssaice; // (nsize_ice, nbnd, nrghice)
  real3d lut_asyice; // (nsize_ice, nbnd, nrghice)
  // Pade approximant coefficients
  real3d pade_extliq; // (nbnd, nsizereg, ncoeff_ext)
  real3d pade_ssaliq; // (nbnd, nsizereg, ncoeff_ssa_g)
  real3d pade_asyliq; // (nbnd, nsizereg, ncoeff_ssa_g)
  real4d pade_extice; // (nbnd, nsizereg, ncoeff_ext, nrghice)
  real4d pade_ssaice; // (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  real4d pade_asyice; // (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  // Particle size regimes for Pade formulations
  real1d pade_sizreg_extliq; // (nbound)
  real1d pade_sizreg_ssaliq; // (nbound)
  real1d pade_sizreg_asyliq; // (nbound)
  real1d pade_sizreg_extice; // (nbound)
  real1d pade_sizreg_ssaice; // (nbound)
  real1d pade_sizreg_asyice; // (nbound)

  CloudOptics() = default;

  // Routines to load data needed for cloud optics calculations. Two routines: one to load
  //    lookup-tables and one for coefficients for Pade approximates
  void load(real2d const &band_lims_wvn, real radliq_lwr, real radliq_upr, real radliq_fac, real radice_lwr, real radice_upr,
           real radice_fac, real2d const &lut_extliq, real2d const &lut_ssaliq, real2d const &lut_asyliq,
           real3d const &lut_extice, real3d const &lut_ssaice, real3d const &lut_asyice) {
    using yakl::intrinsics::size;

    // Local variables
    this->init(band_lims_wvn, int2d(), "RRTMGP cloud optics");
    // LUT coefficient dimensions
    int nsize_liq = size(lut_extliq,1);
    int nsize_ice = size(lut_extice,1);
    int nbnd      = size(lut_extliq,2);
    int nrghice   = size(lut_extice,3);
    //
    // Error checking
    //   Can we check for consistency between table bounds and _fac?
    //
    if (nbnd != this->get_nband()) {
      stoprun("cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization");
    }
    if (size(lut_extice,2) != nbnd) { stoprun("cloud_optics%init(): array lut_extice has the wrong number of bands"); }
    if (size(lut_ssaliq,1) != nsize_liq || size(lut_ssaliq,2) != nbnd) {
      stoprun("cloud_optics%init(): array lut_ssaliq isn't consistently sized");
    }
    if (size(lut_asyliq,1) != nsize_liq || size(lut_asyliq,2) != nbnd) {
      stoprun("cloud_optics%init(): array lut_asyliq isn't consistently sized");
    }
    if (size(lut_ssaice,1) != nsize_ice || size(lut_ssaice,2) != nbnd || size(lut_ssaice,3) != nrghice) {
      stoprun("cloud_optics%init(): array lut_ssaice  isn't consistently sized");
    }
    if (size(lut_asyice,1) != nsize_ice || size(lut_asyice,2) != nbnd || size(lut_asyice,3) != nrghice) {
      stoprun("cloud_optics%init(): array lut_asyice  isn't consistently sized");
    }

    this->liq_nsteps = nsize_liq;
    this->ice_nsteps = nsize_ice;
    this->liq_step_size = (radliq_upr - radliq_lwr) / (nsize_liq-1.);
    this->ice_step_size = (radice_upr - radice_lwr) / (nsize_ice-1.);
    // Load LUT constants
    this->radliq_lwr = radliq_lwr;
    this->radliq_upr = radliq_upr;
    this->radice_lwr = radice_lwr;
    this->radice_upr = radice_upr;
    // Load LUT coefficients
    this->lut_extliq = lut_extliq;
    this->lut_ssaliq = lut_ssaliq;
    this->lut_asyliq = lut_asyliq;
    this->lut_extice = lut_extice;
    this->lut_ssaice = lut_ssaice;
    this->lut_asyice = lut_asyice;
    // Set default ice roughness - min values
    this->set_ice_roughness(1);
  }



  // Cloud optics initialization function - Pade
  void load(real2d const &band_lims_wvn, real3d const &pade_extliq, real3d const &pade_ssaliq, real3d const &pade_asyliq,
            real4d const &pade_extice, real4d const &pade_ssaice, real4d const &pade_asyice,
            real1d const &pade_sizreg_extliq, real1d const &pade_sizreg_ssaliq, real1d const &pade_sizreg_asyliq,
            real1d const &pade_sizreg_extice, real1d const &pade_sizreg_ssaice, real1d const &pade_sizreg_asyice) {
    using yakl::intrinsics::size;

    // Pade coefficient dimensions
    int nbnd         = size(pade_extliq,1);
    int nsizereg     = size(pade_extliq,2);
    int ncoeff_ext   = size(pade_extliq,3);
    int ncoeff_ssa_g = size(pade_ssaliq,3);
    int nrghice      = size(pade_extice,4);
    int nbound       = size(pade_sizreg_extliq);
    // The number of size regimes is assumed in the Pade evaluations
    if (nsizereg != 3) {
      stoprun("cloud optics: code assumes exactly three size regimes for Pade approximants but data is otherwise");
    }
    this->init(band_lims_wvn, int2d(), "RRTMGP cloud optics");
    // Error checking
    if (nbnd != this->get_nband()) {
      stoprun("cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization");
    }
    if (size(pade_ssaliq,1) != nbnd || size(pade_ssaliq,2) != nsizereg || size(pade_ssaliq,3) != ncoeff_ssa_g) {
      stoprun("cloud_optics%init(): array pade_ssaliq isn't consistently sized");
    }
    if (size(pade_asyliq,1) != nbnd || size(pade_asyliq,2) != nsizereg || size(pade_asyliq,3) != ncoeff_ssa_g) {
      stoprun("cloud_optics%init(): array pade_asyliq isn't consistently sized");
    }
    if (size(pade_extice,1) != nbnd || size(pade_extice,2) != nsizereg || size(pade_extice,3) != ncoeff_ext || size(pade_extice,4) != nrghice) {
      stoprun("cloud_optics%init(): array pade_extice isn't consistently sized");
    }
    if (size(pade_ssaice,1) != nbnd || size(pade_ssaice,2) != nsizereg || size(pade_ssaice,3) != ncoeff_ssa_g || size(pade_ssaice,4) != nrghice) {
      stoprun("cloud_optics%init(): array pade_ssaice isn't consistently sized");
    }
    if (size(pade_asyice,1) != nbnd || size(pade_asyice,2) != nsizereg || size(pade_asyice,3) != ncoeff_ssa_g || size(pade_asyice,4) != nrghice) {
      stoprun("cloud_optics%init(): array pade_asyice isn't consistently sized");
    }
    if (size(pade_sizreg_ssaliq,1) != nbound || size(pade_sizreg_asyliq,1) != nbound || size(pade_sizreg_extice,1) != nbound ||
        size(pade_sizreg_ssaice,1) != nbound || size(pade_sizreg_asyice,1) != nbound ) {
      stoprun("cloud_optics%init(): one or more Pade size regime arrays are inconsistently sized");
    }
    if (nsizereg != 3) { stoprun("cloud_optics%init(): Expecting precisely three size regimes for Pade approximants"); }
    // Consistency among size regimes
    this->radliq_lwr = pade_sizreg_extliq(1);
    this->radliq_upr = pade_sizreg_extliq(nbound);
    this->radice_lwr = pade_sizreg_extice(1);
    this->radice_upr = pade_sizreg_extice(nbound);
    if (pade_sizreg_ssaliq(1) < this->radliq_lwr || pade_sizreg_asyliq(1) < this->radliq_lwr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have inconsistent lowest values");
    }
    if (pade_sizreg_ssaice(1) < this->radice_lwr || pade_sizreg_asyice(1) < this->radice_lwr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have inconsistent lower values");
    }
    if (pade_sizreg_ssaliq(nbound) > this->radliq_upr || pade_sizreg_asyliq(nbound) > this->radliq_upr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_upr");
    }
    if (pade_sizreg_ssaice(nbound) > this->radice_upr || pade_sizreg_asyice(nbound) > this->radice_upr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_upr");
    }
    // Load data
    this->pade_extliq        = pade_extliq       ;
    this->pade_ssaliq        = pade_ssaliq       ;
    this->pade_asyliq        = pade_asyliq       ;
    this->pade_extice        = pade_extice       ;
    this->pade_ssaice        = pade_ssaice       ;
    this->pade_asyice        = pade_asyice       ;
    this->pade_sizreg_extliq = pade_sizreg_extliq;
    this->pade_sizreg_ssaliq = pade_sizreg_ssaliq;
    this->pade_sizreg_asyliq = pade_sizreg_asyliq;
    this->pade_sizreg_extice = pade_sizreg_extice;
    this->pade_sizreg_ssaice = pade_sizreg_ssaice;
    this->pade_sizreg_asyice = pade_sizreg_asyice;
    // Set default ice roughness - min values
    this->set_ice_roughness(1);
  }



  // Finalize
  void finalize() {
    icergh = 0;
    radliq_lwr = 0;
    radliq_upr = 0;
    radice_lwr = 0;
    radice_upr = 0;
    liq_nsteps = 0;
    ice_nsteps = 0;
    liq_step_size = 0;
    ice_step_size = 0;
    lut_extliq.deallocate();
    lut_ssaliq.deallocate();
    lut_asyliq.deallocate();
    lut_extice.deallocate();
    lut_ssaice.deallocate();
    lut_asyice.deallocate();
    pade_extliq.deallocate();
    pade_ssaliq.deallocate();
    pade_asyliq.deallocate();
    pade_extice.deallocate();
    pade_ssaice.deallocate();
    pade_asyice.deallocate();
    pade_sizreg_extliq.deallocate();
    pade_sizreg_ssaliq.deallocate();
    pade_sizreg_asyliq.deallocate();
    pade_sizreg_extice.deallocate();
    pade_sizreg_ssaice.deallocate();
    pade_sizreg_asyice.deallocate();

    // Base class finalize
    OpticalProps::finalize();
  }



  // Derive cloud optical properties from provided cloud physical properties
  // Compute single-scattering properties
  template <class T>  // T is a template for a child class of OpticalPropsArry
  void cloud_optics(const int ncol, const int nlay, real2d const &clwp, real2d const &ciwp, real2d const &reliq, real2d const &reice, T &optical_props) {
    using yakl::COLON;
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::any;
    using yakl::componentwise::operator>;
    using yakl::componentwise::operator<;
    using yakl::componentwise::operator&&;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    int nbnd = this->get_nband();
    // Error checking
    if (! (allocated(this->lut_extliq) || allocated(this->pade_extliq))) { stoprun("cloud optics: no data has been initialized"); }
    // Array sizes
    bool2d liqmsk("liqmsk",ncol, nlay);
    bool2d icemsk("icemsk",ncol, nlay);

    // Spectral consistency
    if (! this->bands_are_equal(optical_props)) { stoprun("cloud optics: optical properties don't have the same band structure"); }
    if (optical_props.get_nband() != optical_props.get_ngpt() ) {
      stoprun("cloud optics: optical properties must be requested by band not g-points");
    }

    // Cloud masks; don't need value re values if there's no cloud
    // do ilay = 1, nlay
    //   do icol = 1, ncol
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      liqmsk(icol,ilay) = clwp(icol,ilay) > 0.;
      icemsk(icol,ilay) = ciwp(icol,ilay) > 0.;
    }));

    #ifdef RRTMGP_EXPENSIVE_CHECKS
      // Particle size, liquid/ice water paths
      if ( any(liqmsk && (reliq < this->radliq_lwr)) || any(liqmsk && (reliq > this->radliq_upr)) ) {
        stoprun("cloud optics: liquid effective radius is out of bounds");
      }
      if ( any(icemsk && (reice < this->radice_lwr)) || any(icemsk && (reice > this->radice_upr)) ) {
        stoprun("cloud optics: ice effective radius is out of bounds");
      }
      if ( any(liqmsk && (clwp < 0.)) || any(icemsk && (ciwp < 0.)) ) {
        stoprun("cloud optics: negative clwp or ciwp where clouds are supposed to be");
      }
    #endif

    // The tables and Pade coefficients determing extinction coeffient, single-scattering albedo,
    //   and asymmetry parameter g as a function of effective raduis
    // We compute the optical depth tau (=exintinction coeff * condensed water path)
    //   and the products tau*ssa and tau*ssa*g for liquid and ice cloud separately.
    // These are used to determine the optical properties of ice and water cloud together.
    // We could compute the properties for liquid and ice separately and
    //    use ty_optical_props_arry.increment but this involves substantially more division.
    real3d ltau    ("ltau    ",size(clwp,1), size(clwp,2), this->get_nband());
    real3d ltaussa ("ltaussa ",size(clwp,1), size(clwp,2), this->get_nband());
    real3d ltaussag("ltaussag",size(clwp,1), size(clwp,2), this->get_nband());
    real3d itau    ("itau    ",size(clwp,1), size(clwp,2), this->get_nband());
    real3d itaussa ("itaussa ",size(clwp,1), size(clwp,2), this->get_nband());
    real3d itaussag("itaussag",size(clwp,1), size(clwp,2), this->get_nband());
    if (allocated(this->lut_extliq)) {
      // Liquid
      compute_all_from_table(ncol, nlay, nbnd, liqmsk, clwp, reliq, this->liq_nsteps,this->liq_step_size,this->radliq_lwr,
                             this->lut_extliq, this->lut_ssaliq, this->lut_asyliq, ltau, ltaussa, ltaussag);
      // Ice
      compute_all_from_table(ncol, nlay, nbnd, icemsk, ciwp, reice, this->ice_nsteps,this->ice_step_size,this->radice_lwr,
                             this->lut_extice.slice<2>(COLON,COLON,this->icergh),
                             this->lut_ssaice.slice<2>(COLON,COLON,this->icergh),
                             this->lut_asyice.slice<2>(COLON,COLON,this->icergh),
                             itau, itaussa, itaussag);
    } else {
      // Cloud optical properties from Pade coefficient method
      //   Hard coded assumptions: order of approximants, three size regimes
      int nsizereg = size(this->pade_extliq,2);
      compute_all_from_pade(ncol, nlay, nbnd, nsizereg, liqmsk, clwp, reliq,
                            2, 3, this->pade_sizreg_extliq, this->pade_extliq,
                            2, 2, this->pade_sizreg_ssaliq, this->pade_ssaliq,
                            2, 2, this->pade_sizreg_asyliq, this->pade_asyliq,
                            ltau, ltaussa, ltaussag);
      compute_all_from_pade(ncol, nlay, nbnd, nsizereg, icemsk, ciwp, reice,
                           2, 3, this->pade_sizreg_extice, this->pade_extice.slice<3>(COLON,COLON,COLON,this->icergh),
                           2, 2, this->pade_sizreg_ssaice, this->pade_ssaice.slice<3>(COLON,COLON,COLON,this->icergh),
                           2, 2, this->pade_sizreg_asyice, this->pade_asyice.slice<3>(COLON,COLON,COLON,this->icergh),
                           itau, itaussa, itaussag);
    }

    // Combine liquid and ice contributions into total cloud optical properties
    //   See also the increment routines in mo_optical_props_kernels
    combine( nbnd, nlay, ncol, ltau, itau, ltaussa, itaussa, ltaussag, itaussag, optical_props );
  }



  void combine( int nbnd, int nlay, int ncol, real3d const &ltau, real3d const &itau, real3d const &ltaussa, real3d const &itaussa,
                real3d const &ltaussag, real3d const &itaussag, OpticalProps1scl &optical_props ) {
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    auto &optical_props_tau = optical_props.tau;
    // do ibnd = 1, nbnd
    //   do ilay = 1, nlay
    //     do icol = 1,ncol
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlay,ncol) , YAKL_LAMBDA (int ibnd, int ilay, int icol) {
      // Absorption optical depth  = (1-ssa) * tau = tau - taussa
      optical_props_tau(icol,ilay,ibnd) = (ltau(icol,ilay,ibnd) - ltaussa(icol,ilay,ibnd)) +
                                          (itau(icol,ilay,ibnd) - itaussa(icol,ilay,ibnd));
    }));
  }



  void combine( int nbnd, int nlay, int ncol, real3d const &ltau, real3d const &itau, real3d const &ltaussa, real3d const &itaussa,
                real3d const &ltaussag, real3d const &itaussag, OpticalProps2str &optical_props ) {
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;
    using yakl::intrinsics::epsilon;

    auto &optical_props_g   = optical_props.g  ;
    auto &optical_props_ssa = optical_props.ssa;
    auto &optical_props_tau = optical_props.tau;
    // do ibnd = 1, nbnd
    //   do ilay = 1, nlay
    //     do icol = 1,ncol
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlay,ncol) , YAKL_LAMBDA (int ibnd, int ilay, int icol) {
      real tau    = ltau   (icol,ilay,ibnd) + itau   (icol,ilay,ibnd);
      real taussa = ltaussa(icol,ilay,ibnd) + itaussa(icol,ilay,ibnd);
      optical_props_g  (icol,ilay,ibnd) = (ltaussag(icol,ilay,ibnd) + itaussag(icol,ilay,ibnd)) / std::max(epsilon(tau), taussa);
      optical_props_ssa(icol,ilay,ibnd) = taussa / std::max(epsilon(tau), tau);
      optical_props_tau(icol,ilay,ibnd) = tau;
    }));
  }



  void set_ice_roughness(int icergh) {
    using yakl::intrinsics::allocated;

    if (! allocated(this->pade_extice) && ! allocated(this->lut_extice)) {
      stoprun("cloud_optics%set_ice_roughness(): can't set before initialization");
    }
    if (icergh < 1 || icergh > this->get_num_ice_roughness_types()) {
      stoprun("cloud optics: cloud ice surface roughness flag is out of bounds");
    }
    this->icergh = icergh;
  }



  int get_num_ice_roughness_types() {
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;

    int i;
    if (allocated(this->pade_extice)) { i = size(this->pade_extice,4); }
    if (allocated(this->lut_extice )) { i = size(this->lut_extice ,3); }
    return i;
  }



  real get_min_radius_liq() { return this->radliq_lwr; }



  real get_max_radius_liq() { return this->radliq_upr; }



  real get_min_radius_ice() { return this->radice_lwr; }



  real get_max_radius_ice() { return this->radice_upr; }



  // Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
  //   elements starting at "offset." The table's second dimension is band.
  // Returns 0 where the mask is false.
  // We could also try gather/scatter for efficiency
  void compute_all_from_table(int ncol, int nlay, int nbnd, bool2d const &mask, real2d const &lwp, real2d const &re,
                              int nsteps, real step_size, real offset, real2d const &tau_table, real2d const &ssa_table,
                              real2d const &asy_table, real3d &tau, real3d &taussa, real3d &taussag) {
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    // do ibnd = 1, nbnd
    //   do ilay = 1,nlay
    //     do icol = 1, ncol
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlay,ncol) , YAKL_LAMBDA (int ibnd, int ilay, int icol) {
      if (mask(icol,ilay)) {
        int index = std::min( floor( (re(icol,ilay) - offset) / step_size)+1, nsteps-1.);
        real fint = (re(icol,ilay) - offset)/step_size - (index-1);
        real t   = lwp(icol,ilay)    * (tau_table(index,  ibnd) + fint * (tau_table(index+1,ibnd) - tau_table(index,ibnd)));
        real ts  = t                 * (ssa_table(index,  ibnd) + fint * (ssa_table(index+1,ibnd) - ssa_table(index,ibnd)));
        taussag(icol,ilay,ibnd) = ts * (asy_table(index,  ibnd) + fint * (asy_table(index+1,ibnd) - asy_table(index,ibnd)));
        taussa (icol,ilay,ibnd) = ts;
        tau    (icol,ilay,ibnd) = t;
      } else {
        tau    (icol,ilay,ibnd) = 0;
        taussa (icol,ilay,ibnd) = 0;
        taussag(icol,ilay,ibnd) = 0;
      }
    }));
  }



  // Pade functions
  void compute_all_from_pade(int ncol, int nlay, int nbnd, int nsizes, bool2d const &mask, real2d const &lwp, real2d const &re,
                             int m_ext, int n_ext, real1d const &re_bounds_ext, real3d const &coeffs_ext,
                             int m_ssa, int n_ssa, real1d const &re_bounds_ssa, real3d const &coeffs_ssa,
                             int m_asy, int n_asy, real1d const &re_bounds_asy, real3d const &coeffs_asy,
                             real3d &tau, real3d &taussa, real3d &taussag) {
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    // do ibnd = 1, nbnd
    //   do ilay = 1, nlay
    //     do icol = 1, ncol
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlay,ncol) , YAKL_LAMBDA (int ibnd, int ilay, int icol) {
      if (mask(icol,ilay)) {
        int irad;
        // Finds index into size regime table
        // This works only if there are precisely three size regimes (four bounds) and it's
        //   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
        irad = std::min(floor((re(icol,ilay) - re_bounds_ext(2))/re_bounds_ext(3))+2, 3.);
        real t = lwp(icol,ilay) *         pade_eval(ibnd, nbnd, nsizes, m_ext, n_ext, irad, re(icol,ilay), coeffs_ext);

        irad = std::min(floor((re(icol,ilay) - re_bounds_ssa(2))/re_bounds_ssa(3))+2, 3.);
        // Pade approximants for co-albedo can sometimes be negative
        real ts = t * (1. - std::max(0., pade_eval(ibnd, nbnd, nsizes, m_ssa, n_ssa, irad, re(icol,ilay), coeffs_ssa)));

        irad = std::min(floor((re(icol,ilay) - re_bounds_asy(2))/re_bounds_asy(3))+2, 3.);
        taussag(icol,ilay,ibnd) = ts *    pade_eval(ibnd, nbnd, nsizes, m_asy, n_asy, irad, re(icol,ilay), coeffs_asy);

        taussa (icol,ilay,ibnd) = ts;
        tau    (icol,ilay,ibnd) = t;
      } else {
        tau    (icol,ilay,ibnd) = 0;
        taussa (icol,ilay,ibnd) = 0;
        taussag(icol,ilay,ibnd) = 0;
      }
    }));
  }



  // Evaluate Pade approximant of order [m/n]
  YAKL_INLINE real pade_eval(int iband, int nbnd, int nrads, int m, int n, int irad, real re, real3d const &pade_coeffs) {
    real denom = pade_coeffs(iband,irad,n+m);
    for (int i = n-1+m ; i >= 1+m ; i--) {
      denom = pade_coeffs(iband,irad,i)+re*denom;
    }
    denom = 1+re*denom;

    real numer = pade_coeffs(iband,irad,m);
    for (int i=m-1 ; i >= 1 ; i--) {
      numer = pade_coeffs(iband,irad,i)+re*numer;
    }
    numer = pade_coeffs(iband,irad,0)+re*numer;

    return numer/denom;
  }



  void print_norms(const bool print_prefix=false) const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;
    std::string prefix = print_prefix ? "JGFY" : "";

                                         std::cout << prefix << "name                   : " << name                    << "\n";
                                         std::cout << prefix << "icergh                 : " << icergh                  << "\n";
                                         std::cout << prefix << "radliq_lwr             : " << radliq_lwr              << "\n";
                                         std::cout << prefix << "radliq_upr             : " << radliq_upr              << "\n";
                                         std::cout << prefix << "radice_lwr             : " << radice_lwr              << "\n";
                                         std::cout << prefix << "radice_upr             : " << radice_upr              << "\n";
                                         std::cout << prefix << "liq_nsteps             : " << liq_nsteps              << "\n";
                                         std::cout << prefix << "ice_nsteps             : " << ice_nsteps              << "\n";
                                         std::cout << prefix << "liq_step_size          : " << liq_step_size           << "\n";
                                         std::cout << prefix << "ice_step_size          : " << ice_step_size           << "\n";
    if (allocated(lut_extliq        )) { std::cout << prefix << "sum(lut_extliq        ): " << sum(lut_extliq        ) << "\n"; }
    if (allocated(lut_ssaliq        )) { std::cout << prefix << "sum(lut_ssaliq        ): " << sum(lut_ssaliq        ) << "\n"; }
    if (allocated(lut_asyliq        )) { std::cout << prefix << "sum(lut_asyliq        ): " << sum(lut_asyliq        ) << "\n"; }
    if (allocated(lut_extice        )) { std::cout << prefix << "sum(lut_extice        ): " << sum(lut_extice        ) << "\n"; }
    if (allocated(lut_ssaice        )) { std::cout << prefix << "sum(lut_ssaice        ): " << sum(lut_ssaice        ) << "\n"; }
    if (allocated(lut_asyice        )) { std::cout << prefix << "sum(lut_asyice        ): " << sum(lut_asyice        ) << "\n"; }
    if (allocated(pade_extliq       )) { std::cout << prefix << "sum(pade_extliq       ): " << sum(pade_extliq       ) << "\n"; }
    if (allocated(pade_ssaliq       )) { std::cout << prefix << "sum(pade_ssaliq       ): " << sum(pade_ssaliq       ) << "\n"; }
    if (allocated(pade_asyliq       )) { std::cout << prefix << "sum(pade_asyliq       ): " << sum(pade_asyliq       ) << "\n"; }
    if (allocated(pade_extice       )) { std::cout << prefix << "sum(pade_extice       ): " << sum(pade_extice       ) << "\n"; }
    if (allocated(pade_ssaice       )) { std::cout << prefix << "sum(pade_ssaice       ): " << sum(pade_ssaice       ) << "\n"; }
    if (allocated(pade_asyice       )) { std::cout << prefix << "sum(pade_asyice       ): " << sum(pade_asyice       ) << "\n"; }
    if (allocated(pade_sizreg_extliq)) { std::cout << prefix << "sum(pade_sizreg_extliq): " << sum(pade_sizreg_extliq) << "\n"; }
    if (allocated(pade_sizreg_ssaliq)) { std::cout << prefix << "sum(pade_sizreg_ssaliq): " << sum(pade_sizreg_ssaliq) << "\n"; }
    if (allocated(pade_sizreg_asyliq)) { std::cout << prefix << "sum(pade_sizreg_asyliq): " << sum(pade_sizreg_asyliq) << "\n"; }
    if (allocated(pade_sizreg_extice)) { std::cout << prefix << "sum(pade_sizreg_extice): " << sum(pade_sizreg_extice) << "\n"; }
    if (allocated(pade_sizreg_ssaice)) { std::cout << prefix << "sum(pade_sizreg_ssaice): " << sum(pade_sizreg_ssaice) << "\n"; }
    if (allocated(pade_sizreg_asyice)) { std::cout << prefix << "sum(pade_sizreg_asyice): " << sum(pade_sizreg_asyice) << "\n"; }
    if (allocated(band2gpt          )) { std::cout << prefix << "band2gpt               : " << sum(band2gpt          ) << "\n"; }
    if (allocated(gpt2band          )) { std::cout << prefix << "gpt2band               : " << sum(gpt2band          ) << "\n"; }
    if (allocated(band_lims_wvn     )) { std::cout << prefix << "band_lims_wvn          : " << sum(band_lims_wvn     ) << "\n"; }
  }

};
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
template <typename RealT=double, typename LayoutT=Kokkos::LayoutLeft, typename DeviceT=DefaultDevice>
class CloudOpticsK : public OpticalPropsK<RealT, LayoutT, DeviceT> {
 public:

  using parent_t = OpticalPropsK<RealT, LayoutT, DeviceT>;

  template <typename T>
  using view_t = typename parent_t::template view_t<T>;

  template <typename T>
  using uview_t = conv::Unmanaged<view_t<T>>;

  using mdrp_t = conv::MDRP<LayoutT, DeviceT>;

  // Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  int icergh;  // (1 = none, 2 = medium, 3 = high)
  // Lookup table information
  // Upper and lower limits of the tables
  RealT radliq_lwr;
  RealT radliq_upr;
  RealT radice_lwr;
  RealT radice_upr;
  // How many steps in the table? (for convenience)
  int liq_nsteps;
  int ice_nsteps;
  // How big is each step in the table?
  RealT liq_step_size;
  RealT ice_step_size;
  // The tables themselves.
  view_t<RealT**> lut_extliq; // (nsize_liq, nbnd)
  view_t<RealT**> lut_ssaliq; // (nsize_liq, nbnd)
  view_t<RealT**> lut_asyliq; // (nsize_liq, nbnd)
  view_t<RealT***> lut_extice; // (nsize_ice, nbnd, nrghice)
  view_t<RealT***> lut_ssaice; // (nsize_ice, nbnd, nrghice)
  view_t<RealT***> lut_asyice; // (nsize_ice, nbnd, nrghice)
  // Pade approximant coefficients
  view_t<RealT***> pade_extliq; // (nbnd, nsizereg, ncoeff_ext)
  view_t<RealT***> pade_ssaliq; // (nbnd, nsizereg, ncoeff_ssa_g)
  view_t<RealT***> pade_asyliq; // (nbnd, nsizereg, ncoeff_ssa_g)
  view_t<RealT****> pade_extice; // (nbnd, nsizereg, ncoeff_ext, nrghice)
  view_t<RealT****> pade_ssaice; // (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  view_t<RealT****> pade_asyice; // (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  // Particle size regimes for Pade formulations
  view_t<RealT*> pade_sizreg_extliq; // (nbound)
  view_t<RealT*> pade_sizreg_ssaliq; // (nbound)
  view_t<RealT*> pade_sizreg_asyliq; // (nbound)
  view_t<RealT*> pade_sizreg_extice; // (nbound)
  view_t<RealT*> pade_sizreg_ssaice; // (nbound)
  view_t<RealT*> pade_sizreg_asyice; // (nbound)

  CloudOpticsK() = default;

  // Routines to load data needed for cloud optics calculations. Two routines: one to load
  //    lookup-tables and one for coefficients for Pade approximates
  template <typename BandLimsT, typename LutExtliqT, typename LutSsaliqT, typename LutAsyliqT,
            typename LutExticeT, typename LutSsaiceT, typename LutAsyiceT>
  void load(BandLimsT const &band_lims_wvn, RealT radliq_lwr, RealT radliq_upr, RealT radliq_fac, RealT radice_lwr, RealT radice_upr,
           RealT radice_fac, LutExtliqT const &lut_extliq, LutSsaliqT const &lut_ssaliq, LutAsyliqT const &lut_asyliq,
           LutExticeT const &lut_extice, LutSsaiceT const &lut_ssaice, LutAsyiceT const &lut_asyice) {
    using int2d_t = typename Kokkos::View<int**, LayoutT, DeviceT>;

    // Local variables
    this->init(band_lims_wvn, int2d_t(), "RRTMGP cloud optics");
    // LUT coefficient dimensions
    int nsize_liq = lut_extliq.extent(0);
    int nsize_ice = lut_extice.extent(0);
    int nbnd      = lut_extliq.extent(1);
    int nrghice   = lut_extice.extent(2);
    //
    // Error checking
    //   Can we check for consistency between table bounds and _fac?
    //
    if (nbnd != this->get_nband()) {
      stoprun("cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization");
    }
    if (lut_extice.extent(1) != nbnd) { stoprun("cloud_optics%init(): array lut_extice has the wrong number of bands"); }
    if (lut_ssaliq.extent(0) != nsize_liq || lut_ssaliq.extent(1) != nbnd) {
      stoprun("cloud_optics%init(): array lut_ssaliq isn't consistently sized");
    }
    if (lut_asyliq.extent(0) != nsize_liq || lut_asyliq.extent(1) != nbnd) {
      stoprun("cloud_optics%init(): array lut_asyliq isn't consistently sized");
    }
    if (lut_ssaice.extent(0) != nsize_ice || lut_ssaice.extent(1) != nbnd || lut_ssaice.extent(2) != nrghice) {
      stoprun("cloud_optics%init(): array lut_ssaice  isn't consistently sized");
    }
    if (lut_asyice.extent(0) != nsize_ice || lut_asyice.extent(1) != nbnd || lut_asyice.extent(2) != nrghice) {
      stoprun("cloud_optics%init(): array lut_asyice  isn't consistently sized");
    }

    this->liq_nsteps = nsize_liq;
    this->ice_nsteps = nsize_ice;
    this->liq_step_size = (radliq_upr - radliq_lwr) / (nsize_liq-1.);
    this->ice_step_size = (radice_upr - radice_lwr) / (nsize_ice-1.);
    // Load LUT constants
    this->radliq_lwr = radliq_lwr;
    this->radliq_upr = radliq_upr;
    this->radice_lwr = radice_lwr;
    this->radice_upr = radice_upr;
    // Load LUT coefficients
    this->lut_extliq = lut_extliq;
    this->lut_ssaliq = lut_ssaliq;
    this->lut_asyliq = lut_asyliq;
    this->lut_extice = lut_extice;
    this->lut_ssaice = lut_ssaice;
    this->lut_asyice = lut_asyice;
    // Set default ice roughness - min values
    this->set_ice_roughness(1);
  }

  // Cloud optics initialization function - Pade
  template <typename BandLimsT, typename PadeExtliqT, typename PadeSsaliqT, typename PadeAsyliqT,
            typename PadeExticeT, typename PadeSsaiceT, typename PadeAsyiceT, typename PadeSizeExtliqT,
            typename PadeSizeSsaliqT, typename PadeSizeAsyliqT, typename PadeSizeExticeT,
            typename PadeSizeSsaiceT, typename PadeSizeAsyiceT>
  void load(BandLimsT const &band_lims_wvn, PadeExtliqT const &pade_extliq, PadeSsaliqT const &pade_ssaliq, PadeAsyliqT const &pade_asyliq,
            PadeExticeT const &pade_extice, PadeSsaiceT const &pade_ssaice, PadeAsyiceT const &pade_asyice,
            PadeSizeExtliqT const &pade_sizreg_extliq, PadeSizeSsaliqT const &pade_sizreg_ssaliq, PadeSizeAsyliqT const &pade_sizreg_asyliq,
            PadeSizeExticeT const &pade_sizreg_extice, PadeSizeSsaiceT const &pade_sizreg_ssaice, PadeSizeAsyiceT const &pade_sizreg_asyice) {
    using int2d_t = typename Kokkos::View<int**, LayoutT, DeviceT>;

    // Pade coefficient dimensions
    int nbnd         = pade_extliq.extent(0);
    int nsizereg     = pade_extliq.extent(1);
    int ncoeff_ext   = pade_extliq.extent(2);
    int ncoeff_ssa_g = pade_ssaliq.extent(2);
    int nrghice      = pade_extice.extent(3);
    int nbound       = pade_sizreg_extliq.size();
    // The number of size regimes is assumed in the Pade evaluations
    if (nsizereg != 3) {
      stoprun("cloud optics: code assumes exactly three size regimes for Pade approximants but data is otherwise");
    }
    this->init(band_lims_wvn, int2d_t(), "RRTMGP cloud optics");
    // Error checking
    if (nbnd != this->get_nband()) {
      stoprun("cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization");
    }
    if (pade_ssaliq.extent(0) != nbnd || pade_ssaliq.extent(1) != nsizereg || pade_ssaliq.extent(2) != ncoeff_ssa_g) {
      stoprun("cloud_optics%init(): array pade_ssaliq isn't consistently sized");
    }
    if (pade_asyliq.extent(0) != nbnd || pade_asyliq.extent(1) != nsizereg || pade_asyliq.extent(2) != ncoeff_ssa_g) {
      stoprun("cloud_optics%init(): array pade_asyliq isn't consistently sized");
    }
    if (pade_extice.extent(0) != nbnd || pade_extice.extent(1) != nsizereg || pade_extice.extent(2) != ncoeff_ext || pade_extice.extent(3) != nrghice) {
      stoprun("cloud_optics%init(): array pade_extice isn't consistently sized");
    }
    if (pade_ssaice.extent(0) != nbnd || pade_ssaice.extent(1) != nsizereg || pade_ssaice.extent(2) != ncoeff_ssa_g || pade_ssaice.extent(3) != nrghice) {
      stoprun("cloud_optics%init(): array pade_ssaice isn't consistently sized");
    }
    if (pade_asyice.extent(0) != nbnd || pade_asyice.extent(1) != nsizereg || pade_asyice.extent(2) != ncoeff_ssa_g || pade_asyice.extent(3) != nrghice) {
      stoprun("cloud_optics%init(): array pade_asyice isn't consistently sized");
    }
    if (pade_sizreg_ssaliq.extent(0) != nbound || pade_sizreg_asyliq.extent(0) != nbound || pade_sizreg_extice.extent(0) != nbound ||
        pade_sizreg_ssaice.extent(0) != nbound || pade_sizreg_asyice.extent(0) != nbound ) {
      stoprun("cloud_optics%init(): one or more Pade size regime arrays are inconsistently sized");
    }
    if (nsizereg != 3) { stoprun("cloud_optics%init(): Expecting precisely three size regimes for Pade approximants"); }
    // Consistency among size regimes
    this->radliq_lwr = pade_sizreg_extliq(1);
    this->radliq_upr = pade_sizreg_extliq(nbound);
    this->radice_lwr = pade_sizreg_extice(1);
    this->radice_upr = pade_sizreg_extice(nbound);
    if (pade_sizreg_ssaliq(1) < this->radliq_lwr || pade_sizreg_asyliq(1) < this->radliq_lwr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have inconsistent lowest values");
    }
    if (pade_sizreg_ssaice(1) < this->radice_lwr || pade_sizreg_asyice(1) < this->radice_lwr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have inconsistent lower values");
    }
    if (pade_sizreg_ssaliq(nbound) > this->radliq_upr || pade_sizreg_asyliq(nbound) > this->radliq_upr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_upr");
    }
    if (pade_sizreg_ssaice(nbound) > this->radice_upr || pade_sizreg_asyice(nbound) > this->radice_upr) {
      stoprun("cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_upr");
    }
    // Load data
    this->pade_extliq        = pade_extliq       ;
    this->pade_ssaliq        = pade_ssaliq       ;
    this->pade_asyliq        = pade_asyliq       ;
    this->pade_extice        = pade_extice       ;
    this->pade_ssaice        = pade_ssaice       ;
    this->pade_asyice        = pade_asyice       ;
    this->pade_sizreg_extliq = pade_sizreg_extliq;
    this->pade_sizreg_ssaliq = pade_sizreg_ssaliq;
    this->pade_sizreg_asyliq = pade_sizreg_asyliq;
    this->pade_sizreg_extice = pade_sizreg_extice;
    this->pade_sizreg_ssaice = pade_sizreg_ssaice;
    this->pade_sizreg_asyice = pade_sizreg_asyice;
    // Set default ice roughness - min values
    this->set_ice_roughness(1);
  }

  // Finalize
  void finalize() {
    icergh = 0;
    radliq_lwr = 0;
    radliq_upr = 0;
    radice_lwr = 0;
    radice_upr = 0;
    liq_nsteps = 0;
    ice_nsteps = 0;
    liq_step_size = 0;
    ice_step_size = 0;
    lut_extliq = decltype(lut_extliq)();
    lut_ssaliq = decltype(lut_ssaliq)();
    lut_asyliq = decltype(lut_asyliq)();
    lut_extice = decltype(lut_extice)();
    lut_ssaice = decltype(lut_ssaice)();
    lut_asyice = decltype(lut_asyice)();
    pade_extliq = decltype(pade_extliq)();
    pade_ssaliq = decltype(pade_ssaliq)();
    pade_asyliq = decltype(pade_asyliq)();
    pade_extice = decltype(pade_extice)();
    pade_ssaice = decltype(pade_ssaice)();
    pade_asyice = decltype(pade_asyice)();
    pade_sizreg_extliq = decltype(pade_sizreg_extliq)();
    pade_sizreg_ssaliq = decltype(pade_sizreg_ssaliq)();
    pade_sizreg_asyliq = decltype(pade_sizreg_asyliq)();
    pade_sizreg_extice = decltype(pade_sizreg_extice)();
    pade_sizreg_ssaice = decltype(pade_sizreg_ssaice)();
    pade_sizreg_asyice = decltype(pade_sizreg_asyice)();

    // Base class finalize
    parent_t::finalize();
  }


  // Derive cloud optical properties from provided cloud physical properties
  // Compute single-scattering properties
  template <typename ClwpT, typename CiwpT, typename ReliqT, typename ReiceT, class OpticalPropsT>
  void cloud_optics(const int ncol, const int nlay,
                    ClwpT const &clwp, CiwpT const &ciwp, ReliqT const &reliq, ReiceT const &reice, OpticalPropsT &optical_props) {
    using pool = conv::MemPoolSingleton<RealT, LayoutT, DeviceT>;

    int nbnd = this->get_nband();
    // Error checking
    if (! (this->lut_extliq.is_allocated() || this->pade_extliq.is_allocated())) { stoprun("cloud optics: no data has been initialized"); }
    // Array sizes
    const int num_bools = ncol * nlay;
    bool* bdata         = pool::template alloc_raw<bool>(num_bools * 2), *bcurr = bdata;
    uview_t<bool**> liqmsk(bcurr, ncol, nlay); bcurr += num_bools;
    uview_t<bool**> icemsk(bcurr, ncol, nlay); bcurr += num_bools;

    // Spectral consistency
    if (! this->bands_are_equal(optical_props)) { stoprun("cloud optics: optical properties don't have the same band structure"); }
    if (optical_props.get_nband() != optical_props.get_ngpt() ) {
      stoprun("cloud optics: optical properties must be requested by band not g-points");
    }

    // Cloud masks; don't need value re values if there's no cloud
    // do ilay = 1, nlay
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol, nlay}), KOKKOS_LAMBDA (int icol, int ilay) {
      liqmsk(icol,ilay) = clwp(icol,ilay) > 0.;
      icemsk(icol,ilay) = ciwp(icol,ilay) > 0.;
    }));

    #ifdef RRTMGP_EXPENSIVE_CHECKS
      // Particle size, liquid/ice water paths
      if ( any(liqmsk && (reliq < this->radliq_lwr)) || any(liqmsk && (reliq > this->radliq_upr)) ) {
        stoprun("cloud optics: liquid effective radius is out of bounds");
      }
      if ( any(icemsk && (reice < this->radice_lwr)) || any(icemsk && (reice > this->radice_upr)) ) {
        stoprun("cloud optics: ice effective radius is out of bounds");
      }
      if ( any(liqmsk && (clwp < 0.)) || any(icemsk && (ciwp < 0.)) ) {
        stoprun("cloud optics: negative clwp or ciwp where clouds are supposed to be");
      }
    #endif

    // The tables and Pade coefficients determing extinction coeffient, single-scattering albedo,
    //   and asymmetry parameter g as a function of effective raduis
    // We compute the optical depth tau (=exintinction coeff * condensed water path)
    //   and the products tau*ssa and tau*ssa*g for liquid and ice cloud separately.
    // These are used to determine the optical properties of ice and water cloud together.
    // We could compute the properties for liquid and ice separately and
    //    use ty_optical_props_arry.increment but this involves substantially more division.
    const int num_tau = clwp.extent(0) *  clwp.extent(1) * this->get_nband();
    RealT* data        = pool::template alloc_raw<RealT>(num_tau * 6), *dcurr = data;
    uview_t<RealT***> ltau    (dcurr, clwp.extent(0), clwp.extent(1), this->get_nband()); dcurr += num_tau;
    uview_t<RealT***> ltaussa (dcurr, clwp.extent(0), clwp.extent(1), this->get_nband()); dcurr += num_tau;
    uview_t<RealT***> ltaussag(dcurr, clwp.extent(0), clwp.extent(1), this->get_nband()); dcurr += num_tau;
    uview_t<RealT***> itau    (dcurr, clwp.extent(0), clwp.extent(1), this->get_nband()); dcurr += num_tau;
    uview_t<RealT***> itaussa (dcurr, clwp.extent(0), clwp.extent(1), this->get_nband()); dcurr += num_tau;
    uview_t<RealT***> itaussag(dcurr, clwp.extent(0), clwp.extent(1), this->get_nband()); dcurr += num_tau;
    if (this->lut_extliq.is_allocated()) {
      // Liquid
      compute_all_from_table(ncol, nlay, nbnd, liqmsk, clwp, reliq, this->liq_nsteps,this->liq_step_size,this->radliq_lwr,
                             this->lut_extliq, this->lut_ssaliq, this->lut_asyliq, ltau, ltaussa, ltaussag);
      // Ice
      compute_all_from_table(ncol, nlay, nbnd, icemsk, ciwp, reice, this->ice_nsteps,this->ice_step_size,this->radice_lwr,
                             Kokkos::subview(this->lut_extice, Kokkos::ALL, Kokkos::ALL, this->icergh-1),
                             Kokkos::subview(this->lut_ssaice, Kokkos::ALL, Kokkos::ALL, this->icergh-1),
                             Kokkos::subview(this->lut_asyice, Kokkos::ALL, Kokkos::ALL, this->icergh-1),
                             itau, itaussa, itaussag);
    } else {
      // Cloud optical properties from Pade coefficient method
      //   Hard coded assumptions: order of approximants, three size regimes
      int nsizereg = this->pade_extliq.extent(1);
      compute_all_from_pade(ncol, nlay, nbnd, nsizereg, liqmsk, clwp, reliq,
                            2, 3, this->pade_sizreg_extliq, this->pade_extliq,
                            2, 2, this->pade_sizreg_ssaliq, this->pade_ssaliq,
                            2, 2, this->pade_sizreg_asyliq, this->pade_asyliq,
                            ltau, ltaussa, ltaussag);
      compute_all_from_pade(ncol, nlay, nbnd, nsizereg, icemsk, ciwp, reice,
                            2, 3, this->pade_sizreg_extice, Kokkos::subview(this->pade_extice, Kokkos::ALL,Kokkos::ALL,Kokkos::ALL,this->icergh-1),
                            2, 2, this->pade_sizreg_ssaice, Kokkos::subview(this->pade_ssaice, Kokkos::ALL,Kokkos::ALL,Kokkos::ALL,this->icergh-1),
                            2, 2, this->pade_sizreg_asyice, Kokkos::subview(this->pade_asyice, Kokkos::ALL,Kokkos::ALL,Kokkos::ALL,this->icergh-1),
                           itau, itaussa, itaussag);
    }

    // Combine liquid and ice contributions into total cloud optical properties
    //   See also the increment routines in mo_optical_props_kernels
    combine( nbnd, nlay, ncol, ltau, itau, ltaussa, itaussa, ltaussag, itaussag, optical_props );

    pool::dealloc(bdata, bcurr - bdata);
    pool::dealloc(data , dcurr - data);
  }

  template <typename LtauT, typename ItauT, typename LtaussaT,
            typename ItaussaT, typename LtaussagT, typename ItaussagT>
  void combine( int nbnd, int nlay, int ncol, LtauT const &ltau, ItauT const &itau, LtaussaT const &ltaussa, ItaussaT const &itaussa,
                LtaussagT const &ltaussag, ItaussagT const &itaussag, OpticalProps1sclK<RealT, LayoutT, DeviceT> &optical_props ) {
    auto &optical_props_tau = optical_props.tau;
    // do ibnd = 1, nbnd
    //   do ilay = 1, nlay
    //     do icol = 1,ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ncol,nlay,nbnd}), KOKKOS_LAMBDA (int icol, int ilay, int ibnd) {
      // Absorption optical depth  = (1-ssa) * tau = tau - taussa
      optical_props_tau(icol,ilay,ibnd) = (ltau(icol,ilay,ibnd) - ltaussa(icol,ilay,ibnd)) +
                                          (itau(icol,ilay,ibnd) - itaussa(icol,ilay,ibnd));
    }));
  }

  template <typename LtauT, typename ItauT, typename LtaussaT,
            typename ItaussaT, typename LtaussagT, typename ItaussagT>
  void combine( int nbnd, int nlay, int ncol, LtauT const &ltau, ItauT const &itau, LtaussaT const &ltaussa, ItaussaT const &itaussa,
                LtaussagT const &ltaussag, ItaussagT const &itaussag, OpticalProps2strK<RealT, LayoutT, DeviceT> &optical_props ) {
    auto &optical_props_g   = optical_props.g  ;
    auto &optical_props_ssa = optical_props.ssa;
    auto &optical_props_tau = optical_props.tau;
    // do ibnd = 1, nbnd
    //   do ilay = 1, nlay
    //     do icol = 1,ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ncol,nlay,nbnd}), KOKKOS_LAMBDA (int icol, int ilay, int ibnd) {
      RealT tau    = ltau   (icol,ilay,ibnd) + itau   (icol,ilay,ibnd);
      RealT taussa = ltaussa(icol,ilay,ibnd) + itaussa(icol,ilay,ibnd);
      optical_props_g  (icol,ilay,ibnd) = (ltaussag(icol,ilay,ibnd) + itaussag(icol,ilay,ibnd)) / Kokkos::fmax(conv::epsilon(tau), taussa);
      optical_props_ssa(icol,ilay,ibnd) = taussa / Kokkos::fmax(conv::epsilon(tau), tau);
      optical_props_tau(icol,ilay,ibnd) = tau;
    }));
  }

  void set_ice_roughness(int icergh) {
    if (! this->pade_extice.is_allocated() && ! this->lut_extice.is_allocated() ) {
      stoprun("cloud_optics%set_ice_roughness(): can't set before initialization");
    }
    if (icergh < 1 || icergh > this->get_num_ice_roughness_types()) {
      stoprun("cloud optics: cloud ice surface roughness flag is out of bounds");
    }
    this->icergh = icergh;
  }

  int get_num_ice_roughness_types() {
    int i;
    if (this->pade_extice.is_allocated()) { i = this->pade_extice.extent(3); }
    if (this->lut_extice.is_allocated()) { i = this->lut_extice.extent(2); }
    return i;
  }

  RealT get_min_radius_liq() { return this->radliq_lwr; }
  RealT get_max_radius_liq() { return this->radliq_upr; }
  RealT get_min_radius_ice() { return this->radice_lwr; }
  RealT get_max_radius_ice() { return this->radice_upr; }

  // Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
  //   elements starting at "offset." The table's second dimension is band.
  // Returns 0 where the mask is false.
  // We could also try gather/scatter for efficiency
  template <typename MaskT, typename LwpT, typename ReT, typename StrideViewT, typename TauT, typename TaussaT, typename TaussagT>
  void compute_all_from_table(int ncol, int nlay, int nbnd, MaskT const &mask, LwpT const &lwp, ReT const &re,
                              int nsteps, RealT step_size, RealT offset,
                              StrideViewT const &tau_table, StrideViewT const &ssa_table, StrideViewT const &asy_table,
                              TauT &tau, TaussaT &taussa, TaussagT &taussag) {
    // do ibnd = 1, nbnd
    //   do ilay = 1,nlay
    //     do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ncol,nlay,nbnd}), KOKKOS_LAMBDA (int icol, int ilay, int ibnd) {
      if (mask(icol,ilay)) {
        int index = Kokkos::fmin( floor( (re(icol,ilay) - offset) / step_size)+1, nsteps-1.) - 1;
        RealT fint = (re(icol,ilay) - offset)/step_size - index;
        RealT t   = lwp(icol,ilay)    * (tau_table(index,  ibnd) + fint * (tau_table(index+1,ibnd) - tau_table(index,ibnd)));
        RealT ts  = t                 * (ssa_table(index,  ibnd) + fint * (ssa_table(index+1,ibnd) - ssa_table(index,ibnd)));
        taussag(icol,ilay,ibnd) = ts * (asy_table(index,  ibnd) + fint * (asy_table(index+1,ibnd) - asy_table(index,ibnd)));
        taussa (icol,ilay,ibnd) = ts;
        tau    (icol,ilay,ibnd) = t;
      } else {
        tau    (icol,ilay,ibnd) = 0;
        taussa (icol,ilay,ibnd) = 0;
        taussag(icol,ilay,ibnd) = 0;
      }
    }));
  }

  // Pade functions
  template <typename MaskT, typename LwpT, typename ReT, typename ReBoundsExtT, typename ReBoundsSsaT, typename ReBoundsAsyT, typename StrideViewT,
            typename TauT, typename TaussaT, typename TaussagT>
  void compute_all_from_pade(int ncol, int nlay, int nbnd, int nsizes, MaskT const &mask, LwpT const &lwp, ReT const &re,
                             int m_ext, int n_ext, ReBoundsExtT const &re_bounds_ext, StrideViewT const &coeffs_ext,
                             int m_ssa, int n_ssa, ReBoundsSsaT const &re_bounds_ssa, StrideViewT const &coeffs_ssa,
                             int m_asy, int n_asy, ReBoundsAsyT const &re_bounds_asy, StrideViewT const &coeffs_asy,
                             TauT &tau, TaussaT &taussa, TaussagT &taussag) {
    // do ibnd = 1, nbnd
    //   do ilay = 1, nlay
    //     do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ncol, nlay, nbnd}), KOKKOS_LAMBDA (int icol, int ilay, int ibnd) {
      if (mask(icol,ilay)) {
        int irad;
        // Finds index into size regime table
        // This works only if there are precisely three size regimes (four bounds) and it's
        //   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
        irad = Kokkos::fmin(floor((re(icol,ilay) - re_bounds_ext(2))/re_bounds_ext(3))+2, 3.);
        RealT t = lwp(icol,ilay) *         pade_eval(ibnd, nbnd, nsizes, m_ext, n_ext, irad, re(icol,ilay), coeffs_ext);

        irad = Kokkos::fmin(floor((re(icol,ilay) - re_bounds_ssa(2))/re_bounds_ssa(3))+2, 3.);
        // Pade approximants for co-albedo can sometimes be negative
        RealT ts = t * (1. - Kokkos::fmax(0., pade_eval(ibnd, nbnd, nsizes, m_ssa, n_ssa, irad, re(icol,ilay), coeffs_ssa)));

        irad = Kokkos::fmin(floor((re(icol,ilay) - re_bounds_asy(2))/re_bounds_asy(3))+2, 3.);
        taussag(icol,ilay,ibnd) = ts *    pade_eval(ibnd, nbnd, nsizes, m_asy, n_asy, irad, re(icol,ilay), coeffs_asy);

        taussa (icol,ilay,ibnd) = ts;
        tau    (icol,ilay,ibnd) = t;
      } else {
        tau    (icol,ilay,ibnd) = 0;
        taussa (icol,ilay,ibnd) = 0;
        taussag(icol,ilay,ibnd) = 0;
      }
    }));
  }

  // Evaluate Pade approximant of order [m/n]
  template <typename PadeT>
  KOKKOS_INLINE_FUNCTION
  static RealT pade_eval(int iband, int nbnd, int nrads, int m, int n, int irad, RealT re, PadeT const &pade_coeffs) {
    RealT denom = pade_coeffs(iband,irad,n+m-1);
    for (int i = n-1+m ; i >= 1+m ; i--) {
      denom = pade_coeffs(iband,irad,i-1)+re*denom;
    }
    denom = 1+re*denom;

    RealT numer = pade_coeffs(iband,irad,m-1);
    for (int i=m-1 ; i >= 1 ; i--) {
      numer = pade_coeffs(iband,irad,i-1)+re*numer;
    }
    numer = pade_coeffs(iband,irad,0)+re*numer;

    return numer/denom;
  }

  void print_norms(const bool print_prefix=false) const {
    std::string prefix = print_prefix ? "JGFK" : "";
                                             std::cout << prefix << "name                   : " << this->name                          << "\n";
                                             std::cout << prefix << "icergh                 : " << this->icergh                        << "\n";
                                             std::cout << prefix << "radliq_lwr             : " << this->radliq_lwr                    << "\n";
                                             std::cout << prefix << "radliq_upr             : " << this->radliq_upr                    << "\n";
                                             std::cout << prefix << "radice_lwr             : " << this->radice_lwr                    << "\n";
                                             std::cout << prefix << "radice_upr             : " << this->radice_upr                    << "\n";
                                             std::cout << prefix << "liq_nsteps             : " << this->liq_nsteps                    << "\n";
                                             std::cout << prefix << "ice_nsteps             : " << this->ice_nsteps                    << "\n";
                                             std::cout << prefix << "liq_step_size          : " << this->liq_step_size                 << "\n";
                                             std::cout << prefix << "ice_step_size          : " << this->ice_step_size                 << "\n";
    if (this->lut_extliq.is_allocated()        ) { std::cout << prefix << "sum(lut_extliq        ): " << conv::sum(this->lut_extliq        ) << "\n"; }
    if (this->lut_ssaliq.is_allocated()        ) { std::cout << prefix << "sum(lut_ssaliq        ): " << conv::sum(this->lut_ssaliq        ) << "\n"; }
    if (this->lut_asyliq.is_allocated()        ) { std::cout << prefix << "sum(lut_asyliq        ): " << conv::sum(this->lut_asyliq        ) << "\n"; }
    if (this->lut_extice.is_allocated()        ) { std::cout << prefix << "sum(lut_extice        ): " << conv::sum(this->lut_extice        ) << "\n"; }
    if (this->lut_ssaice.is_allocated()        ) { std::cout << prefix << "sum(lut_ssaice        ): " << conv::sum(this->lut_ssaice        ) << "\n"; }
    if (this->lut_asyice.is_allocated()        ) { std::cout << prefix << "sum(lut_asyice        ): " << conv::sum(this->lut_asyice        ) << "\n"; }
    if (this->pade_extliq.is_allocated()       ) { std::cout << prefix << "sum(pade_extliq       ): " << conv::sum(this->pade_extliq       ) << "\n"; }
    if (this->pade_ssaliq.is_allocated()       ) { std::cout << prefix << "sum(pade_ssaliq       ): " << conv::sum(this->pade_ssaliq       ) << "\n"; }
    if (this->pade_asyliq.is_allocated()       ) { std::cout << prefix << "sum(pade_asyliq       ): " << conv::sum(this->pade_asyliq       ) << "\n"; }
    if (this->pade_extice.is_allocated()       ) { std::cout << prefix << "sum(pade_extice       ): " << conv::sum(this->pade_extice       ) << "\n"; }
    if (this->pade_ssaice.is_allocated()       ) { std::cout << prefix << "sum(pade_ssaice       ): " << conv::sum(this->pade_ssaice       ) << "\n"; }
    if (this->pade_asyice.is_allocated()       ) { std::cout << prefix << "sum(pade_asyice       ): " << conv::sum(this->pade_asyice       ) << "\n"; }
    if (this->pade_sizreg_extliq.is_allocated()) { std::cout << prefix << "sum(pade_sizreg_extliq): " << conv::sum(this->pade_sizreg_extliq) << "\n"; }
    if (this->pade_sizreg_ssaliq.is_allocated()) { std::cout << prefix << "sum(pade_sizreg_ssaliq): " << conv::sum(this->pade_sizreg_ssaliq) << "\n"; }
    if (this->pade_sizreg_asyliq.is_allocated()) { std::cout << prefix << "sum(pade_sizreg_asyliq): " << conv::sum(this->pade_sizreg_asyliq) << "\n"; }
    if (this->pade_sizreg_extice.is_allocated()) { std::cout << prefix << "sum(pade_sizreg_extice): " << conv::sum(this->pade_sizreg_extice) << "\n"; }
    if (this->pade_sizreg_ssaice.is_allocated()) { std::cout << prefix << "sum(pade_sizreg_ssaice): " << conv::sum(this->pade_sizreg_ssaice) << "\n"; }
    if (this->pade_sizreg_asyice.is_allocated()) { std::cout << prefix << "sum(pade_sizreg_asyice): " << conv::sum(this->pade_sizreg_asyice) << "\n"; }
    if (this->band2gpt.is_allocated()          ) { std::cout << prefix << "band2gpt               : " << conv::sum(this->band2gpt          ) << "\n"; }
    if (this->gpt2band.is_allocated()          ) { std::cout << prefix << "gpt2band               : " << conv::sum(this->gpt2band          ) << "\n"; }
    if (this->band_lims_wvn.is_allocated()     ) { std::cout << prefix << "band_lims_wvn          : " << conv::sum(this->band_lims_wvn     ) << "\n"; }
  }

#ifdef RRTMGP_ENABLE_YAKL
  void validate_kokkos(const CloudOptics& orig) const
  {
    parent_t::validate_kokkos(orig);

    RRT_REQUIRE(orig.icergh == this->icergh, "Bad icergh");
    RRT_REQUIRE(orig.radliq_lwr == this->radliq_lwr, "Bad radliq_lwr");
    RRT_REQUIRE(orig.radliq_upr == this->radliq_upr, "Bad radliq_upr");
    RRT_REQUIRE(orig.radice_lwr == this->radice_lwr, "Bad radice_lwr");
    RRT_REQUIRE(orig.radice_upr == this->radice_upr, "Bad radice_upr");
    RRT_REQUIRE(orig.liq_nsteps == this->liq_nsteps, "Bad liq_nsteps");
    RRT_REQUIRE(orig.ice_nsteps == this->ice_nsteps, "Bad ice_nsteps");
    RRT_REQUIRE(orig.liq_step_size == this->liq_step_size, "Bad liq_step_size");
    RRT_REQUIRE(orig.ice_step_size == this->ice_step_size, "Bad ice_step_size");

    conv::compare_yakl_to_kokkos(orig.lut_extliq, this->lut_extliq);
    conv::compare_yakl_to_kokkos(orig.lut_ssaliq, this->lut_ssaliq);
    conv::compare_yakl_to_kokkos(orig.lut_asyliq, this->lut_asyliq);
    conv::compare_yakl_to_kokkos(orig.lut_extice, this->lut_extice);
    conv::compare_yakl_to_kokkos(orig.lut_ssaice, this->lut_ssaice);
    conv::compare_yakl_to_kokkos(orig.lut_asyice, this->lut_asyice);
    conv::compare_yakl_to_kokkos(orig.pade_extliq, this->pade_extliq);
    conv::compare_yakl_to_kokkos(orig.pade_ssaliq, this->pade_ssaliq);
    conv::compare_yakl_to_kokkos(orig.pade_asyliq, this->pade_asyliq);
    conv::compare_yakl_to_kokkos(orig.pade_extice, this->pade_extice);
    conv::compare_yakl_to_kokkos(orig.pade_ssaice, this->pade_ssaice);
    conv::compare_yakl_to_kokkos(orig.pade_asyice, this->pade_asyice);
    conv::compare_yakl_to_kokkos(orig.pade_sizreg_extliq, this->pade_sizreg_extliq);
    conv::compare_yakl_to_kokkos(orig.pade_sizreg_ssaliq, this->pade_sizreg_ssaliq);
    conv::compare_yakl_to_kokkos(orig.pade_sizreg_asyliq, this->pade_sizreg_asyliq);
    conv::compare_yakl_to_kokkos(orig.pade_sizreg_extice, this->pade_sizreg_extice);
    conv::compare_yakl_to_kokkos(orig.pade_sizreg_ssaice, this->pade_sizreg_ssaice);
    conv::compare_yakl_to_kokkos(orig.pade_sizreg_asyice, this->pade_sizreg_asyice);
  }
#endif
};
#endif
