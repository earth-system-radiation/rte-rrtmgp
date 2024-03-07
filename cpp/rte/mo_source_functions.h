
#pragma once

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
// Encapsulate source function arrays for longwave/lw/internal sources
//    and shortwave/sw/external source.
//
// -------------------------------------------------------------------------------------------------

// Type for longwave sources: computed at layer center, at layer edges using
//   spectral mapping in each direction separately, and at the surface
// Not implementing get_subset because it isn't used
#ifdef RRTMGP_ENABLE_YAKL
class SourceFuncLW : public OpticalProps {
public:
  real3d lay_source;
  real3d lev_source_inc;
  real3d lev_source_dec;
  real2d sfc_source;


  bool is_allocated() const  { return this->is_initialized() && yakl::intrinsics::allocated(this->sfc_source); }


  void alloc(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("source_func_lw%alloc: not initialized so can't allocate"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("source_func_lw%alloc: must provide positive extents for ncol, nlay"); }
    int ngpt = this->get_ngpt();
    this->sfc_source     = real2d("sfc_source"    ,ncol,ngpt);
    this->lay_source     = real3d("lay_source"    ,ncol,nlay,ngpt);
    this->lev_source_inc = real3d("lev_source_inc",ncol,nlay,ngpt);
    this->lev_source_dec = real3d("lev_source_dec",ncol,nlay,ngpt);
  }


  void alloc(int ncol, int nlay, OpticalProps const &op) {
    if (! op.is_initialized()) { stoprun("source_func_lw::alloc: op not initialized"); }
    this->finalize();
    this->init(op);
    this->alloc(ncol,nlay);
  }


  void finalize() {
    this->lay_source     = real3d();
    this->lev_source_inc = real3d();
    this->lev_source_dec = real3d();
    this->sfc_source     = real2d();
    OpticalProps::finalize();
  }


  int get_ncol() const {
    if (this->is_allocated()) { return yakl::intrinsics::size(this->lay_source,1); } else { return 0; }
  }


  int get_nlay() const {
    if (this->is_allocated()) { return yakl::intrinsics::size(this->lay_source,2); } else { return 0; }
  }


  void print_norms(const bool print_prefix=false) const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

    std::string prefix = print_prefix ? "JGFY" : "";

                                     std::cout << prefix << "name          : " << name                << "\n";
    if (allocated(lay_source    )) { std::cout << prefix << "lay_source    : " << sum(lay_source    ) << "\n"; }
    if (allocated(lev_source_inc)) { std::cout << prefix << "lev_source_inc: " << sum(lev_source_inc) << "\n"; }
    if (allocated(lev_source_dec)) { std::cout << prefix << "lev_source_dec: " << sum(lev_source_dec) << "\n"; }
    if (allocated(sfc_source    )) { std::cout << prefix << "sfc_source    : " << sum(sfc_source    ) << "\n"; }
    if (allocated(band2gpt      )) { std::cout << prefix << "band2gpt      : " << sum(band2gpt      ) << "\n"; }
    if (allocated(gpt2band      )) { std::cout << prefix << "gpt2band      : " << sum(gpt2band      ) << "\n"; }
    if (allocated(band_lims_wvn )) { std::cout << prefix << "band_lims_wvn : " << sum(band_lims_wvn ) << "\n"; }
  }

};

// Type for shortave sources: top-of-domain spectrally-resolved flux
// Not implementing get_subset because it isn't used
class SourceFuncSW : public OpticalProps {
public:
  real2d toa_source;


  bool is_allocated() const { return this->is_initialized() && yakl::intrinsics::allocated(this->toa_source); }


  void alloc(int ncol) {
    if (! this->is_initialized()) { stoprun("source_func_sw%alloc: not initialized so can't allocate"); }
    if (ncol <= 0) { stoprun("source_func_sw%alloc: must provide positive extents for ncol"); }
    this->toa_source = real2d("toa_source",ncol,this->get_ngpt());
  }


  void alloc(int ncol, OpticalProps const &op) {
    if (! op.is_initialized()) { stoprun("source_func_sw::alloc: op not initialized"); }
    this->init(op);
    this->alloc(ncol);
  }


  void finalize() {
    this->toa_source = real2d();
    OpticalProps::finalize();
  }


  int get_ncol() const {
    if (this->is_allocated()) { return yakl::intrinsics::size(this->toa_source,1); } else { return 0; }
  }


  void print_norms(const bool print_prefix=false) const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

    std::string prefix = print_prefix ? "JGFY" : "";
                                     std::cout << prefix << "name          : " << name                << "\n";
    if (allocated(toa_source    )) { std::cout << prefix << "toa_source    : " << sum(toa_source    ) << "\n"; }
    if (allocated(band2gpt      )) { std::cout << prefix << "band2gpt      : " << sum(band2gpt      ) << "\n"; }
    if (allocated(gpt2band      )) { std::cout << prefix << "gpt2band      : " << sum(gpt2band      ) << "\n"; }
    if (allocated(band_lims_wvn )) { std::cout << prefix << "band_lims_wvn : " << sum(band_lims_wvn ) << "\n"; }
  }

};

#endif

#ifdef RRTMGP_ENABLE_KOKKOS
class SourceFuncLWK : public OpticalPropsK {
public:
  real3dk lay_source;
  real3dk lev_source_inc;
  real3dk lev_source_dec;
  real2dk sfc_source;


  bool is_allocated() const  { return this->is_initialized() && this->sfc_source.is_allocated(); }


  void alloc(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("source_func_lw%alloc: not initialized so can't allocate"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("source_func_lw%alloc: must provide positive extents for ncol, nlay"); }
    int ngpt = this->get_ngpt();
    this->sfc_source     = real2dk("sfc_source"    ,ncol,ngpt);
    this->lay_source     = real3dk("lay_source"    ,ncol,nlay,ngpt);
    this->lev_source_inc = real3dk("lev_source_inc",ncol,nlay,ngpt);
    this->lev_source_dec = real3dk("lev_source_dec",ncol,nlay,ngpt);
  }


  void alloc(int ncol, int nlay, OpticalPropsK const &op) {
    if (! op.is_initialized()) { stoprun("source_func_lw::alloc: op not initialized"); }
    this->finalize();
    this->init(op);
    this->alloc(ncol,nlay);
  }


  void finalize() {
    this->lay_source     = real3dk();
    this->lev_source_inc = real3dk();
    this->lev_source_dec = real3dk();
    this->sfc_source     = real2dk();
    OpticalPropsK::finalize();
  }


  int get_ncol() const {
    if (this->is_allocated()) { return this->lay_source.extent(0); } else { return 0; }
  }


  int get_nlay() const {
    if (this->is_allocated()) { return this->lay_source.extent(1); } else { return 0; }
  }

  void print_norms(const bool print_prefix=false) const {
    std::string prefix = print_prefix ? "JGFK" : "";
                                         std::cout << prefix << "name          : " << name                << "\n";
    if (lay_source.is_allocated()    ) { std::cout << prefix << "lay_source    : " << conv::sum(lay_source    ) << "\n"; }
    if (lev_source_inc.is_allocated()) { std::cout << prefix << "lev_source_inc: " << conv::sum(lev_source_inc) << "\n"; }
    if (lev_source_dec.is_allocated()) { std::cout << prefix << "lev_source_dec: " << conv::sum(lev_source_dec) << "\n"; }
    if (sfc_source.is_allocated()    ) { std::cout << prefix << "sfc_source    : " << conv::sum(sfc_source    ) << "\n"; }
    if (band2gpt.is_allocated()      ) { std::cout << prefix << "band2gpt      : " << conv::sum(band2gpt      ) << "\n"; }
    if (gpt2band.is_allocated()      ) { std::cout << prefix << "gpt2band      : " << conv::sum(gpt2band      ) << "\n"; }
    if (band_lims_wvn.is_allocated() ) { std::cout << prefix << "band_lims_wvn : " << conv::sum(band_lims_wvn ) << "\n"; }
  }

#ifdef RRTMGP_ENABLE_YAKL
  void validate_kokkos(const SourceFuncLW& orig)
  {
    OpticalPropsK::validate_kokkos(orig);
    conv::compare_yakl_to_kokkos(orig.lay_source, lay_source);
    conv::compare_yakl_to_kokkos(orig.lev_source_inc, lev_source_inc);
    conv::compare_yakl_to_kokkos(orig.lev_source_dec, lev_source_dec);
    conv::compare_yakl_to_kokkos(orig.sfc_source, sfc_source);
  }
#endif
};

// Type for shortave sources: top-of-domain spectrally-resolved flux
// Not implementing get_subset because it isn't used
class SourceFuncSWK : public OpticalPropsK {
public:
  real2dk toa_source;


  bool is_allocated() const { return this->is_initialized() && this->toa_source.is_allocated(); }


  void alloc(int ncol) {
    if (! this->is_initialized()) { stoprun("source_func_sw%alloc: not initialized so can't allocate"); }
    if (ncol <= 0) { stoprun("source_func_sw%alloc: must provide positive extents for ncol"); }
    this->toa_source = real2dk("toa_source",ncol,this->get_ngpt());
  }


  void alloc(int ncol, OpticalPropsK const &op) {
    if (! op.is_initialized()) { stoprun("source_func_sw::alloc: op not initialized"); }
    this->init(op);
    this->alloc(ncol);
  }


  void finalize() {
    this->toa_source = real2dk();
    OpticalPropsK::finalize();
  }


  int get_ncol() const {
    if (this->is_allocated()) { return this->toa_source.extent(0); } else { return 0; }
  }


  void print_norms(const bool print_prefix=false) const {
    std::string prefix = print_prefix ? "JGFK" : "";
                                         std::cout << prefix << "name          : " << name                << "\n";
    if (toa_source.is_allocated()    ) { std::cout << prefix << "toa_source    : " << conv::sum(toa_source    ) << "\n"; }
    if (band2gpt.is_allocated()      ) { std::cout << prefix << "band2gpt      : " << conv::sum(band2gpt      ) << "\n"; }
    if (gpt2band.is_allocated()      ) { std::cout << prefix << "gpt2band      : " << conv::sum(gpt2band      ) << "\n"; }
    if (band_lims_wvn.is_allocated() ) { std::cout << prefix << "band_lims_wvn : " << conv::sum(band_lims_wvn ) << "\n"; }
  }

#ifdef RRTMGP_ENABLE_YAKL
  void validate_kokkos(const SourceFuncSW& orig)
  {
    OpticalPropsK::validate_kokkos(orig);
    conv::compare_yakl_to_kokkos(orig.toa_source, toa_source);
  }
#endif

};
#endif
