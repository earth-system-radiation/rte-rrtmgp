
#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props_kernels.h"

// Base class for optical properties
//   Describes the spectral discretization including the wavenumber limits
//   of each band (spectral region) and the mapping between g-points and bands
class OpticalProps {
public:
  int2d  band2gpt;       // (begin g-point, end g-point) = band2gpt(2,band)
  int1d  gpt2band;       // band = gpt2band(g-point)
  int    ngpt;
  real2d band_lims_wvn;  // (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  std::string name;


  // Base class: Initialization
  //   Values are assumed to be defined in bands a mapping between bands and g-points is provided
  void init( real2d const &band_lims_wvn , int2d const &band_lims_gpt=int2d() , std::string name="" ) {
    using yakl::intrinsics::size;
    using yakl::intrinsics::any;
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::maxval;
    using yakl::componentwise::operator<;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    int2d band_lims_gpt_lcl("band_lims_gpt_lcl",2,size(band_lims_wvn,2));
    if (size(band_lims_wvn,1) != 2) { stoprun("optical_props::init(): band_lims_wvn 1st dim should be 2"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(band_lims_wvn < 0._wp)) { stoprun("optical_props::init(): band_lims_wvn has values <  0."); }
    #endif
    if (allocated(band_lims_gpt)) {
      if (size(band_lims_gpt,2) != size(band_lims_wvn,2)) {
        stoprun("optical_props::init(): band_lims_gpt size inconsistent with band_lims_wvn");
      }
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (any(band_lims_gpt < 1) ) { stoprun("optical_props::init(): band_lims_gpt has values < 1"); }
      #endif
      // for (int j=1; j <= size(band_lims_gpt,2); j++) {
      //   for (int i=1; i <= size(band_lims_gpt,1); i++) {
      parallel_for( SimpleBounds<2>(size(band_lims_gpt,2),size(band_lims_gpt,1)) , YAKL_LAMBDA (int j, int i) {
        band_lims_gpt_lcl(i,j) = band_lims_gpt(i,j);
      });
    } else {
      // Assume that values are defined by band, one g-point per band
      // for (int iband = 1; iband <= size(band_lims_wvn, 2); iband++) {
      parallel_for( SimpleBounds<1>(size(band_lims_wvn, 2)) , YAKL_LAMBDA (int iband) {
        band_lims_gpt_lcl(2,iband) = iband;
        band_lims_gpt_lcl(1,iband) = iband;
      });
    }
    // Assignment
    this->band2gpt      = band_lims_gpt_lcl;
    this->band_lims_wvn = band_lims_wvn;
    this->name          = name;
    this->ngpt          = maxval(this->band2gpt);

    // Make a map between g-points and bands
    //   Efficient only when g-point indexes start at 1 and are contiguous.
    this->gpt2band = int1d("gpt2band",maxval(band_lims_gpt_lcl));
    // TODO: I didn't want to bother with race conditions at the moment, so it's an entirely serialized kernel for now
    auto &this_gpt2band = this->gpt2band;
    parallel_for( SimpleBounds<1>(1) , YAKL_LAMBDA (int dummy) {
      for (int iband=1; iband <= size(band_lims_gpt_lcl,2); iband++) {
        for (int i=band_lims_gpt_lcl(1,iband); i <= band_lims_gpt_lcl(2,iband); i++) {
          this_gpt2band(i) = iband;
        }
      }
    });
  }


  void init(OpticalProps const &in) {
    if ( ! in.is_initialized() ) {
      stoprun("optical_props::init(): can't initialize based on un-initialized input");
    } else {
      this->init( in.get_band_lims_wavenumber() , in.get_band_lims_gpoint() );
    }
  }


  bool is_initialized() const { return yakl::intrinsics::allocated(this->band2gpt); }


  // Base class: finalize (deallocate memory)
  void finalize() {
    this->band2gpt      = int2d();
    this->gpt2band      = int1d();
    this->band_lims_wvn = real2d();
    this->name          = "";
  }


  // Number of bands
  int get_nband() const {
    using yakl::intrinsics::size;
    if (this->is_initialized()) { return yakl::intrinsics::size(this->band2gpt,2); }
    return 0;
  }


  // Number of g-points
  int get_ngpt() const {
    if (this->is_initialized()) { return this->ngpt; }
    return 0;
  }


  // Bands for all the g-points at once;  dimension (ngpt)
  int1d get_gpoint_bands() const { return gpt2band; }


  // First and last g-point of a specific band
  int1d convert_band2gpt(int band) const {
    int1d ret("band2gpt",2);
    if (this->is_initialized()) {
      ret(1) = this->band2gpt(1,band);
      ret(2) = this->band2gpt(2,band);
    } else {
      ret(1) = 0;
      ret(2) = 0;
    }
    return ret;
  }


  // Band associated with a specific g-point
  int convert_gpt2band(int gpt) const {
    if (this->is_initialized()) { return this->gpt2band(gpt); }
    return 0;
  }


  // The first and last g-point of all bands at once;  dimension (2, nbands)
  int2d get_band_lims_gpoint() const { return this->band2gpt; }


  // Lower and upper wavenumber of all bands
  // (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  real2d get_band_lims_wavenumber() const { return this->band_lims_wvn; }


  // Lower and upper wavelength of all bands
  real2d get_band_lims_wavelength() const {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    real2d ret("band_lim_wavelength",size(band_lims_wvn,1),size(band_lims_wvn,2));
    // for (int j = 1; j <= size(band_lims_wvn,2); j++) {
    //   for (int i = 1; i <= size(band_lims_wvn,1); i++) {
    auto &this_band_lims_wvn = this->band_lims_wvn;
    parallel_for( SimpleBounds<2>( size(band_lims_wvn,2) , size(band_lims_wvn,1) ) , YAKL_LAMBDA (int j, int i) {
      if (this->is_initialized()) {
        ret(i,j) = 1._wp / this_band_lims_wvn(i,j);
      } else {
        ret(i,j) = 0._wp;
      }
    });
    return ret;
  }


  // Are the bands of two objects the same? (same number, same wavelength limits)
  bool bands_are_equal(OpticalProps const &rhs) const {
    using yakl::intrinsics::size;
    using yakl::intrinsics::epsilon;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    // if ( this->get_nband() != rhs.get_nband() || this->get_nband() == 0) { return false; }
    // yakl::ScalarLiveOut<bool> ret(true);
    // // for (int j=1 ; j <= size(this->band_lims_wvn,2); j++) {
    // //   for (int i=1 ; i <= size(this->band_lims_wvn,1); i++) {
    // auto &this_band_lims_wvn = this->band_lims_wvn;
    // auto &rhs_band_lims_wvn  = rhs.band_lims_wvn;
    // parallel_for( Bounds<2>( size(this->band_lims_wvn,2) , size(this->band_lims_wvn,1) ) , YAKL_LAMBDA (int j, int i) {
    //   if ( abs( this_band_lims_wvn(i,j) - rhs_band_lims_wvn(i,j) ) > 5*epsilon(this_band_lims_wvn) ) {
    //     ret = false;
    //   }
    // });
    // return ret.hostRead();


    // This is working around an issue that arises in E3SM's rrtmgpxx integration.
    // Previously the code failed in the creation of the ScalarLiveOut variable, but only for higher optimizations
    bool ret = true;
    auto this_band_lims_wvn = this->band_lims_wvn.createHostCopy();
    auto rhs_band_lims_wvn  = rhs.band_lims_wvn  .createHostCopy();
    for (int j=1 ; j <= size(this->band_lims_wvn,2); j++) {
      for (int i=1 ; i <= size(this->band_lims_wvn,1); i++) {
        if ( abs( this_band_lims_wvn(i,j) - rhs_band_lims_wvn(i,j) ) > 5*epsilon(this_band_lims_wvn) ) {
          ret = false;
        }
      }
    }
    return ret;
  }


  // Is the g-point structure of two objects the same?
  //   (same bands, same number of g-points, same mapping between bands and g-points)
  bool gpoints_are_equal(OpticalProps const &rhs) const {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    if ( ! this->bands_are_equal(rhs) || this->get_ngpt() != rhs.get_ngpt() ) { return false; }
    // yakl::ScalarLiveOut<bool> ret(true);
    // // for (int i=1; i <= size(this->gpt2bnd,1); i++) {
    // auto &this_gpt2band = this->gpt2band;
    // auto &rhs_gpt2band  = rhs.gpt2band;
    // parallel_for( Bounds<1>(size(this->gpt2band,1)) , YAKL_LAMBDA (int i) {
    //   if ( this_gpt2band(i) != rhs_gpt2band(i) ) { ret = false; }
    // });
    // return ret.hostRead();


    // This is working around an issue that arises in E3SM's rrtmgpxx integration.
    // Previously the code failed in the creation of the ScalarLiveOut variable, but only for higher optimizations
    bool ret = true;
    auto this_gpt2band = this->gpt2band.createHostCopy();
    auto rhs_gpt2band  = rhs.gpt2band  .createHostCopy();
    for (int i=1; i <= size(this->gpt2band,1); i++) {
      if ( this_gpt2band(i) != rhs_gpt2band(i) ) { ret = false; }
    }
    return ret;
  }


  // Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  // real1d expand(real1d const &arr_in) const {
  //   real1d ret("arr_out",size(this->gpt2band,1));
  //   // do iband=1,this->get_nband()
  //   // TODO: I don't know if this needs to be serialize or not at first glance. Need to look at it more.
  //   auto &this_band2gpt = this->gpt2band;
  //   int nband = get_nband();
  //   parallel_for( Bounds<1>(1) , YAKL_LAMBDA (int dummy) {
  //     for (int iband = 1 ; iband <= nband ; iband++) {
  //       for (int i=this_band2gpt(1,iband) ; i <= this_band2gpt(2,iband) ; i++) {
  //         ret(i) = arr_in(iband);
  //       }
  //     }
  //   });
  //   return ret;
  // }


  void set_name( std::string name ) { this->name = name; }


  std::string get_name() const { return this->name; }
};



class OpticalPropsArry : public OpticalProps {
public:
  real3d tau; // optical depth (ncol, nlay, ngpt)

  int get_ncol() const { if (yakl::intrinsics::allocated(tau)) { return yakl::intrinsics::size(this->tau,1); } else { return 0; } }
  int get_nlay() const { if (yakl::intrinsics::allocated(tau)) { return yakl::intrinsics::size(this->tau,2); } else { return 0; } }
};



// We need to know about 2str's existence because it is referenced in 1scl
class OpticalProps2str;



// Not implementing get_subset because it isn't used
class OpticalProps1scl : public OpticalPropsArry {
public:
  void validate() const {
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::any;
    using yakl::componentwise::operator<;

    if (! allocated(this->tau)) { stoprun("validate: tau not allocated/initialized"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(this->tau < 0)) { stoprun("validate: tau values out of range"); }
    #endif
  }


  void delta_scale(real3d const &dummy) const { }


  void alloc_1scl(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("OpticalProps1scl::alloc_1scl: spectral discretization hasn't been provided"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("OpticalProps1scl::alloc_1scl: must provide > 0 extents for ncol, nlay"); }
    this->tau = real3d("tau",ncol,nlay,this->get_ngpt());
    memset(tau,0._wp);
  }


  // Initialization by specifying band limits and possibly g-point/band mapping
  void alloc_1scl(int ncol, int nlay, real2d const &band_lims_wvn, int2d const &band_lims_gpt=int2d(), std::string name="") {
    this->init(band_lims_wvn, band_lims_gpt, name);
    this->alloc_1scl(ncol, nlay);
  }


  void alloc_1scl(int ncol, int nlay, OpticalProps const &opIn, std::string name="") {
    if (this->is_initialized()) { this->finalize(); }
    this->init(opIn.get_band_lims_wavenumber(), opIn.get_band_lims_gpoint(), name);
    this->alloc_1scl(ncol, nlay);
  }


  void increment(OpticalProps1scl &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol();
    int nlay = that.get_nlay();
    int ngpt = that.get_ngpt();
    if (this->gpoints_are_equal(that)) {
      increment_1scalar_by_1scalar(ncol, nlay, ngpt, that.tau, this->tau);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt, that.tau, this->tau, that.get_nband(), that.get_band_lims_gpoint());
    }
  }


  // Implemented later because OpticalProps2str hasn't been created yet
  inline void increment(OpticalProps2str &that);


  void print_norms() const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

                                    std::cout << "name         : " << name               << "\n";
    if (allocated(band2gpt     )) { std::cout << "band2gpt     : " << sum(band2gpt     ) << "\n"; }
    if (allocated(gpt2band     )) { std::cout << "gpt2band     : " << sum(gpt2band     ) << "\n"; }
    if (allocated(band_lims_wvn)) { std::cout << "band_lims_wvn: " << sum(band_lims_wvn) << "\n"; }
    if (allocated(tau          )) { std::cout << "tau          : " << sum(tau          ) << "\n"; }
  }

};



// Not implementing get_subset because it isn't used
class OpticalProps2str : public OpticalPropsArry {
public:
  real3d ssa; // single-scattering albedo (ncol, nlay, ngpt)
  real3d g;   // asymmetry parameter (ncol, nlay, ngpt)


  void validate() const {
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::any;
    using yakl::componentwise::operator<;
    using yakl::componentwise::operator>;

    if ( ! allocated(this->tau) || ! allocated(this->ssa) || ! allocated(this->g) ) {
      stoprun("validate: arrays not allocated/initialized");
    }
    int d1 = size(this->tau,1);
    int d2 = size(this->tau,2);
    int d3 = size(this->tau,3);
    if ( d1 != size(this->ssa,1) || d2 != size(this->ssa,2) || d3 != size(this->ssa,3) || 
         d1 != size(this->g  ,1) || d2 != size(this->g  ,2) || d3 != size(this->g  ,3) ) {
      stoprun("validate: arrays not sized consistently");
    }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(this->tau <  0)                      ) { stoprun("validate: tau values out of range"); }
      if (any(this->ssa <  0) || any(this->ssa > 1)) { stoprun("validate: ssa values out of range"); }
      if (any(this->g   < -1) || any(this->g   > 1)) { stoprun("validate: g   values out of range"); }
    #endif
  }


  void delta_scale(real3d const &forward=real3d()) {
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::any;
    using yakl::componentwise::operator<;
    using yakl::componentwise::operator>;

    // Forward scattering fraction; g**2 if not provided
    int ncol = this->get_ncol();
    int nlay = this->get_nlay();
    int ngpt = this->get_ngpt();
    if (allocated(forward)) {
      if (size(forward,1) != ncol || size(forward,2) != nlay || size(forward,3) != ngpt) {
        stoprun("delta_scale: dimension of 'forward' don't match optical properties arrays");
      }
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (any(forward < 0) || any(forward > 1)) { stoprun("delta_scale: values of 'forward' out of bounds [0,1]"); }
      #endif
      delta_scale_2str_kernel(ncol, nlay, ngpt, this->tau, this->ssa, this->g, forward);
    } else {
      delta_scale_2str_kernel(ncol, nlay, ngpt, this->tau, this->ssa, this->g);
    }
  }


  void alloc_2str(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("optical_props::alloc: spectral discretization hasn't been provided"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("optical_props::alloc: must provide positive extents for ncol, nlay"); }
    this->tau = real3d("tau",ncol,nlay,this->get_ngpt());
    this->ssa = real3d("ssa",ncol,nlay,this->get_ngpt());
    this->g   = real3d("g  ",ncol,nlay,this->get_ngpt());
    memset(tau,0._wp);
    memset(ssa,0._wp);
    memset(g  ,0._wp);
  }


  void alloc_2str(int ncol, int nlay, real2d const &band_lims_wvn, int2d const &band_lims_gpt=int2d(), std::string name="") {
    this->init(band_lims_wvn, band_lims_gpt, name);
    this->alloc_2str(ncol, nlay);
  }


  void alloc_2str(int ncol, int nlay, OpticalProps const &opIn, std::string name="") {
    if (this->is_initialized()) { this->finalize(); }
    this->init(opIn.get_band_lims_wavenumber(), opIn.get_band_lims_gpoint(), name);
    this->alloc_2str(ncol, nlay);
  }


  void increment(OpticalProps1scl &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol();
    int nlay = that.get_nlay();
    int ngpt = that.get_ngpt();
    if (this->gpoints_are_equal(that)) {
      increment_1scalar_by_2stream(ncol, nlay, ngpt, that.tau, this->tau, this->ssa);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt, that.tau, this->tau, this->ssa, that.get_nband(), that.get_band_lims_gpoint());
    }
  }


  void increment(OpticalProps2str &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol();
    int nlay = that.get_nlay();
    int ngpt = that.get_ngpt();
    if (this->gpoints_are_equal(that)) {
      increment_2stream_by_2stream(ncol, nlay, ngpt, that.tau, that.ssa, that.g, this->tau, this->ssa, this->g);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt, that.tau, that.ssa, that.g, this->tau, this->ssa, this->g, that.get_nband(), that.get_band_lims_gpoint());
    }
  }


  void print_norms() const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

                                    std::cout << "name         : " << name               << "\n";
    if (allocated(band2gpt     )) { std::cout << "band2gpt     : " << sum(band2gpt     ) << "\n"; }
    if (allocated(gpt2band     )) { std::cout << "gpt2band     : " << sum(gpt2band     ) << "\n"; }
    if (allocated(band_lims_wvn)) { std::cout << "band_lims_wvn: " << sum(band_lims_wvn) << "\n"; }
    if (allocated(tau          )) { std::cout << "tau          : " << sum(tau          ) << "\n"; }
    if (allocated(ssa          )) { std::cout << "ssa          : " << sum(ssa          ) << "\n"; }
    if (allocated(g            )) { std::cout << "g            : " << sum(g            ) << "\n"; }  
  }

};



inline void OpticalProps1scl::increment(OpticalProps2str &that) {
  if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
  int ncol = that.get_ncol();
  int nlay = that.get_nlay();
  int ngpt = that.get_ngpt();
  if (this->gpoints_are_equal(that)) {
    increment_2stream_by_1scalar(ncol, nlay, ngpt, that.tau, that.ssa, this->tau);
  } else {
    if (this->get_ngpt() != that.get_nband()) {
      stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
    }
    inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt, that.tau, that.ssa, this->tau, that.get_nband(), that.get_band_lims_gpoint());
  }
}



