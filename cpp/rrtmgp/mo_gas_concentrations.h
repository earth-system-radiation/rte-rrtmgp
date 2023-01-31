
#pragma once

#include "rrtmgp_const.h"
#include "mo_rrtmgp_util_string.h"

// This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
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
//  Stores, sets, and gets volume mixing ratios for a set of gasses
//
//  Volume Mixing Ratios passed to set_vmr are expected to be on the device
//
//  All loops involving strings are on the host. All loops involving concs are on the device
//
// -------------------------------------------------------------------------------------------------


class GasConcs {
public:
  static int constexpr GAS_NOT_IN_LIST = -1;
  string1d gas_name;  // List of gas names defined upon init
  real3d   concs;     // List of gas concentrations (ngas,ncol,nlay)
  int      ncol;
  int      nlay;
  int      ngas;


  GasConcs() {
    ncol = 0;
    nlay = 0;
    ngas = 0;
  }


  ~GasConcs() {
    reset();
  }


  void reset() {
    gas_name = string1d();  // Dealloc
    concs    = real3d();    // Dealloc
    ncol = 0;
    nlay = 0;
    ngas = 0;
  }


  void init(string1d const &gas_names , int ncol , int nlay) {
    using yakl::intrinsics::size;

    this->reset();
    this->ngas = size(gas_names,1);
    this->ncol = ncol;
    this->nlay = nlay;

    // Transform gas names to lower case, check for empty strings, check for duplicates
    for (int i=1; i<=ngas; i++) {
      gas_names(i) = lower_case( gas_names(i) );
      // Empty string
      if (gas_names(i) == "") { stoprun("ERROR: GasConcs::init(): must provide non-empty gas names"); }
      // Duplicate gas name
      for (int j=i+1; j<=ngas; j++) {
        if ( gas_names(i) == gas_names(j) ) { stoprun("GasConcs::init(): duplicate gas names aren't allowed"); }
      }
    }
    
    // Allocate
    this->gas_name = string1d("gas_name",ngas);
    this->concs    = real3d  ("concs"   ,ncol,nlay,ngas);

    // Assign gas names
    for (int i=1; i<=ngas; i++) {
      this->gas_name(i) = gas_names(i);
    }
  }


  // Set concentration as a scalar copied to every column and level
  void set_vmr(std::string gas, real w) {
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name was not provided at initialization");
    }
    if (w < 0. || w > 1.) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    YAKL_SCOPE( this_concs , this->concs );
    // for (int ilay=1; ilay<=this->nlay; ilay++) {
    //   for (int icol=1; icol<=this->ncol; icol++) {
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      this_concs(icol,ilay,igas) = w;
    });
  }


  // Set concentration as a single column copied to all other columns
  // w is expected to be in device memory
  void set_vmr(std::string gas, real1d const &w) {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    if (size(w,1) != this->nlay) { stoprun("GasConcs::set_vmr: different dimension (nlay)"); }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name not provided at initialization");
    }
    // Check for bad values in w
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      yakl::ScalarLiveOut<bool> badVal(false); // Scalar that must exist in device memory (equiv: bool badVal = false;)
      // for (int i=1; i<=size(w,1); i++) {
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>(size(w,1)) , YAKL_LAMBDA (int i) {
        if (w(i) < 0. || w(i) > 1.) { badVal = true; }
      });
      if (badVal.hostRead()) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    #endif
    YAKL_SCOPE( this_concs , this->concs );
    // for (int ilay=1; ilay<=this->nlay; ilay++) {
    //   for (int icol=1; icol<=this->ncol; icol++) {
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      this_concs(icol,ilay,igas) = w(ilay);
    });
  }
  

  // Set concentration as a 2-D field of columns and levels
  // w is expected to be in device memory
  void set_vmr(std::string gas, real2d const &w) {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    if (size(w,1) != this->ncol) { stoprun("GasConcs::set_vmr: different dimension (ncol)" ); }
    if (size(w,2) != this->nlay) { stoprun("GasConcs::set_vmr: different dimension (nlay)" ); }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name not provided at initialization" );
    }
    // Check for bad values in w
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      yakl::ScalarLiveOut<bool> badVal(false); // Scalar that must exist in device memory (equiv: bool badVal = false;)
      // for (int j=1; j<=size(w,2); j++) {
      //   for (int i=1; i<=size(w,1); i++) {
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(size(w,2),size(w,1)) , YAKL_LAMBDA (int j, int i) {
        if (w(i,j) < 0. || w(i,j) > 1.) { badVal = true;}
      });
      if (badVal.hostRead()) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    #endif
    YAKL_SCOPE( this_concs , this->concs );
    // for (int ilay=1; ilay<=this->nlay; ilay++) {
    //   for (int icol=1; icol<=this->ncol; icol++) {
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      this_concs(icol,ilay,igas) = w(icol,ilay);
    });
  }


  // Get concentration as a 2-D field of columns and levels
  // array is expected to be in devide memory
  void get_vmr(std::string gas, real2d const &array) const {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    if (this->ncol != size(array,1)) { stoprun("ty_gas_concs->get_vmr; gas array is wrong size (ncol)" ); }
    if (this->nlay != size(array,2)) { stoprun("ty_gas_concs->get_vmr; gas array is wrong size (nlay)" ); }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) { stoprun("GasConcs::get_vmr; gas not found" ); }
    // for (int ilay=1; ilay<=size(array,2); ilay++) {
    //   for (int icol=1; icol<=size(array,1); icol++) {
    YAKL_SCOPE( this_concs , this->concs );
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(size(array,2),size(array,1)) , YAKL_LAMBDA (int ilay, int icol) {
      array(icol,ilay) = this_concs(icol,ilay,igas);
    });
  }


  int get_num_gases() const { return yakl::intrinsics::size(gas_name,1); }
  

  string1d get_gas_names() const { return gas_name; }


  // find gas in list; GAS_NOT_IN_LIST if not found
  int find_gas(std::string gas) const {
    if (ngas == 0) { return GAS_NOT_IN_LIST; }
    for (int igas=1; igas<=ngas; igas++) {
      if ( lower_case(gas) == this->gas_name(igas) ) {
        return igas;
      }
    }
    return GAS_NOT_IN_LIST;
  }


  void print_norms() const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

                               std::cout << "ncol      : " << ncol       << "\n";
                               std::cout << "nlay      : " << nlay       << "\n";
                               std::cout << "ngas      : " << ngas       << "\n";
    if (allocated(gas_name)) { std::cout << "gas_name  : " << gas_name   << "\n"; }
    if (allocated(concs   )) { std::cout << "sum(concs): " << sum(concs) << "\n"; }
  }


};


