
#pragma once

#include "rrtmgp_const.h"

inline std::string lower_case( std::string in ) {
  std::for_each( in.begin() , in.end() , [] (char & c) { c = ::tolower(c); } );
  return in;
}


inline bool string_in_array(std::string str, string1d const &arr) {
  for (int i=1; i <= yakl::intrinsics::size(arr); i++) {
    if ( lower_case(str) == lower_case(arr(i)) ) { return true; }
  }
  return false;
}


inline int string_loc_in_array(std::string str, string1d const &arr) {
  for (int i=1; i <= yakl::intrinsics::size(arr); i++) {
    if ( lower_case(str) == lower_case(arr(i)) ) { return i; }
  }
  return -1;
}

inline bool string_in_array(std::string str, string1dv const &arr) {
  for (const auto& item : arr) {
    if ( lower_case(str) == lower_case(item) ) { return true; }
  }
  return false;
}


inline int string_loc_in_array(std::string str, string1dv const &arr) {
  for (int i=0; i < arr.size(); i++) {
    if ( lower_case(str) == lower_case(arr[i]) ) { return i+1; } // use 1-based idx for now
  }
  return -1;
}


inline string1d char2d_to_string1d( charHost2d &in , std::string label="") {
  int nstr  = yakl::intrinsics::size(in,2);
  string1d out(label.c_str(),nstr);
  for (int j=1 ; j <= nstr ; j++) {
    out(j) = "";
    for (int i=1 ; i <= yakl::intrinsics::size(in,1) ; i++) {
      if ( ! isspace(in(i,j)) ) { out(j) += in(i,j); }
    }
  }
  return out;
}

#ifdef RRTMGP_ENABLE_KOKKOS
inline string1dv char2d_to_string1d( charHost2dk &in , std::string label="") {
  int nstr  = in.extent(1);
  string1dv out(nstr, "");
  for (int j=0 ; j < nstr ; ++j) {
    for (int i=0 ; i < in.extent(0); ++i) {
      if ( ! isspace(in(i,j)) ) { out[j] += in(i,j); }
    }
  }
  return out;
}
#endif
