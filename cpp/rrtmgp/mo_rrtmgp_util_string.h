
#pragma once

#include "const.h"


inline std::string lower_case( std::string in ) {
  std::for_each( in.begin() , in.end() , [] (char & c) { c = ::tolower(c); } );
  return in;
}


inline bool string_in_array(std::string str, string1d const &arr) {
  for (int i=1; i <= size(arr); i++) {
    if ( lower_case(str) == lower_case(arr(i)) ) { return true; }
  }
  return false;
}


inline int string_loc_in_array(std::string str, string1d const &arr) {
  for (int i=1; i <= size(arr); i++) {
    if ( lower_case(str) == lower_case(arr(i)) ) { return i; }
  }
  return -1;
}


inline string1d char2d_to_string1d( charHost2d &in , std::string label="") {
  int nstr  = size(in,2);
  string1d out(label.c_str(),nstr);
  for (int j=1 ; j <= nstr ; j++) {
    out(j) = "";
    for (int i=1 ; i <= size(in,1) ; i++) {
      if ( ! isspace(in(i,j)) ) { out(j) += in(i,j); }
    }
  }
  return out;
}


