
#pragma once

#include "YAKL.h"
#include "Array.h"
#include <iostream>
#include <cmath>

template <class T, int rank, int myMem> using FArray = yakl::Array<T,rank,myMem,yakl::styleFortran>;

typedef double real;

YAKL_INLINE real constexpr operator"" _wp( long double x ) {
  return static_cast<real>(x);
}


using yakl::fortran::parallel_for;
using yakl::fortran::Bounds;
using std::max;
using std::min;
using std::abs;
using yakl::intrinsics::mod;
using yakl::intrinsics::merge;
using yakl::intrinsics::size;
using yakl::intrinsics::sum;
using yakl::intrinsics::count;
using yakl::intrinsics::allocated;
using yakl::memHost;
using yakl::memDevice;


typedef FArray<real,1,yakl::memDevice> umgReal1d;
typedef FArray<real,2,yakl::memDevice> umgReal2d;
typedef FArray<real,3,yakl::memDevice> umgReal3d;
typedef FArray<real,4,yakl::memDevice> umgReal4d;
typedef FArray<real,5,yakl::memDevice> umgReal5d;
typedef FArray<real,6,yakl::memDevice> umgReal6d;
typedef FArray<real,7,yakl::memDevice> umgReal7d;

typedef FArray<real,1,yakl::memDevice> real1d;
typedef FArray<real,2,yakl::memDevice> real2d;
typedef FArray<real,3,yakl::memDevice> real3d;
typedef FArray<real,4,yakl::memDevice> real4d;
typedef FArray<real,5,yakl::memDevice> real5d;
typedef FArray<real,6,yakl::memDevice> real6d;
typedef FArray<real,7,yakl::memDevice> real7d;

typedef FArray<real,1,yakl::memHost> realHost1d;
typedef FArray<real,2,yakl::memHost> realHost2d;
typedef FArray<real,3,yakl::memHost> realHost3d;
typedef FArray<real,4,yakl::memHost> realHost4d;
typedef FArray<real,5,yakl::memHost> realHost5d;
typedef FArray<real,6,yakl::memHost> realHost6d;
typedef FArray<real,7,yakl::memHost> realHost7d;

typedef FArray<int,1,yakl::memDevice> umgInt1d;
typedef FArray<int,2,yakl::memDevice> umgInt2d;
typedef FArray<int,3,yakl::memDevice> umgInt3d;
typedef FArray<int,4,yakl::memDevice> umgInt4d;
typedef FArray<int,5,yakl::memDevice> umgInt5d;
typedef FArray<int,6,yakl::memDevice> umgInt6d;
typedef FArray<int,7,yakl::memDevice> umgInt7d;

typedef FArray<int,1,yakl::memDevice> int1d;
typedef FArray<int,2,yakl::memDevice> int2d;
typedef FArray<int,3,yakl::memDevice> int3d;
typedef FArray<int,4,yakl::memDevice> int4d;
typedef FArray<int,5,yakl::memDevice> int5d;
typedef FArray<int,6,yakl::memDevice> int6d;
typedef FArray<int,7,yakl::memDevice> int7d;

typedef FArray<int,1,yakl::memHost> intHost1d;
typedef FArray<int,2,yakl::memHost> intHost2d;
typedef FArray<int,3,yakl::memHost> intHost3d;
typedef FArray<int,4,yakl::memHost> intHost4d;
typedef FArray<int,5,yakl::memHost> intHost5d;
typedef FArray<int,6,yakl::memHost> intHost6d;
typedef FArray<int,7,yakl::memHost> intHost7d;

typedef FArray<bool,1,yakl::memDevice> umgBool1d;
typedef FArray<bool,2,yakl::memDevice> umgBool2d;
typedef FArray<bool,3,yakl::memDevice> umgBool3d;
typedef FArray<bool,4,yakl::memDevice> umgBool4d;
typedef FArray<bool,5,yakl::memDevice> umgBool5d;
typedef FArray<bool,6,yakl::memDevice> umgBool6d;
typedef FArray<bool,7,yakl::memDevice> umgBool7d;

typedef FArray<bool,1,yakl::memDevice> bool1d;
typedef FArray<bool,2,yakl::memDevice> bool2d;
typedef FArray<bool,3,yakl::memDevice> bool3d;
typedef FArray<bool,4,yakl::memDevice> bool4d;
typedef FArray<bool,5,yakl::memDevice> bool5d;
typedef FArray<bool,6,yakl::memDevice> bool6d;
typedef FArray<bool,7,yakl::memDevice> bool7d;

typedef FArray<bool,1,yakl::memHost> boolHost1d;
typedef FArray<bool,2,yakl::memHost> boolHost2d;
typedef FArray<bool,3,yakl::memHost> boolHost3d;
typedef FArray<bool,4,yakl::memHost> boolHost4d;
typedef FArray<bool,5,yakl::memHost> boolHost5d;
typedef FArray<bool,6,yakl::memHost> boolHost6d;
typedef FArray<bool,7,yakl::memHost> boolHost7d;

typedef FArray<char,2,yakl::memHost> charHost2d;

typedef FArray<std::string,1,yakl::memHost> string1d;


inline void stoprun( std::string str ) {
  std::cout << "FATAL ERROR:\n";
  std::cout << str << "\n" << std::endl;
  throw str;
}


template <class T , int rank , int myStyle> inline void zero_array( yakl::Array<T,rank,yakl::memHost,myStyle> &arr ) {
  for (int i=0; i < arr.totElems() ; i++) { arr.myData[i] = 0; }
}


template <class T , int rank , int myStyle> inline void zero_array( yakl::Array<T,rank,yakl::memDevice,myStyle> &arr ) {
  yakl::c::parallel_for( yakl::c::Bounds<1>(arr.totElems()) , YAKL_LAMBDA (int i) {
    arr.myData[i] = 0;
  });
}


template <class T , int rank , int myStyle> inline void memset( yakl::Array<T,rank,yakl::memHost,myStyle> &arr , T val) {
  for (int i=0; i < arr.totElems() ; i++) { arr.myData[i] = val; }
}


template <class T , int rank , int myStyle> inline void memset( yakl::Array<T,rank,yakl::memDevice,myStyle> &arr , T val) {
  yakl::c::parallel_for( yakl::c::Bounds<1>(arr.totElems()) , YAKL_LAMBDA (int i) {
    arr.myData[i] = val;
  });
}



