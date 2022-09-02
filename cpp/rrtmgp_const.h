
#pragma once

#include "YAKL.h"
#include <iostream>
#include <cmath>

template <class T, int rank, int myMem> using FArray = yakl::Array<T,rank,myMem,yakl::styleFortran>;

typedef double real;

YAKL_INLINE real constexpr operator"" _wp( long double x ) {
  return static_cast<real>(x);
}

typedef FArray<real,1,yakl::memDevice> real1d;
typedef FArray<real,2,yakl::memDevice> real2d;
typedef FArray<real,3,yakl::memDevice> real3d;
typedef FArray<real,4,yakl::memDevice> real4d;
typedef FArray<real,5,yakl::memDevice> real5d;
typedef FArray<real,6,yakl::memDevice> real6d;
typedef FArray<real,7,yakl::memDevice> real7d;

typedef FArray<real const,1,yakl::memDevice> realConst1d;
typedef FArray<real const,2,yakl::memDevice> realConst2d;
typedef FArray<real const,3,yakl::memDevice> realConst3d;
typedef FArray<real const,4,yakl::memDevice> realConst4d;
typedef FArray<real const,5,yakl::memDevice> realConst5d;
typedef FArray<real const,6,yakl::memDevice> realConst6d;
typedef FArray<real const,7,yakl::memDevice> realConst7d;

typedef FArray<real,1,yakl::memHost> realHost1d;
typedef FArray<real,2,yakl::memHost> realHost2d;
typedef FArray<real,3,yakl::memHost> realHost3d;
typedef FArray<real,4,yakl::memHost> realHost4d;
typedef FArray<real,5,yakl::memHost> realHost5d;
typedef FArray<real,6,yakl::memHost> realHost6d;
typedef FArray<real,7,yakl::memHost> realHost7d;

typedef FArray<int,1,yakl::memDevice> int1d;
typedef FArray<int,2,yakl::memDevice> int2d;
typedef FArray<int,3,yakl::memDevice> int3d;
typedef FArray<int,4,yakl::memDevice> int4d;
typedef FArray<int,5,yakl::memDevice> int5d;
typedef FArray<int,6,yakl::memDevice> int6d;
typedef FArray<int,7,yakl::memDevice> int7d;

typedef FArray<int const,1,yakl::memDevice> intConst1d;
typedef FArray<int const,2,yakl::memDevice> intConst2d;
typedef FArray<int const,3,yakl::memDevice> intConst3d;
typedef FArray<int const,4,yakl::memDevice> intConst4d;
typedef FArray<int const,5,yakl::memDevice> intConst5d;
typedef FArray<int const,6,yakl::memDevice> intConst6d;
typedef FArray<int const,7,yakl::memDevice> intConst7d;

typedef FArray<int,1,yakl::memHost> intHost1d;
typedef FArray<int,2,yakl::memHost> intHost2d;
typedef FArray<int,3,yakl::memHost> intHost3d;
typedef FArray<int,4,yakl::memHost> intHost4d;
typedef FArray<int,5,yakl::memHost> intHost5d;
typedef FArray<int,6,yakl::memHost> intHost6d;
typedef FArray<int,7,yakl::memHost> intHost7d;

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


