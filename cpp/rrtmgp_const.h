
#pragma once

#ifdef RRTMGP_ENABLE_YAKL
#include "YAKL.h"
#endif

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

#ifdef RRTMGP_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>

using DefaultDevice =
  Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;
using HostDevice =
  Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space>;

template <typename T, typename Device=DefaultDevice>
using FView = Kokkos::View<T, Kokkos::LayoutLeft, Device>;

template <typename T, typename Device=DefaultDevice>
using CView = Kokkos::View<T, Kokkos::LayoutRight, Device>;

template <typename T, typename Device=DefaultDevice>
using FOView = Kokkos::Experimental::OffsetView<T, Kokkos::LayoutLeft, Device>;

template <typename T, typename Device=DefaultDevice>
using COView = Kokkos::Experimental::OffsetView<T, Kokkos::LayoutRight, Device>;

template <int Rank, typename ExecutionSpace=Kokkos::DefaultExecutionSpace>
#ifdef KOKKOS_ENABLE_CUDA
using MDRangeP = Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<Rank, Kokkos::Iterate::Left, Kokkos::Iterate::Right> >;
#else
using MDRangeP = Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<Rank> >;
#endif
#endif

typedef double real;

#ifdef RRTMGP_ENABLE_YAKL
template <class T, int rank, int myMem> using FArray = yakl::Array<T,rank,myMem,yakl::styleFortran>;

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
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
typedef FView<real*>       real1dk;
typedef FView<real**>      real2dk;
typedef FView<real***>     real3dk;
typedef FView<real****>    real4dk;
typedef FView<real*****>   real5dk;
typedef FView<real******>  real6dk;
typedef FView<real*******> real7dk;

typedef FView<const real*>       realConst1dk;
typedef FView<const real**>      realConst2dk;
typedef FView<const real***>     realConst3dk;
typedef FView<const real****>    realConst4dk;
typedef FView<const real*****>   realConst5dk;
typedef FView<const real******>  realConst6dk;
typedef FView<const real*******> realConst7dk;

typedef FView<real*, HostDevice>       realHost1dk;
typedef FView<real**, HostDevice>      realHost2dk;
typedef FView<real***, HostDevice>     realHost3dk;
typedef FView<real****, HostDevice>    realHost4dk;
typedef FView<real*****, HostDevice>   realHost5dk;
typedef FView<real******, HostDevice>  realHost6dk;
typedef FView<real*******, HostDevice> realHost7dk;

typedef FView<int*>       int1dk;
typedef FView<int**>      int2dk;
typedef FView<int***>     int3dk;
typedef FView<int****>    int4dk;
typedef FView<int*****>   int5dk;
typedef FView<int******>  int6dk;
typedef FView<int*******> int7dk;

typedef FView<const int*>       intConst1dk;
typedef FView<const int**>      intConst2dk;
typedef FView<const int***>     intConst3dk;
typedef FView<const int****>    intConst4dk;
typedef FView<const int*****>   intConst5dk;
typedef FView<const int******>  intConst6dk;
typedef FView<const int*******> intConst7dk;

typedef FView<int*, HostDevice>       intHost1dk;
typedef FView<int**, HostDevice>      intHost2dk;
typedef FView<int***, HostDevice>     intHost3dk;
typedef FView<int****, HostDevice>    intHost4dk;
typedef FView<int*****, HostDevice>   intHost5dk;
typedef FView<int******, HostDevice>  intHost6dk;
typedef FView<int*******, HostDevice> intHost7dk;

typedef FView<bool*>       bool1dk;
typedef FView<bool**>      bool2dk;
typedef FView<bool***>     bool3dk;
typedef FView<bool****>    bool4dk;
typedef FView<bool*****>   bool5dk;
typedef FView<bool******>  bool6dk;
typedef FView<bool*******> bool7dk;

typedef FView<bool*, HostDevice>       boolHost1dk;
typedef FView<bool**, HostDevice>      boolHost2dk;
typedef FView<bool***, HostDevice>     boolHost3dk;
typedef FView<bool****, HostDevice>    boolHost4dk;
typedef FView<bool*****, HostDevice>   boolHost5dk;
typedef FView<bool******, HostDevice>  boolHost6dk;
typedef FView<bool*******, HostDevice> boolHost7dk;

typedef FView<char**, HostDevice> charHost2dk;

// this is useful in a couple situations
typedef FOView<real***>             realOff3dk;
typedef FOView<real***, HostDevice> realOffHost3dk;

typedef CView<real*>       real1dkc;
typedef CView<real**>      real2dkc;
typedef CView<real***>     real3dkc;
typedef CView<real****>    real4dkc;
typedef CView<real*****>   real5dkc;
typedef CView<real******>  real6dkc;
typedef CView<real*******> real7dkc;

typedef CView<const real*>       realConst1dkc;
typedef CView<const real**>      realConst2dkc;
typedef CView<const real***>     realConst3dkc;
typedef CView<const real****>    realConst4dkc;
typedef CView<const real*****>   realConst5dkc;
typedef CView<const real******>  realConst6dkc;
typedef CView<const real*******> realConst7dkc;

typedef CView<real*, HostDevice>       realHost1dkc;
typedef CView<real**, HostDevice>      realHost2dkc;
typedef CView<real***, HostDevice>     realHost3dkc;
typedef CView<real****, HostDevice>    realHost4dkc;
typedef CView<real*****, HostDevice>   realHost5dkc;
typedef CView<real******, HostDevice>  realHost6dkc;
typedef CView<real*******, HostDevice> realHost7dkc;

typedef CView<int*>       int1dkc;
typedef CView<int**>      int2dkc;
typedef CView<int***>     int3dkc;
typedef CView<int****>    int4dkc;
typedef CView<int*****>   int5dkc;
typedef CView<int******>  int6dkc;
typedef CView<int*******> int7dkc;

typedef CView<const int*>       intConst1dkc;
typedef CView<const int**>      intConst2dkc;
typedef CView<const int***>     intConst3dkc;
typedef CView<const int****>    intConst4dkc;
typedef CView<const int*****>   intConst5dkc;
typedef CView<const int******>  intConst6dkc;
typedef CView<const int*******> intConst7dkc;

typedef CView<int*, HostDevice>       intHost1dkc;
typedef CView<int**, HostDevice>      intHost2dkc;
typedef CView<int***, HostDevice>     intHost3dkc;
typedef CView<int****, HostDevice>    intHost4dkc;
typedef CView<int*****, HostDevice>   intHost5dkc;
typedef CView<int******, HostDevice>  intHost6dkc;
typedef CView<int*******, HostDevice> intHost7dkc;

typedef CView<bool*>       bool1dkc;
typedef CView<bool**>      bool2dkc;
typedef CView<bool***>     bool3dkc;
typedef CView<bool****>    bool4dkc;
typedef CView<bool*****>   bool5dkc;
typedef CView<bool******>  bool6dkc;
typedef CView<bool*******> bool7dkc;

typedef CView<bool*, HostDevice>       boolHost1dkc;
typedef CView<bool**, HostDevice>      boolHost2dkc;
typedef CView<bool***, HostDevice>     boolHost3dkc;
typedef CView<bool****, HostDevice>    boolHost4dkc;
typedef CView<bool*****, HostDevice>   boolHost5dkc;
typedef CView<bool******, HostDevice>  boolHost6dkc;
typedef CView<bool*******, HostDevice> boolHost7dkc;

typedef CView<char**, HostDevice> charHost2dkc;

typedef COView<real***>             realOff3dkc;
#endif

typedef std::vector<std::string> string1dv;

inline void stoprun( std::string str ) {
  std::cout << "FATAL ERROR:\n";
  std::cout << str << "\n" << std::endl;
  throw str;
}

inline
std::ostream& operator<<(std::ostream& out, const string1dv& names)
{
  for (const auto& name : names) {
    out << name << " ";
  }
  out << std::endl;
  return out;
}
