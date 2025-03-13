#pragma once

#include "rrtmgp_const.h"

#include <stdexcept>
#include <chrono>

// Validate if both enabled?
#ifdef RRTMGP_ENABLE_KOKKOS
#include <netcdf.h>
#ifdef RRTMGP_ENABLE_YAKL
// Both are on, validate
#define COMPUTE_SWITCH(yimpl, kimpl) \
  (kimpl); RRT_REQUIRE((yimpl) == (kimpl), "Bad COMPUTE_SWITCH")
#else
#define COMPUTE_SWITCH(yimpl, kimpl) kimpl
#endif
#else
#define COMPUTE_SWITCH(yimpl, kimpl) yimpl
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
#define GENERIC_INLINE KOKKOS_INLINE_FUNCTION
#define KERNEL_FENCE Kokkos::fence()
#else
#define GENERIC_INLINE YAKL_INLINE
#define KERNEL_FENCE yakl::fence()
#endif

//#define ENABLE_TIMING
// Macro for timing kernels
#ifdef ENABLE_TIMING
#define TIMED_KERNEL(kernel)                                            \
{                                                                       \
  KERNEL_FENCE;                                                         \
  auto start_t = std::chrono::high_resolution_clock::now();             \
  kernel;                                                               \
  KERNEL_FENCE;                                                         \
  auto stop_t = std::chrono::high_resolution_clock::now();              \
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_t - start_t); \
  static double total_s = 0.;                                           \
  total_s += duration.count() / 1000000.0;                              \
  std::cout << "TIMING For func " << __func__ << " file " << __FILE__ << " line " << __LINE__ << " total " << total_s << " s" << std::endl; \
}
#else
#define TIMED_KERNEL(kernel) kernel
#endif

// Macro for timing kernels that does not make a code block. Requires
// a name to disambiguate variable names. Useful when code block defines
// variables that need to persist.
#ifdef ENABLE_TIMING
#define TIMED_INLINE_KERNEL(name, kernel)                               \
  KERNEL_FENCE;                                                         \
  auto start_t##name = std::chrono::high_resolution_clock::now();       \
  kernel;                                                               \
  KERNEL_FENCE;                                                         \
  auto stop_t##name = std::chrono::high_resolution_clock::now();        \
  auto duration##name = std::chrono::duration_cast<std::chrono::microseconds>(stop_t##name - start_t##name); \
  static double total_s##name = 0.;                                     \
  total_s##name += duration##name.count() / 1000000.0;                  \
  std::cout << "TIMING For func " << __func__ << " file " << __FILE__ << " line " << __LINE__ << " total " << total_s##name << " s" << std::endl
#else
#define TIMED_INLINE_KERNEL(name, kernel) kernel
#endif


/**
 * Helper functions for the conversion to Kokkos
 */

namespace conv {

// Print an entire view + sentinel
template <typename KView>
void printk(const std::string& name, const KView& view)
{
  for (size_t i = 0; i < view.size(); ++i) {
    std::cout << "JGFK " << name << "(" << i << ") = " << view.data()[i] << std::endl;
  }
}

// Print an entire yakl array + sentinel
template <typename YArray>
void printy(const std::string& name, const YArray& array)
{
  for (size_t i = 0; i < array.totElems(); ++i) {
    std::cout << "JGFY " << name << "(" << i << ") = " << array.data()[i] << std::endl;
  }
}

// Copied from YAKL
template <class T1, class T2,
          typename std::enable_if<std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,bool>::type=false>
GENERIC_INLINE decltype(T1()+T2()) merge(T1 const t, T2 const f, bool cond) noexcept { return cond ? t : f; }

// A meta function that will return true if T is a Kokkos::View
template <typename T>
struct is_view
{
#ifdef RRTMGP_ENABLE_KOKKOS
  static constexpr bool value = Kokkos::is_view<T>::value;
#else
  static constexpr bool value = false;
#endif
};

// Convenient way of using is_view meta function

template <class T>
inline constexpr bool is_view_v = is_view<T>::value;

// Copied with minor mods from YAKL intrinsics. Gets epsilon for a primitive
// or View.
template <class T,
          typename std::enable_if<!is_view_v<T>>::type* = nullptr>
GENERIC_INLINE
T constexpr epsilon(T) { return std::numeric_limits<T>::epsilon(); }

template <class View,
          typename std::enable_if<is_view_v<View>>::type* = nullptr>
GENERIC_INLINE
typename View::non_const_value_type constexpr epsilon(const View& arr)
{ return std::numeric_limits<typename View::non_const_value_type>::epsilon(); }

//
// These are for debugging. They print values of kviews/yarrays preceeded by a
// sentinel string that can be grepped for. This is how the compare_yk.sh script
// can be used to debug differences between YAKL and Kokkos. You put the pNd calls
// in the equivalent place for both the YAKL and Kokkos version of the code.
//
template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p1d(const KView& view, const std::string& name, int idx)
{ std::cout << "JGFK " << name << "(" << idx << ") = " << view(idx) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p1d(const YArray& array, const std::string& name, int idx)
{
  const int adjust_val = std::is_same_v<typename YArray::non_const_value_type, int> ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx-1 << ") = " << array(idx) - adjust_val << std::endl; }

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p2d(const KView& view, const std::string& name, int idx1, int idx2)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ") = " << view(idx1, idx2) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p2d(const YArray& array, const std::string& name, int idx1, int idx2)
{
  const int adjust_val = std::is_same_v<typename YArray::non_const_value_type, int> ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx1-1 << ", " << idx2-1 << ") = " << array(idx1, idx2) - adjust_val << std::endl; }

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p3d(const KView& view, const std::string& name, int idx1, int idx2, int idx3)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ", " << idx3 << ") = " << view(idx1, idx2, idx3) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p3d(const YArray& array, const std::string& name, int idx1, int idx2, int idx3)
{
  const int adjust_val = std::is_same_v<typename YArray::non_const_value_type, int> ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx1-1 << ", " << idx2-1 << ", " << idx3-1 << ") = " << array(idx1, idx2, idx3) - adjust_val << std::endl;
}

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p4d(const KView& view, const std::string& name, int idx1, int idx2, int idx3, int idx4)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ", " << idx3 << ", " << idx4 << ") = " << view(idx1, idx2, idx3, idx4) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p4d(const YArray& array, const std::string& name, int idx1, int idx2, int idx3, int idx4)
{
  const int adjust_val = std::is_same_v<typename YArray::non_const_value_type, int> ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx1-1 << ", " << idx2-1 << ", " << idx3-1 << ", " << idx4-1 << ") = " << array(idx1, idx2, idx3, idx4) - adjust_val << std::endl;
}

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p5d(const KView& view, const std::string& name, int idx1, int idx2, int idx3, int idx4, int idx5)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ", " << idx3 << ", " << idx4 << ", " << idx5 << ") = " << view(idx1, idx2, idx3, idx4, idx5) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p5d(const YArray& array, const std::string& name, int idx1, int idx2, int idx3, int idx4, int idx5)
{
  const int adjust_val = std::is_same_v<typename YArray::non_const_value_type, int> ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx1-1 << ", " << idx2-1 << ", " << idx3-1 << ", " << idx4-1 << ", " << idx5-1 << ") = " << array(idx1, idx2, idx3, idx4, idx5) - adjust_val << std::endl;
}

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p6d(const KView& view, const std::string& name, int idx1, int idx2, int idx3, int idx4, int idx5, int idx6)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ", " << idx3 << ", " << idx4 << ", " << idx5 << ", " << idx6 << ") = " << view(idx1, idx2, idx3, idx4, idx5, idx6) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p6d(const YArray& array, const std::string& name, int idx1, int idx2, int idx3, int idx4, int idx5, int idx6)
{
  const int adjust_val = std::is_same_v<typename YArray::non_const_value_type, int> ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx1-1 << ", " << idx2-1 << ", " << idx3-1 << ", " << idx4-1 << ", " << idx5-1 << ", " << idx6-1 << ") = " << array(idx1, idx2, idx3, idx4, idx5, idx6) - adjust_val << std::endl;
}

// Copied from EKAT
#define IMPL_THROW_RRT(condition, msg, exception_type)    \
do {                                                      \
  if ( ! (condition) ) {                                  \
    std::stringstream _ss_;                               \
    _ss_ << "\n FAIL:\n" << #condition  << "\n";          \
    _ss_ << "\n" << msg;                                  \
    throw exception_type(_ss_.str());                     \
  }                                                       \
} while(0)

// Just a convenient require macro for checking things. Note this is
// not an assert macro, so it's always on.
#define RRT_REQUIRE(condition, msg) IMPL_THROW_RRT(condition, msg, std::runtime_error)

#ifdef RRTMGP_ENABLE_KOKKOS

// Macros for validating kokkos using yakl. These all do nothing if
// YAKL is not enabled.
#ifdef RRTMGP_ENABLE_YAKL
#define VALIDATE_KOKKOS(yobj, kobj) kobj.validate_kokkos(yobj)
#else
#define VALIDATE_KOKKOS(yobj, kobj) (void)0
#endif

#ifdef RRTMGP_ENABLE_YAKL
#define COMPARE_ALL_WRAP(yobjs, kobjs) conv::compare_all_yakl_to_kokkos(yobjs, kobjs)
#else
#define COMPARE_ALL_WRAP(yobjs, kobjs) (void)0
#endif

#ifdef RRTMGP_ENABLE_YAKL
#define COMPARE_WRAP(yobj, kobj) conv::compare_yakl_to_kokkos(yobj, kobj)
#else
#define COMPARE_WRAP(yobj, kobj) (void)0
#endif

// Copied from EKAT. Get Kokkos view template type
template<typename T, int N>
struct DataND {
  using type = typename DataND<T,N-1>::type*;
};
template<typename T>
struct DataND<T,0> {
  using type = T;
};

// Copied from EKAT: Create a kokkos layout from a container of dims
template <typename Layout, typename Container>
Layout get_layout(const Container& dims)
{
  Layout result;
  int dimitr = 0;
  for (const auto& dim : dims) {
    result.dimension[dimitr++] = dim;
  }
  return result;
}

// Copied from EKAT
template <typename View>
struct MemoryTraitsMask {
  enum : unsigned int {
    value = ((View::traits::memory_traits::is_random_access ? Kokkos::RandomAccess : 0) |
             (View::traits::memory_traits::is_atomic ? Kokkos::Atomic : 0) |
             (View::traits::memory_traits::is_restrict ? Kokkos::Restrict : 0) |
             (View::traits::memory_traits::is_aligned ? Kokkos::Aligned : 0) |
             (View::traits::memory_traits::is_unmanaged ? Kokkos::Unmanaged : 0))
      };
};

// Copied from EKAT
template <typename View>
using Unmanaged =
  // Provide a full View type specification, augmented with Unmanaged.
  Kokkos::View<typename View::traits::scalar_array_type,
               typename View::traits::array_layout,
               typename View::traits::device_type,
               Kokkos::MemoryTraits<
                 // All the current values...
                 MemoryTraitsMask<View>::value |
                 // ... |ed with the one we want, whether or not it's
                 // already there.
                 Kokkos::Unmanaged> >;

// approx eq. Compares values with a tolerance if reals.
template <typename T,
          typename std::enable_if<std::is_integral_v<T>>::type* = nullptr>
inline
bool approx_eq(const T lhs, const T rhs)
{
  return lhs == rhs;
}

template <typename T,
          typename std::enable_if<std::is_floating_point_v<T>>::type* = nullptr>
inline
bool approx_eq(const T lhs, const T rhs)
{
  constexpr T tol = 1e-12;
  if (std::isnan(lhs) && std::isnan(rhs)) {
    return true;
  }
  return std::abs(lhs - rhs) < tol;
}

// get_mdrp, returns multidimensional range policy. Start is assumed to be {0}*N
template <int N> struct DefaultTile;

#define OMEGA_TILE_LENGTH 64

template <> struct DefaultTile<1> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<2> {
   static constexpr int value[] = {1, OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<3> {
   static constexpr int value[] = {1, 1, OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<4> {
   static constexpr int value[] = {1, 1, 1, OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<5> {
   static constexpr int value[] = {1, 1, 1, 1, OMEGA_TILE_LENGTH};
};

// template <> struct DefaultTile<1> {
//    static constexpr int value[] = {OMEGA_TILE_LENGTH};
// };

// template <> struct DefaultTile<2> {
//    static constexpr int value[] = {OMEGA_TILE_LENGTH, 1};
// };

// template <> struct DefaultTile<3> {
//    static constexpr int value[] = {OMEGA_TILE_LENGTH, 1, 1};
// };

// template <> struct DefaultTile<4> {
//    static constexpr int value[] = {OMEGA_TILE_LENGTH, 1, 1, 1};
// };

// template <> struct DefaultTile<5> {
//    static constexpr int value[] = {OMEGA_TILE_LENGTH, 1, 1, 1, 1};
// };

template <typename LayoutT, typename DeviceT=DefaultDevice>
struct MDRP
{
  // By default, follow the Layout's fast index
  static constexpr Kokkos::Iterate LeftI = std::is_same_v<LayoutT, Kokkos::LayoutRight>
    ? Kokkos::Iterate::Left
    : Kokkos::Iterate::Right;
  static constexpr Kokkos::Iterate RightI = std::is_same_v<LayoutT, Kokkos::LayoutRight>
    ? Kokkos::Iterate::Right
    : Kokkos::Iterate::Left;

  using exe_space_t = typename DeviceT::execution_space;

  template <int Rank>
  using MDRP_t = Kokkos::MDRangePolicy<exe_space_t, Kokkos::Rank<Rank, LeftI, RightI> >;

  // Force a right to left (right fast index) loop order
  template <int Rank>
  using MDRPRL_t = Kokkos::MDRangePolicy<exe_space_t, Kokkos::Rank<Rank, Kokkos::Iterate::Left, Kokkos::Iterate::Right> >;

  // Force a left to right (left fast index) loop order
  template <int Rank>
  using MDRPLR_t = Kokkos::MDRangePolicy<exe_space_t, Kokkos::Rank<Rank, Kokkos::Iterate::Right, Kokkos::Iterate::Left> >;

  template <int N, typename IntT>
  static inline
  MDRP_t<N> get(const IntT (&upper_bounds)[N])
  {
    assert(N > 1);
    const IntT lower_bounds[N] = {0};
    return MDRP_t<N>(lower_bounds, upper_bounds); //, DefaultTile<N>::value);
  }

  template <int N, typename IntT>
  static inline
  MDRPLR_t<N> getlr(const IntT (&upper_bounds)[N])
  {
    assert(N > 1);
    const IntT lower_bounds[N] = {0};
    return MDRPLR_t<N>(lower_bounds, upper_bounds); //, DefaultTile<N>::value);
  }


  template <int N, typename IntT>
  static inline
  MDRPRL_t<N> getrl(const IntT (&upper_bounds)[N])
  {
    assert(N > 1);
    const IntT lower_bounds[N] = {0};
    return MDRPRL_t<N>(lower_bounds, upper_bounds); //, DefaultTile<N>::value);
  }

  template <int N, typename IntT>
  static inline
  MDRP_t<N> get(const IntT (&lower_bounds)[N], const IntT (&upper_bounds)[N])
  {
    assert(N > 1);
    return MDRP_t<N>(lower_bounds, upper_bounds); //, DefaultTile<N>::value);
  }
};

KOKKOS_INLINE_FUNCTION
void unflatten_idx_left(const int idx, const Kokkos::Array<int, 2>& dims, int& i, int& j)
{
  i = idx % dims[0];
  j = idx / dims[0];
}

KOKKOS_INLINE_FUNCTION
void unflatten_idx_left(const int idx, const Kokkos::Array<int, 3>& dims, int& i, int& j, int& k)
{
  i = idx % dims[0];
  j = (idx / dims[0]) % dims[1];
  k = idx / (dims[0] * dims[1]);
}

KOKKOS_INLINE_FUNCTION
void unflatten_idx_left(const int idx, const Kokkos::Array<int, 4>& dims, int& i, int& j, int& k, int& l)
{
  i = idx % dims[0];
  j = (idx / dims[0]) % dims[1];
  k = (idx / (dims[0]*dims[1])) % dims[2];
  l = idx / (dims[0]*dims[1]*dims[2]);
}

KOKKOS_INLINE_FUNCTION
void unflatten_idx_right(const int idx, const Kokkos::Array<int, 2>& dims, int& i, int& j)
{
  i = idx / dims[1];
  j = idx % dims[1];
}

KOKKOS_INLINE_FUNCTION
void unflatten_idx_right(const int idx, const Kokkos::Array<int, 3>& dims, int& i, int& j, int& k)
{
  i = idx / (dims[2] * dims[1]);
  j = (idx / dims[2]) % dims[1];
  k =  idx % dims[2];
}

KOKKOS_INLINE_FUNCTION
void unflatten_idx_right(const int idx, const Kokkos::Array<int, 4>& dims, int& i, int& j, int& k, int& l)
{
  i = idx / (dims[3]*dims[2]*dims[1]);
  j = (idx / (dims[3]*dims[2])) % dims[1];
  k = (idx / dims[3]) % dims[2];
  l = idx % dims[3];
}

template <typename LayoutT>
KOKKOS_INLINE_FUNCTION
void unflatten_idx(const int idx, const Kokkos::Array<int, 2>& dims, int& i, int& j)
{
  if constexpr (std::is_same_v<LayoutT, Kokkos::LayoutLeft>) {
    unflatten_idx_left(idx, dims, i, j);
  }
  else {
    unflatten_idx_right(idx, dims, i, j);
  }
}

template <typename LayoutT>
KOKKOS_INLINE_FUNCTION
void unflatten_idx(const int idx, const Kokkos::Array<int, 3>& dims, int& i, int& j, int& k)
{
  if constexpr (std::is_same_v<LayoutT, Kokkos::LayoutLeft>) {
    unflatten_idx_left(idx, dims, i, j, k);
  }
  else {
    unflatten_idx_right(idx, dims, i, j, k);
  }
}

template <typename LayoutT>
KOKKOS_INLINE_FUNCTION
void unflatten_idx(const int idx, const Kokkos::Array<int, 4>& dims, int& i, int& j, int& k, int& l)
{
  if constexpr (std::is_same_v<LayoutT, Kokkos::LayoutLeft>) {
    unflatten_idx_left(idx, dims, i, j, k, l);
  }
  else {
    unflatten_idx_right(idx, dims, i, j, k, l);
  }
}

// The FLATTEN* macros expect LayoutT to be defined.

#define FLATTEN_MD_KERNEL2(n1, n2, i1, i2, kernel)     \
  {                                                     \
    Kokkos::Array<int, 2> dims_fmk_internal = {n1, n2};         \
    const int dims_fmk_internal_tot = (n1)*(n2);                        \
    Kokkos::parallel_for(dims_fmk_internal_tot, KOKKOS_LAMBDA (int idx_fmk_internal) { \
      int i1, i2;                                                     \
      conv::unflatten_idx<LayoutT>(idx_fmk_internal, dims_fmk_internal, i1, i2); \
      kernel;                                                           \
    });                                                               \
  }

#define FLATTEN_MD_KERNEL3(n1, n2, n3, i1, i2, i3, kernel)       \
  {                                                              \
    Kokkos::Array<int, 3> dims_fmk_internal = {n1, n2, n3};      \
    const int dims_fmk_internal_tot = (n1)*(n2)*(n3);                   \
    Kokkos::parallel_for(dims_fmk_internal_tot, KOKKOS_LAMBDA (int idx_fmk_internal) { \
      int i1, i2, i3;                                                 \
      conv::unflatten_idx<LayoutT>(idx_fmk_internal, dims_fmk_internal, i1, i2, i3); \
      kernel;                                                           \
    });                                                               \
  }

#define FLATTEN_MD_KERNEL4(n1, n2, n3, n4, i1, i2, i3, i4, kernel)       \
  {                                                                     \
    Kokkos::Array<int, 4> dims_fmk_internal = {n1, n2, n3, n4};         \
    const int dims_fmk_internal_tot = (n1)*(n2)*(n3)*(n4);              \
    Kokkos::parallel_for(dims_fmk_internal_tot, KOKKOS_LAMBDA (int idx_fmk_internal) { \
      int i1, i2, i3, i4;                                             \
      conv::unflatten_idx<LayoutT>(idx_fmk_internal, dims_fmk_internal, i1, i2, i3, i4); \
      kernel;                                                           \
    });                                                               \
  }


#ifdef RRTMGP_ENABLE_YAKL
// Compare a yakl array to a kokkos view, checking they are functionally
// identical (same rank, dims, and values).
template <typename YArray, typename KView>
void compare_yakl_to_kokkos(const YArray& yarray, const KView& kview, bool index_data=false)
{
  using yakl::intrinsics::size;
  using LeftHostView = Kokkos::View<typename KView::non_const_data_type, Kokkos::LayoutLeft, HostDevice>;

  constexpr auto krank = KView::rank;
  const auto yrank = yarray.get_rank();

  RRT_REQUIRE(krank == yrank, "Rank mismatch for: " << kview.label());

  Kokkos::LayoutLeft llayout;
  for (auto r = 0; r < krank; ++r) {
    llayout.dimension[r] = kview.layout().dimension[r];
  }
  LeftHostView hkview("read_data", llayout);
  Kokkos::deep_copy(hkview, kview);

  RRT_REQUIRE(kview.is_allocated() == yakl::intrinsics::allocated(yarray),
              "Allocation status mismatch for: " << kview.label());
  if (!kview.is_allocated()) {
    // we're done
    return;
  }

  auto hyarray = yarray.createHostCopy();

  for (auto r = 0; r < krank; ++r) {
    RRT_REQUIRE(kview.extent(r) == size(yarray,r+1), "Dim mismatch for: " << kview.label() << ", rank: " << r << ", " << kview.extent(r) << " != " <<  size(yarray,r+1));
  }

  auto total_size = kview.size();
  for (auto i = 0; i < total_size; ++i) {
    const auto kdata = hkview.data()[i];
    const auto ydata = hyarray.data()[i];
    if (index_data) {
      if (kdata < 0 && ydata <= 0) {
        // pass
      }
      else {
        RRT_REQUIRE((kdata + (index_data ? 1 : 0)) == ydata, "Data mismatch for: " << kview.label() << ", i: " << i << ", " << kdata << " != " << ydata);
      }
    }
    else {
      RRT_REQUIRE(approx_eq(kdata, ydata), "Data mismatch for: " << kview.label() << ", i: " << i << ", " << kdata << " != " << ydata);
    }

  }
}

inline
void compare_yakl_to_kokkos_str(const string1d& yarray, const string1dv& kstrs)
{
  using yakl::intrinsics::size;

  RRT_REQUIRE(size(yarray, 1) == kstrs.size(), "Dim mistmatch for: " << yarray.label());
  for (auto i = 0; i < kstrs.size(); ++i) {
    RRT_REQUIRE(yarray(i+1) == kstrs[i], "Data mismatch for: " << yarray.label());
  }
}

template <typename YArray, typename KView>
void compare_all_yakl_to_kokkos(const std::vector<YArray>& yarrays, const std::vector<KView>& kviews)
{
  RRT_REQUIRE(yarrays.size() == kviews.size(), "Mismatched vector lengths");
  for (size_t i = 0; i < yarrays.size(); ++i) {
    compare_yakl_to_kokkos(yarrays[i], kviews[i]);
  }
}

template <typename KView>
struct ToYakl
{
  using scalar_t = typename KView::value_type;
  static constexpr auto yakl_mem = std::is_same_v<typename KView::device_type, HostDevice> ? yakl::memHost : yakl::memDevice;
  using type = FArray<scalar_t, KView::rank, yakl_mem>;
};

// Allocate and populate a yakl array from a kokkos view
template <typename KView>
typename ToYakl<KView>::type to_yakl(const KView& view)
{
  using yarray_t    = typename ToYakl<KView>::type;
  using exe_space_t = typename KView::execution_space;

  std::vector<int> dims(KView::rank); // List of dimensions for this variable
  for (auto i = 0; i < KView::rank; ++i) {
    dims[i] = view.extent(i);
  }
  yarray_t rv(view.name(), dims);
  Kokkos::parallel_for( Kokkos::RangePolicy<exe_space_t>(0, view.size()),
                        KOKKOS_LAMBDA(size_t i) {
    rv.data()[i] = view.data()[i];
  });
  return rv;
}
#endif

// A < functor
template <typename T>
struct LTFunc
{
  T m_val;
  LTFunc(T val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(const T& val) const noexcept
  {
    return val < m_val;
  }
};

// Return true if any item in view returns true when passed to functor
template <typename KView, typename Functor>
bool any(const KView& view, const Functor& functor)
{
  using exe_space_t = typename KView::execution_space;

  bool rv = false;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exe_space_t>(0, view.size()),
    KOKKOS_LAMBDA(size_t i, bool& val) {
      val |= functor(view.data()[i]);
    }, Kokkos::LOr<bool>(rv));

  return rv;
}

// Get max val of View
template <typename KView>
typename KView::non_const_value_type maxval(const KView& view)
{
  using scalar_t    = typename KView::non_const_value_type;
  using exe_space_t = typename KView::execution_space;

  scalar_t rv;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exe_space_t>(0, view.size()),
    KOKKOS_LAMBDA(size_t i, scalar_t& lmax) {
      const scalar_t val = view.data()[i];
      if (val > lmax) lmax = val;
    }, Kokkos::Max<scalar_t>(rv));
  return rv;
}

// Get min val of View
template <typename KView>
typename KView::non_const_value_type minval(const KView& view)
{
  using scalar_t    = typename KView::non_const_value_type;
  using exe_space_t = typename KView::execution_space;

  scalar_t rv;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exe_space_t>(0, view.size()),
    KOKKOS_LAMBDA(size_t i, scalar_t& lmin) {
      const scalar_t val = view.data()[i];
      if (val < lmin) lmin = val;
    }, Kokkos::Min<scalar_t>(rv));
  return rv;
}

// Get sum of view
template <typename KView>
std::conditional_t<std::is_same_v<typename KView::non_const_value_type, bool>, int, typename KView::non_const_value_type>
sum(const KView& view)
{
  using scalar_t    = typename KView::non_const_value_type;
  using exe_space_t = typename KView::execution_space;
  using sum_t       = std::conditional_t<std::is_same_v<scalar_t, bool>, int, scalar_t>;

  // If comparing sums of ints against f90, the values will need to be adjusted if
  // the sums are going to match. This would only be done during debugging, so disable
  // this for now.
  // auto adjust = std::is_same_v<scalar_t, int> ? 1 : 0;

  sum_t rv;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exe_space_t>(0, view.size()),
    KOKKOS_LAMBDA(size_t i, sum_t& lsum) {
      lsum += (view.data()[i] /*+ adjust*/);
    }, Kokkos::Sum<sum_t>(rv));
  return rv;
}

// MemPool singleton. Stack allocation pattern only.
template <typename RealT, typename LayoutT, typename DeviceT>
struct MemPoolSingleton
{
 public:
  using memview_t = Kokkos::View<RealT*, LayoutT, DeviceT>;

  template <typename T>
  using view_t = Kokkos::View<T, LayoutT, DeviceT>;

  static inline memview_t  s_mem;
  static inline int64_t s_curr_used;
  static inline int64_t s_high_water;

  static void init(const int64_t capacity)
  {
    static bool is_init = false;
    RRT_REQUIRE(!is_init, "Multiple MemPoolSingleton inits");
    s_mem = memview_t("s_mem", capacity);
    s_curr_used = 0;
    s_high_water = 0;
  }

  template <typename T>
  static inline
  int64_t get_num_reals(const int64_t num) noexcept
  {
    assert(sizeof(T) <= sizeof(RealT));
    static constexpr int64_t CACHE_LINE_SIZE = 64;
    static constexpr int64_t reals_per_cache_line = CACHE_LINE_SIZE / sizeof(RealT);
    const int64_t num_reals = ((num * sizeof(T) + (sizeof(RealT) - 1)) / sizeof(RealT));
    // + reals_per_cache_line; // pad. This didn't seem to help at all
    const int64_t num_reals_cache_aligned = ((num_reals + reals_per_cache_line - 1) / reals_per_cache_line) * reals_per_cache_line;
    return num_reals_cache_aligned;
  }

  /**
   * Allocate and return raw memory. This is useful for when you
   * want to batch several view allocations at once.
   */
  template <typename T>
  static inline
  T* alloc_raw(const int64_t num) noexcept
  {
    const int64_t num_reals = get_num_reals<T>(num);
    T* rv = reinterpret_cast<T*>(s_mem.data() + s_curr_used);
    s_curr_used += num_reals;
    assert(s_curr_used <= s_mem.size());
    if (s_curr_used > s_high_water) {
      s_high_water = s_curr_used;
    }
    return rv;
  }

  template <typename T>
  static inline
  auto alloc(const int64_t dim1) noexcept
  {
    using uview_t = Unmanaged<view_t<T*>>;
    return uview_t(alloc_raw<T>(dim1), dim1);
  }

  template <typename T>
  static inline
  auto alloc(const int64_t dim1, const int64_t dim2) noexcept
  {
    using uview_t = Unmanaged<view_t<T**>>;
    return uview_t(alloc_raw<T>(dim1*dim2), dim1, dim2);
  }

  template <typename T>
  static inline
  auto alloc(const int64_t dim1, const int64_t dim2, const int64_t dim3) noexcept
  {
    using uview_t = Unmanaged<view_t<T***>>;
    return uview_t(alloc_raw<T>(dim1*dim2*dim3), dim1, dim2, dim3);
  }

  template <typename T>
  static inline
  auto alloc(const int64_t dim1, const int64_t dim2, const int64_t dim3, const int64_t dim4) noexcept
  {
    using uview_t = Unmanaged<view_t<T****>>;
    return uview_t(alloc_raw<T>(dim1*dim2*dim3*dim4), dim1, dim2, dim3, dim4);
  }

  template <typename T>
  static inline
  auto alloc(const int64_t dim1, const int64_t dim2, const int64_t dim3, const int64_t dim4, const int dim5) noexcept
  {
    using uview_t = Unmanaged<view_t<T*****>>;
    return uview_t(alloc_raw<T>(dim1*dim2*dim3*dim4*dim5), dim1, dim2, dim3, dim4, dim5);
  }

  template <typename T>
  static inline
  auto alloc(const int64_t dim1, const int64_t dim2, const int64_t dim3, const int64_t dim4, const int dim5, const int dim6) noexcept
  {
    using uview_t = Unmanaged<view_t<T******>>;
    return uview_t(alloc_raw<T>(dim1*dim2*dim3*dim4*dim5*dim6), dim1, dim2, dim3, dim4, dim5, dim6);
  }

  template <typename T>
  static inline
  auto alloc_and_init(const int64_t dim1) noexcept
  {
    using uview_t = Unmanaged<view_t<T*>>;
    uview_t rv(alloc_raw<T>(dim1), dim1);
    Kokkos::deep_copy(rv, 0);
    return rv;
  }

  template <typename T>
  static inline
  auto alloc_and_init(const int64_t dim1, const int64_t dim2) noexcept
  {
    using uview_t = Unmanaged<view_t<T**>>;
    uview_t rv(alloc_raw<T>(dim1*dim2), dim1, dim2);
    Kokkos::deep_copy(rv, 0);
    return rv;
  }

  template <typename T>
  static inline
  auto alloc_and_init(const int64_t dim1, const int64_t dim2, const int64_t dim3) noexcept
  {
    using uview_t = Unmanaged<view_t<T***>>;
    uview_t rv(alloc_raw<T>(dim1*dim2*dim3), dim1, dim2, dim3);
    Kokkos::deep_copy(rv, 0);
    return rv;
  }

  template <typename T>
  static inline
  auto alloc_and_init(const int64_t dim1, const int64_t dim2, const int64_t dim3, const int64_t dim4) noexcept
  {
    using uview_t = Unmanaged<view_t<T****>>;
    uview_t rv(alloc_raw<T>(dim1*dim2*dim3*dim4), dim1, dim2, dim3, dim4);
    Kokkos::deep_copy(rv, 0);
    return rv;
  }

  template <typename T>
  static inline
  auto alloc_and_init(const int64_t dim1, const int64_t dim2, const int64_t dim3, const int64_t dim4, const int dim5) noexcept
  {
    using uview_t = Unmanaged<view_t<T*****>>;
    uview_t rv(alloc_raw<T>(dim1*dim2*dim3*dim4*dim5), dim1, dim2, dim3, dim4, dim5);
    Kokkos::deep_copy(rv, 0);
    return rv;
  }

  template <typename T>
  static inline
  auto alloc_and_init(const int64_t dim1, const int64_t dim2, const int64_t dim3, const int64_t dim4, const int dim5, const int dim6) noexcept
  {
    using uview_t = Unmanaged<view_t<T******>>;
    uview_t rv(alloc_raw<T>(dim1*dim2*dim3*dim4*dim5*dim6), dim1, dim2, dim3, dim4, dim5, dim6);
    Kokkos::deep_copy(rv, 0);
    return rv;
  }

  template <typename T>
  static inline
  void dealloc(const T*, const int64_t num) noexcept
  {
    const int64_t num_reals = get_num_reals<T>(num);
    s_curr_used -= num_reals;
    assert(s_curr_used >= 0);
  }

  template <typename View>
  static inline
  void dealloc(const View& view) noexcept
  {
    dealloc(view.data(), view.size());
  }

  static inline
  void finalize(bool verbose=true)
  {
    if (verbose) print_state();
    assert(s_curr_used == 0); // !=0 indicates we may have forgotten a dealloc
    s_mem = memview_t();
  }

  static inline
  void print_state()
  {
    std::cout << "rrtmgp_conversion MemPoolSingleton used " << s_curr_used << " out of " << s_mem.size() << "; high_water was " << s_high_water << std::endl;
  }
};

//Error reporting routine for the PNetCDF I/O
/** @private */
inline void ncwrap( int ierr , int line ) {
  if (ierr != NC_NOERR) {
    printf("NetCDF Error at line: %d\n", line);
    printf("%s\n",nc_strerror(ierr));
    throw std::runtime_error("");
  }
}

/** @brief Tells NetCDF the opened file should be opened for reading only. */
int constexpr NETCDF_MODE_READ    = NC_NOWRITE;
/** @brief Tells NetCDF the opened file should be opened for reading and writing. */
int constexpr NETCDF_MODE_WRITE   = NC_WRITE;
/** @brief Tells NetCDF the created file should overwite a file of the same name. */
int constexpr NETCDF_MODE_REPLACE = NC_CLOBBER;
/** @brief Tells NetCDF the created file should not overwite a file of the same name. */
int constexpr NETCDF_MODE_NEW     = NC_NOCLOBBER;

// Copied with minor modifications from YAKL_netcdf.h. This is a useful utility that
// should go somewhere.. maybe EKAT?
class SimpleNetCDF {
public:
  /** @private */
  class NcDim {
  public:
    std::string name;
    size_t      len;
    int         id;
    bool        is_unlim;
    NcDim() {
      name = "";
      id = -999;
      len = 0;
      is_unlim = false;
    }
    ~NcDim() {}
    NcDim(std::string name, size_t len, int id, bool is_unlim) {
      this->name     = name;
      this->len      = len;
      this->id       = id;
      this->is_unlim = is_unlim;
    }
    NcDim(NcDim &&in) {
      this->name     = in.name;
      this->len      = in.len;
      this->id       = in.id;
      this->is_unlim = in.is_unlim;
    }
    NcDim(NcDim const &in) {
      this->name     = in.name;
      this->len      = in.len;
      this->id       = in.id;
      this->is_unlim = in.is_unlim;
    }
    NcDim &operator=(NcDim &&in) {
      this->name     = in.name;
      this->len      = in.len;
      this->id       = in.id;
      this->is_unlim = in.is_unlim;
      return *this;
    }
    NcDim &operator=(NcDim const &in) {
      this->name     = in.name;
      this->len      = in.len;
      this->id       = in.id;
      this->is_unlim = in.is_unlim;
      return *this;
    }
    std::string getName()                    const { return name; }
    size_t      getSize()                    const { return len; }
    int         getId()                      const { return id; }
    bool        isNull()                     const { return id == -999; }
    bool        operator==(NcDim const &rhs) const { return this->name == rhs.name && !isNull(); }
    bool        operator!=(NcDim const &rhs) const { return this->name != rhs.name || isNull(); }
    bool        isUnlimited()                const { return is_unlim; }
  };


  /** @private */
  class NcVar {
  public:
    int                ncid;
    std::string        name;
    std::vector<NcDim> dims;
    int                id;
    int                type;
    NcVar() {
      ncid = -999;
      name = "";
      dims = std::vector<NcDim>(0);
      id   = -999;
      type = -999;
    }
    ~NcVar() {}
    NcVar(int ncid , std::string name, std::vector<NcDim> dims, int id, int type) {
      this->ncid = ncid;
      this->name = name;
      this->dims = dims;
      this->id   = id;
      this->type = type;
    }
    NcVar(NcVar &&in) {
      this->ncid = in.ncid;
      this->name = in.name;
      this->dims = in.dims;
      this->id   = in.id;
      this->type = in.type;
    }
    NcVar(NcVar const &in) {
      this->ncid = in.ncid;
      this->name = in.name;
      this->dims = in.dims;
      this->id   = in.id;
      this->type = in.type;
    }
    NcVar &operator=(NcVar &&in) {
      this->ncid = in.ncid;
      this->name = in.name;
      this->dims = in.dims;
      this->id   = in.id;
      this->type = in.type;
      return *this;
    }
    NcVar &operator=(NcVar const &in) {
      this->ncid = in.ncid;
      this->name = in.name;
      this->dims = in.dims;
      this->id   = in.id;
      this->type = in.type;
      return *this;
    }
    std::string        getName()                    const { return name; }
    std::vector<NcDim> getDims()                    const { return dims; }
    int                getDimCount()                const { return dims.size(); }
    int                getId()                      const { return id; }
    int                getType()                    const { return type; }
    bool               isNull ()                    const { return id == -999; }
    bool               operator==(NcDim const &rhs) const { return this->name == rhs.name && !isNull(); }
    bool               operator!=(NcDim const &rhs) const { return this->name != rhs.name || isNull(); }
    NcDim getDim(int i) const {
      if (isNull() || dims.size() <= i) {
        return NcDim();
      } else {
        return dims[i];
      }
    }

    void putVar(double             const *data) { ncwrap( nc_put_var_double   ( ncid , id , data ) , __LINE__ ); }
    void putVar(float              const *data) { ncwrap( nc_put_var_float    ( ncid , id , data ) , __LINE__ ); }
    void putVar(int                const *data) { ncwrap( nc_put_var_int      ( ncid , id , data ) , __LINE__ ); }
    void putVar(long               const *data) { ncwrap( nc_put_var_long     ( ncid , id , data ) , __LINE__ ); }
    void putVar(long long          const *data) { ncwrap( nc_put_var_longlong ( ncid , id , data ) , __LINE__ ); }
    void putVar(signed char        const *data) { ncwrap( nc_put_var_schar    ( ncid , id , data ) , __LINE__ ); }
    void putVar(short              const *data) { ncwrap( nc_put_var_short    ( ncid , id , data ) , __LINE__ ); }
    void putVar(unsigned char      const *data) { ncwrap( nc_put_var_uchar    ( ncid , id , data ) , __LINE__ ); }
    void putVar(unsigned int       const *data) { ncwrap( nc_put_var_uint     ( ncid , id , data ) , __LINE__ ); }
    void putVar(unsigned long      const *data) { ncwrap( nc_put_var_uint     ( ncid , id , (unsigned int const *) data ) , __LINE__ ); }
    void putVar(unsigned long long const *data) { ncwrap( nc_put_var_ulonglong( ncid , id , data ) , __LINE__ ); }
    void putVar(unsigned short     const *data) { ncwrap( nc_put_var_ushort   ( ncid , id , data ) , __LINE__ ); }
    void putVar(char               const *data) { ncwrap( nc_put_var_text     ( ncid , id , data ) , __LINE__ ); }
    void putVar(bool               const *data) { throw std::runtime_error("ERROR: Cannot write bools to netCDF file"); }

    void putVar(std::vector<size_t> start , std::vector<size_t> count, double             const *data) { ncwrap( nc_put_vara_double   ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, float              const *data) { ncwrap( nc_put_vara_float    ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, int                const *data) { ncwrap( nc_put_vara_int      ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, long               const *data) { ncwrap( nc_put_vara_long     ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, long long          const *data) { ncwrap( nc_put_vara_longlong ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, signed char        const *data) { ncwrap( nc_put_vara_schar    ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, short              const *data) { ncwrap( nc_put_vara_short    ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, unsigned char      const *data) { ncwrap( nc_put_vara_uchar    ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, unsigned int       const *data) { ncwrap( nc_put_vara_uint     ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, unsigned long      const *data) { ncwrap( nc_put_vara_uint     ( ncid , id , start.data() , count.data(), (unsigned int const *) data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, unsigned long long const *data) { ncwrap( nc_put_vara_ulonglong( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, unsigned short     const *data) { ncwrap( nc_put_vara_ushort   ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, char               const *data) { ncwrap( nc_put_vara_text     ( ncid , id , start.data() , count.data(), data ) , __LINE__ ); }
    void putVar(std::vector<size_t> start , std::vector<size_t> count, bool               const *data) { throw std::runtime_error("ERROR: Cannot write bools to netCDF file"); }

    void getVar(double             *data) const { ncwrap( nc_get_var_double   ( ncid , id , data ) , __LINE__ ); }
    void getVar(float              *data) const { ncwrap( nc_get_var_float    ( ncid , id , data ) , __LINE__ ); }
    void getVar(int                *data) const { ncwrap( nc_get_var_int      ( ncid , id , data ) , __LINE__ ); }
    void getVar(long               *data) const { ncwrap( nc_get_var_long     ( ncid , id , data ) , __LINE__ ); }
    void getVar(long long          *data) const { ncwrap( nc_get_var_longlong ( ncid , id , data ) , __LINE__ ); }
    void getVar(signed char        *data) const { ncwrap( nc_get_var_schar    ( ncid , id , data ) , __LINE__ ); }
    void getVar(short              *data) const { ncwrap( nc_get_var_short    ( ncid , id , data ) , __LINE__ ); }
    void getVar(unsigned char      *data) const { ncwrap( nc_get_var_uchar    ( ncid , id , data ) , __LINE__ ); }
    void getVar(unsigned int       *data) const { ncwrap( nc_get_var_uint     ( ncid , id , data ) , __LINE__ ); }
    void getVar(unsigned long      *data) const { ncwrap( nc_get_var_uint     ( ncid , id , (unsigned int *) data ) , __LINE__ ); }
    void getVar(unsigned long long *data) const { ncwrap( nc_get_var_ulonglong( ncid , id , data ) , __LINE__ ); }
    void getVar(unsigned short     *data) const { ncwrap( nc_get_var_ushort   ( ncid , id , data ) , __LINE__ ); }
    void getVar(char               *data) const { ncwrap( nc_get_var_text     ( ncid , id , data ) , __LINE__ ); }
    void getVar(bool               *data) const { throw std::runtime_error("ERROR: Cannot read bools directly from netCDF file. This should've been intercepted and changed to int."); }

    void print() {
      std::cout << "Variable Name: " << name << "\n";
      std::cout << "Dims: \n";
      for (int i=0; i < dims.size(); i++) {
        std::cout << "  " << dims[i].getName() << ";  Size: " << dims[i].getSize() << "\n\n";
      }
    }
  };


  /** @private */
  class NcFile {
  public:
    int ncid;
    NcFile() { ncid = -999; }
    ~NcFile() {}
    NcFile(int ncid) { this->ncid = ncid; }
    NcFile(NcFile &&in) {
      this->ncid = in.ncid;
    }
    NcFile(NcFile const &in) {
      this->ncid = in.ncid;
    }
    NcFile &operator=(NcFile &&in) {
      this->ncid = in.ncid;
      return *this;
    }
    NcFile &operator=(NcFile const &in) {
      this->ncid = in.ncid;
      return *this;
    }

    bool isNull() { return ncid == -999; }

    void open( const std::string& fname , int mode ) {
      close();
      if (! (mode == NETCDF_MODE_READ || mode == NETCDF_MODE_WRITE) ) {
        throw std::runtime_error("ERROR: open mode can be NETCDF_MODE_READ or NETCDF_MODE_WRITE");
      }
      ncwrap( nc_open( fname.c_str() , mode , &ncid ) , __LINE__ );
    }

    void create( const std::string& fname , int mode ) {
      close();
      if (! (mode == NETCDF_MODE_NEW || mode == NETCDF_MODE_REPLACE) ) {
        throw std::runtime_error("ERROR: open mode can be NETCDF_MODE_NEW or NETCDF_MODE_REPLACE");
      }
      ncwrap( nc_create( fname.c_str() , mode | NC_NETCDF4 , &ncid ) , __LINE__ );
    }

    void close() {
      if (ncid != -999) ncwrap( nc_close( ncid ) , __LINE__ );
      ncid = -999;
    }

    NcVar getVar( std::string varName ) const {
      int varid;
      int ierr = nc_inq_varid( ncid , varName.c_str() , &varid);
      if (ierr != NC_NOERR) return NcVar();
      char vname[NC_MAX_NAME+1];
      int  type;
      int  ndims;
      int  dimids[NC_MAX_VAR_DIMS];
      int  natts;
      // Get variable information
      ncwrap( nc_inq_var(ncid , varid , vname , &type , &ndims , dimids , &natts ) , __LINE__ );
      // Accumulate the dimensions
      std::vector<NcDim> dims(ndims);
      for (int i=0; i < ndims; i++) {
        dims[i] = getDim( dimids[i] );
      }
      return NcVar( ncid , varName , dims , varid , type );
    }

    NcDim getDim( std::string dimName ) const {
      int dimid;
      int ierr = nc_inq_dimid( ncid , dimName.c_str() , &dimid);
      if (ierr != NC_NOERR) return NcDim();
      return getDim( dimid );
    }

    NcDim getDim( int dimid ) const {
      char   dname[NC_MAX_NAME+1];
      size_t len;
      int    unlim_dimid;
      ncwrap( nc_inq_dim( ncid , dimid , dname , &len ) , __LINE__ );
      ncwrap( nc_inq_unlimdim( ncid , &unlim_dimid ) , __LINE__ );
      return NcDim( std::string(dname) , len , dimid , dimid == unlim_dimid );
    }

    NcVar addVar( std::string varName , int type , std::vector<NcDim> &dims ) {
      std::vector<int> dimids(dims.size());
      for (int i=0; i < dims.size(); i++) { dimids[i] = dims[i].getId(); }
      int varid;
      ncwrap( nc_def_var(ncid , varName.c_str() , type , dims.size() , dimids.data() , &varid) , __LINE__ );
      return NcVar( ncid , varName , dims , varid , type );
    }

    NcVar addVar( std::string varName , int type ) {
      int varid;
      int *dummy = nullptr;
      ncwrap( nc_def_var(ncid , varName.c_str() , type , 0 , dummy , &varid) , __LINE__ );
      return NcVar( ncid , varName , std::vector<NcDim>(0) , varid , type );
    }

    NcDim addDim( std::string dimName , size_t len ) {
      int dimid;
      ncwrap( nc_def_dim(ncid , dimName.c_str() , len , &dimid ) , __LINE__ );
      return NcDim( dimName , len , dimid , false );
    }

    NcDim addDim( std::string dimName ) {
      int dimid;
      ncwrap( nc_def_dim(ncid , dimName.c_str() , NC_UNLIMITED , &dimid ) , __LINE__ );
      return NcDim( dimName , 0 , dimid , true );
    }

  };


  /** @private */
  NcFile file;


  SimpleNetCDF() { }


  /** @brief Files are automatically closed when SimpleNetCDF objects are destroyed */
  ~SimpleNetCDF() { close(); }


  /** @brief Open a netcdf file
   * @param mode Can be NETCDF_MODE_READ or NETCDF_MODE_WRITE */
  void open(const std::string& fname , int mode = NETCDF_MODE_READ) { file.open(fname,mode); }


  /** @brief Create a netcdf file
   * @param mode Can be NETCDF_MODE_CLOBBER or NETCDF_MODE_NOCLOBBER */
  void create(const std::string& fname , int mode = NC_CLOBBER) { file.create(fname,mode); }


  /** @brief Close the netcdf file */
  void close() { file.close(); }


  /** @brief Determine if a variable name exists */
  bool varExists( std::string varName ) const { return ! file.getVar(varName).isNull(); }


  /** @brief Determine if a dimension name exists */
  bool dimExists( std::string dimName ) const { return ! file.getDim(dimName).isNull(); }


  /** @brief Determine the size of a dimension name */
  size_t getDimSize( std::string dimName ) const { return file.getDim(dimName).getSize(); }


  /** @brief Create a dimension of the given length */
  void createDim( std::string dimName , size_t len ) { file.addDim( dimName , len ); }

  /** @brief Create an unlimited dimension */
  void createDim( std::string dimName ) { file.addDim( dimName ); }


  /** @brief Write an entire Array at once */
  template <class View>
  void write(const View& arr, std::string varName , std::vector<std::string> dimNames) {
    using myStyle = typename View::array_layout;
    using myMem   = typename View::memory_space;
    using T       = typename View::non_const_value_type;
    constexpr bool is_c_layout   = std::is_same_v<myStyle, Kokkos::LayoutRight>;
    constexpr bool is_device_mem = !std::is_same_v<myMem, Kokkos::DefaultHostExecutionSpace::memory_space>;
    constexpr auto rank = View::rank;

    if (rank != dimNames.size()) { throw std::runtime_error("dimNames.size() != Array's rank"); }
    std::vector<NcDim> dims(rank); // List of dimensions for this variable
    // Make sure the dimensions are in there and are the right sizes
    for (int i=0; i<rank; i++) {
      auto dimLoc = file.getDim( dimNames[i] );
      // If dimension doesn't exist, create it; otherwise, make sure it's the right size
      NcDim tmp;
      if ( dimLoc.isNull() ) {
        tmp = file.addDim( dimNames[i] , arr.extent(i) );
      } else {
        if (dimLoc.getSize() != arr.extent(i)) {
          throw std::runtime_error("dimension size differs from the file");
        }
        tmp = dimLoc;
      }
      if (is_c_layout) {
        dims[i] = tmp;
      } else {
        dims[rank-1-i] = tmp;
      }
    }
    // Make sure the variable is there and is the right dimension
    auto var = file.getVar(varName);
    if ( var.isNull() ) {
      var = file.addVar( varName , getType<T>() , dims );
    } else {
      if ( var.getType() != getType<T>() ) { throw std::runtime_error("Existing variable's type != array's type"); }
      auto varDims = var.getDims();
      if (varDims.size() != rank) { throw std::runtime_error("Existing variable's rank != array's rank"); }
      for (int i=0; i < varDims.size(); i++) {
        if (is_c_layout) {
          if (varDims[i].getSize() != arr.extent(i)) {
            throw std::runtime_error("Existing variable's dimension sizes are not the same as the array's");
          }
        } else {
          if (varDims[rank-1-i].getSize() != arr.extent(i)) {
            throw std::runtime_error("Existing variable's dimension sizes are not the same as the array's");
          }
        }
      }
    }

    if (is_device_mem) {
      auto hv = Kokkos::create_mirror_view_and_copy(HostDevice(), arr);
      var.putVar(hv.data());
    } else {
      var.putVar(arr.data());
    }
  }


  /** @brief Write one entry of a scalar into the unlimited index */
  template <class T, typename std::enable_if<std::is_arithmetic<T>::value,int>::type = 0 >
  void write1(T val , std::string varName , int ind , std::string ulDimName="unlim" ) {
    // Get the unlimited dimension or create it if it doesn't exist
    auto ulDim = file.getDim( ulDimName );
    if ( ulDim.isNull() ) {
      ulDim = file.addDim( ulDimName );
    }
    // Make sure the variable is there and is the right dimension
    auto var = file.getVar(varName);
    if ( var.isNull() ) {
      std::vector<NcDim> dims(1);
      dims[0] = ulDim;
      var = file.addVar( varName , getType<T>() , dims );
    }
    std::vector<size_t> start(1);
    std::vector<size_t> count(1);
    start[0] = ind;
    count[0] = 1;
    var.putVar(start,count,&val);
  }


  /** @brief Write one entry of an Array into the unlimited index */
  template <typename View>
  void write1(const View& arr , std::string varName , std::vector<std::string> dimNames ,
              int ind , std::string ulDimName="unlim" ) {
    using myStyle = typename View::array_layout;
    using myMem   = typename View::memory_space;
    using T       = typename View::non_const_value_type;
    constexpr bool is_c_layout = std::is_same_v<myStyle, Kokkos::LayoutRight>;
    constexpr bool is_device_mem = !std::is_same_v<myMem, Kokkos::DefaultHostExecutionSpace::memory_space>;
    constexpr auto rank = View::rank;

    if (rank != dimNames.size()) { throw std::runtime_error("dimNames.size() != Array's rank"); }
    std::vector<NcDim> dims(rank+1); // List of dimensions for this variable
    // Get the unlimited dimension or create it if it doesn't exist
    dims[0] = file.getDim( ulDimName );
    if ( dims[0].isNull() ) {
      dims[0] = file.addDim( ulDimName );
    }
    // Make sure the dimensions are in there and are the right sizes
    for (int i=0; i<rank; i++) {
      auto dimLoc = file.getDim( dimNames[i] );
      // If dimension doesn't exist, create it; otherwise, make sure it's the right size
      NcDim tmp;
      if ( dimLoc.isNull() ) {
        tmp = file.addDim( dimNames[i] , arr.extent(i) );
      } else {
        if (dimLoc.getSize() != arr.extent(i)) {
          throw std::runtime_error("dimension size differs from the file");
        }
        tmp = dimLoc;
      }
      if (is_c_layout) {
        dims[i] = tmp;
      } else {
        dims[rank-1-i] = tmp;
      }
    }
    // Make sure the variable is there and is the right dimension
    auto var = file.getVar(varName);
    if ( var.isNull() ) {
      var = file.addVar( varName , getType<T>() , dims );
    } else {
      if ( var.getType() != getType<T>() ) { throw std::runtime_error("Existing variable's type != array's type"); }
      auto varDims = var.getDims();
      if (varDims.size() != rank) {
        throw std::runtime_error("Existing variable's rank != array's rank");
      }
      for (int i=0; i < varDims.size(); i++) {
        if (is_c_layout) {
          if (varDims[i].getSize() != arr.extent(i)) {
            throw std::runtime_error("Existing variable's dimension sizes are not the same as the array's");
          }
        } else {
          if (varDims[rank-1-i].getSize() != arr.extent(i)) {
            throw std::runtime_error("Existing variable's dimension sizes are not the same as the array's");
          }
        }
      }
    }

    std::vector<size_t> start(rank+1);
    std::vector<size_t> count(rank+1);
    start[0] = ind;
    count[0] = 1;
    for (int i=1; i < rank+1; i++) {
      start[i] = 0;
      count[i] = dims[i].getSize();
    }
    if (is_device_mem) {
      auto hv = Kokkos::create_mirror_view_and_copy(HostDevice(), arr);
      var.putVar(start,count,hv.data());
    } else {
      var.putVar(start,count,arr.data());
    }
  }


  /** @brief Read an entire Array */
  template <typename View,
            typename std::enable_if<is_view_v<View>>::type* = nullptr>
  void read(View& arr , std::string varName) {
    using myStyle = typename View::array_layout;
    using T       = typename View::non_const_value_type;

    using LeftHostView = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft, HostDevice>;
    constexpr auto rank = View::rank;

    // Make sure the variable is there and is the right dimension
    auto var = file.getVar(varName);
    std::vector<int> dimSizes(rank);
    if ( ! var.isNull() ) {
      auto varDims = var.getDims();
      if (varDims.size() != rank) { throw std::runtime_error("Existing variable's rank != array's rank"); }
      for (int i=0; i < varDims.size(); i++) { dimSizes[i] = varDims[varDims.size()-1-i].getSize(); }
      bool createArr = false;
      for (int i=0; i < dimSizes.size(); i++) {
        if (dimSizes[i] != arr.extent(i)) {
          createArr = true;
          break;
        }
      }
      if (createArr) {
        arr = View("read arr", get_layout<myStyle>(dimSizes));
      }
    } else { throw std::runtime_error("Variable does not exist"); }

    Kokkos::LayoutLeft llayout;
    for (auto r = 0; r < rank; ++r) {
      llayout.dimension[r] = arr.layout().dimension[r];
    }
    LeftHostView read_data("read_data", llayout);

    if (std::is_same_v<T,bool>) {
      int* tmp = new int[arr.size()];
      var.getVar(tmp);
      for (size_t i=0; i < arr.size(); ++i) { read_data.data()[i] = (tmp[i] == 1); }
      delete[] tmp;
    }
    else {
      var.getVar(read_data.data());
      // integer data is nearly always idx data, so adjust it to 0-based
      if (std::is_same_v<T,int>) {
        for (size_t i=0; i < arr.size(); ++i) { read_data.data()[i] -= 1; }
      }
    }
    auto arr_mirror = Kokkos::create_mirror_view(arr);
    Kokkos::deep_copy(arr_mirror, read_data);
    Kokkos::deep_copy(arr, arr_mirror);
  }


  /** @brief Read a single scalar value */
  template <class T,
            typename std::enable_if<!is_view_v<T>>::type* = nullptr>
  void read(T &arr , std::string varName) {
    auto var = file.getVar(varName);
    if ( var.isNull() ) { throw std::runtime_error("Variable does not exist"); }
    var.getVar(&arr);
  }


  /** @brief Write a single scalar value */
  template <class T,
            typename std::enable_if<!is_view_v<T>>::type* = nullptr>
  void write(T arr , std::string varName) {
    auto var = file.getVar(varName);
    if ( var.isNull() ) {
      var = file.addVar( varName , getType<T>() );
    }
    var.putVar(&arr);
  }


  /** @private */
  template <class T> int getType() const {
    if ( std::is_same_v<typename std::remove_cv<T>::type,signed        char> ) { return NC_BYTE;   }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,unsigned      char> ) { return NC_UBYTE;  }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,             short> ) { return NC_SHORT;  }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,unsigned     short> ) { return NC_USHORT; }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,               int> ) { return NC_INT;    }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,unsigned       int> ) { return NC_UINT;   }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,              long> ) { return NC_INT;    }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,unsigned      long> ) { return NC_UINT;   }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,         long long> ) { return NC_INT64;  }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,unsigned long long> ) { return NC_UINT64; }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,             float> ) { return NC_FLOAT;  }
    else if ( std::is_same_v<typename std::remove_cv<T>::type,            double> ) { return NC_DOUBLE; }
    if ( std::is_same_v<typename std::remove_cv<T>::type,              char> ) { return NC_CHAR;   }
    else { throw std::runtime_error("Invalid type"); }
    return -1;
  }
};

// Copied with minor modes from YAKL_random.h

class Random {
 protected:
  /** @private */
  typedef unsigned long long u8;
  /** @private */
  u8 static constexpr rot(u8 x, u8 k) { return (((x)<<(k))|((x)>>(64-(k)))); }
  /** @private */
  struct State { u8 a, b, c, d; };
  /** @private */
  State state;

 public:

  /** @brief Initializes a prng object with the seed 1368976481. Warm-up of 20 iterations. */
  KOKKOS_INLINE_FUNCTION Random()                            { set_seed(1368976481L); } // I made up this number
  /** @brief Initializes a prng object with the specified seed. Warm-up of 20 iterations. */
  KOKKOS_INLINE_FUNCTION Random(u8 seed)                     { set_seed(seed); }
  /** @brief Copies a Random object */
  KOKKOS_INLINE_FUNCTION Random(Random const            &in) { this->state = in.state; }
  /** @brief Moves a Random object */
  KOKKOS_INLINE_FUNCTION Random(Random                 &&in) { this->state = in.state; }
  /** @brief Copies a Random object */
  KOKKOS_INLINE_FUNCTION Random &operator=(Random const &in) { this->state = in.state; return *this; }
  /** @brief Moves a Random object */
  KOKKOS_INLINE_FUNCTION Random &operator=(Random      &&in) { this->state = in.state; return *this; }

  /** @brief Assigns a seed. Warm-up of 20 iterations. */
  KOKKOS_INLINE_FUNCTION void set_seed(u8 seed) {
    state.a = 0xf1ea5eed;  state.b = seed;  state.c = seed;  state.d = seed;
    for (int i=0; i<20; ++i) { gen(); }
  }

  /** @brief Generates a random unsigned integer between zero and `std::numeric_limits<u8>::max() - 1` */
  KOKKOS_INLINE_FUNCTION u8 gen() {
    u8 e    = state.a - rot(state.b, 7);
    state.a = state.b ^ rot(state.c,13);
    state.b = state.c + rot(state.d,37);
    state.c = state.d + e;
    state.d = e       + state.a;
    return state.d;
  }

  /** @brief Generates a random floating point value between `0` and `1`
   * @param T The type of the floating point number */
  template <class T> KOKKOS_INLINE_FUNCTION T genFP() {
    return static_cast<T>(gen()) / static_cast<T>(std::numeric_limits<u8>::max());
  }

  /** @brief Generates a random floating point value between `lb` and `ub`
   * @param T  The type of the floating point number
   * @param lb Lower bound of the random number
   * @param ub Upper bound of the random number*/
  template <class T> KOKKOS_INLINE_FUNCTION T genFP(T lb, T ub) {
    return genFP<T>() * (ub-lb) + lb;
  }

};

#endif

}
