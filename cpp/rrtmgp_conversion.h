#pragma once

#include "rrtmgp_const.h"

#include <netcdf.h>
#include <stdexcept>

#ifdef RRTMGP_ENABLE_KOKKOS

/**
 * Helper functions for the conversion to Kokkos
 */

namespace conv {

template <typename T>
struct is_view
{
  static constexpr bool value = Kokkos::is_view<T>::value || Kokkos::Experimental::is_offset_view<T>::value;
};

template <class T>
inline constexpr bool is_view_v = is_view<T>::value;

// Copied with minor mods from YAKL intrinsics
template <class T,
          typename std::enable_if<!is_view_v<T>>::type* = nullptr>
KOKKOS_INLINE_FUNCTION
T constexpr epsilon(T) { return std::numeric_limits<T>::epsilon(); }

template <class View,
          typename std::enable_if<is_view_v<View>>::type* = nullptr>
KOKKOS_INLINE_FUNCTION
typename View::non_const_value_type constexpr epsilon(const View& arr) { return std::numeric_limits<typename View::non_const_value_type>::epsilon(); }

// These are for debugging
template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p1d(const KView& view, const std::string& name, int idx)
{ std::cout << "JGFK " << name << "(" << idx << ") = " << view(idx) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p1d(const YArray& array, const std::string& name, int idx)
{
  const int adjust_val = std::is_same<typename YArray::non_const_value_type, int>::value ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx-1 << ") = " << array(idx) - adjust_val << std::endl; }

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p2d(const KView& view, const std::string& name, int idx1, int idx2)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ") = " << view(idx1, idx2) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p2d(const YArray& array, const std::string& name, int idx1, int idx2)
{
  const int adjust_val = std::is_same<typename YArray::non_const_value_type, int>::value ? 1 : 0;
  std::cout << "JGFY " << name << "(" << idx1-1 << ", " << idx2-1 << ") = " << array(idx1, idx2) - adjust_val << std::endl; }

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void p3d(const KView& view, const std::string& name, int idx1, int idx2, int idx3)
{ std::cout << "JGFK " << name << "(" << idx1 << ", " << idx2 << ", " << idx3 << ") = " << view(idx1, idx2, idx3) << std::endl; }

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void p3d(const YArray& array, const std::string& name, int idx1, int idx2, int idx3)
{
  const int adjust_val = std::is_same<typename YArray::non_const_value_type, int>::value ? 1 : 0;
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
  const int adjust_val = std::is_same<typename YArray::non_const_value_type, int>::value ? 1 : 0;
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
  const int adjust_val = std::is_same<typename YArray::non_const_value_type, int>::value ? 1 : 0;
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
  const int adjust_val = std::is_same<typename YArray::non_const_value_type, int>::value ? 1 : 0;
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

#define RRT_REQUIRE(condition, msg) IMPL_THROW_RRT(condition, msg, std::runtime_error)

// Copied from EKAT
template<typename T, int N>
struct DataND {
  using type = typename DataND<T,N-1>::type*;
};
template<typename T>
struct DataND<T,0> {
  using type = T;
};

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

template <typename T>
inline
bool approx_eq(const T lhs, const T rhs)
{
  return lhs == rhs;
}

template <>
inline
bool approx_eq<real>(const real lhs, const real rhs)
{
  constexpr real tol = 1e-12;
  return std::abs(lhs - rhs) < tol;
}

template <typename YArray, typename KView>
void compare_yakl_to_kokkos(const YArray& yarray, const KView& kview, bool index_data=false)
{
  using yakl::intrinsics::size;

  constexpr auto krank = KView::rank;
  const auto yrank = yarray.get_rank();

  RRT_REQUIRE(krank == yrank, "Rank mismatch for: " << kview.label());

  auto hkview = Kokkos::create_mirror_view(kview);
  Kokkos::deep_copy(hkview, kview);
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
  static constexpr auto yakl_mem = std::is_same<typename KView::device_type, HostDevice>::value ? yakl::memHost : yakl::memDevice;
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

template <typename T>
struct LTFunc
{
  T m_val;
  LTFunc(T val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(const T& val) const
  {
    return val < m_val;
  }
};

template <typename KView, typename Functor>
bool any(const KView& view, const Functor& functor)
{
  using exe_space_t = typename KView::execution_space;

  bool rv = false;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exe_space_t>(0, view.size()),
    KOKKOS_LAMBDA(size_t i, bool& val) {
      val = functor(view.data()[i]);
    }, Kokkos::BOr<bool>(rv));

  return rv;
}

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

template <typename KView>
typename KView::non_const_value_type sum(const KView& view)
{
  using scalar_t    = typename KView::non_const_value_type;
  using exe_space_t = typename KView::execution_space;
  using sum_t       = std::conditional_t<std::is_same<scalar_t, bool>::value, int, scalar_t>;

  sum_t rv;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<exe_space_t>(0, view.size()),
    KOKKOS_LAMBDA(size_t i, sum_t& lsum) {
      lsum += view.data()[i];
    }, Kokkos::Sum<sum_t>(rv));
  return rv;
}

template <typename KView,
          typename std::enable_if<is_view_v<KView>>::type* = nullptr>
void print(const std::string& name, const KView& view)
{
  for (size_t i = 0; i < view.size(); ++i) {
    std::cout << "JGFK " << name << "(" << i << ") = " << view.data()[i] << std::endl;
  }
}

template <typename YArray,
          typename std::enable_if<!is_view_v<YArray>>::type* = nullptr>
void print(const std::string& name, const YArray& array)
{
  for (size_t i = 0; i < array.totElems(); ++i) {
    std::cout << "JGFY " << name << "(" << i << ") = " << array.data()[i] << std::endl;
  }
}


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

    void open( std::string fname , int mode ) {
      close();
      if (! (mode == NETCDF_MODE_READ || mode == NETCDF_MODE_WRITE) ) {
        throw std::runtime_error("ERROR: open mode can be NETCDF_MODE_READ or NETCDF_MODE_WRITE");
      }
      ncwrap( nc_open( fname.c_str() , mode , &ncid ) , __LINE__ );
    }

    void create( std::string fname , int mode ) {
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
  void open(std::string fname , int mode = NETCDF_MODE_READ) { file.open(fname,mode); }


  /** @brief Create a netcdf file
   * @param mode Can be NETCDF_MODE_CLOBBER or NETCDF_MODE_NOCLOBBER */
  void create(std::string fname , int mode = NC_CLOBBER) { file.create(fname,mode); }


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
    constexpr bool is_c_layout   = std::is_same<myStyle, Kokkos::LayoutRight>::value;
    constexpr bool is_device_mem = !std::is_same<myMem, Kokkos::DefaultHostExecutionSpace::memory_space>::value;
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
      auto hv = Kokkos::create_mirror_view(arr);
      Kokkos::deep_copy(hv, arr);
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
    constexpr bool is_c_layout = std::is_same<myStyle, Kokkos::LayoutRight>::value;
    constexpr bool is_device_mem = !std::is_same<myMem, Kokkos::DefaultHostExecutionSpace::memory_space>::value;
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
      auto hv = Kokkos::create_mirror_view(arr);
      Kokkos::deep_copy(hv, arr);
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
    using myMem   = typename View::memory_space;
    using T       = typename View::non_const_value_type;
    constexpr bool is_c_layout   = std::is_same<myStyle, Kokkos::LayoutRight>::value;
    constexpr bool is_device_mem = !std::is_same<myMem, Kokkos::DefaultHostExecutionSpace::memory_space>::value;
    constexpr auto rank = View::rank;

    // Make sure the variable is there and is the right dimension
    auto var = file.getVar(varName);
    std::vector<int> dimSizes(rank);
    if ( ! var.isNull() ) {
      auto varDims = var.getDims();
      if (varDims.size() != rank) { throw std::runtime_error("Existing variable's rank != array's rank"); }
      if (is_c_layout) {
        for (int i=0; i < varDims.size(); i++) { dimSizes[i] = varDims[i].getSize(); }
      } else {
        for (int i=0; i < varDims.size(); i++) { dimSizes[i] = varDims[varDims.size()-1-i].getSize(); }
      }
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

    if (is_device_mem) {
      auto arrHost = Kokkos::create_mirror_view(arr);
      if (std::is_same<T,bool>::value) {
        int* tmp = new int[arr.size()];
        var.getVar(tmp);
        for (size_t i=0; i < arr.size(); ++i) { arrHost.data()[i] = (tmp[i] == 1); }
        delete[] tmp;
      }
      else {
        var.getVar(arrHost.data());
        // integer data is nearly always idx data, so adjust it to 0-based
        if (std::is_same<T,int>::value) {
          for (size_t i=0; i < arr.size(); ++i) { arrHost.data()[i] -= 1; }
        }
      }
      Kokkos::deep_copy(arr, arrHost);
    } else {
      if (std::is_same<T,bool>::value) {
        int* tmp = new int[arr.size()];
        var.getVar(tmp);
        for (size_t i=0; i < arr.size(); ++i) { arr.data()[i] = (tmp[i] == 1); }
        delete[] tmp;
      } else {
        var.getVar(arr.data());
        // integer data is nearly always idx data, so adjust it to 0-based
        if (std::is_same<T,int>::value) {
          for (size_t i=0; i < arr.size(); ++i) { arr.data()[i] -= 1; }
        }
      }
    }
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
    if ( std::is_same<typename std::remove_cv<T>::type,signed        char>::value ) { return NC_BYTE;   }
    else if ( std::is_same<typename std::remove_cv<T>::type,unsigned      char>::value ) { return NC_UBYTE;  }
    else if ( std::is_same<typename std::remove_cv<T>::type,             short>::value ) { return NC_SHORT;  }
    else if ( std::is_same<typename std::remove_cv<T>::type,unsigned     short>::value ) { return NC_USHORT; }
    else if ( std::is_same<typename std::remove_cv<T>::type,               int>::value ) { return NC_INT;    }
    else if ( std::is_same<typename std::remove_cv<T>::type,unsigned       int>::value ) { return NC_UINT;   }
    else if ( std::is_same<typename std::remove_cv<T>::type,              long>::value ) { return NC_INT;    }
    else if ( std::is_same<typename std::remove_cv<T>::type,unsigned      long>::value ) { return NC_UINT;   }
    else if ( std::is_same<typename std::remove_cv<T>::type,         long long>::value ) { return NC_INT64;  }
    else if ( std::is_same<typename std::remove_cv<T>::type,unsigned long long>::value ) { return NC_UINT64; }
    else if ( std::is_same<typename std::remove_cv<T>::type,             float>::value ) { return NC_FLOAT;  }
    else if ( std::is_same<typename std::remove_cv<T>::type,            double>::value ) { return NC_DOUBLE; }
    if ( std::is_same<typename std::remove_cv<T>::type,              char>::value ) { return NC_CHAR;   }
    else { throw std::runtime_error("Invalid type"); }
    return -1;
  }
};

}

#endif
