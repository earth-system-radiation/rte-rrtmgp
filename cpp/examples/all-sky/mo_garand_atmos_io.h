
#pragma once
#include "rrtmgp_const.h"
#include "mo_gas_concentrations.h"
#include "rrtmgp_conversion.h"

#ifdef RRTMGP_ENABLE_YAKL
void read_atmos(std::string input_file, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, real2d &col_dry, int ncol);


void write_sw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, real2d const &flux_dir, int ncol);


void write_lw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, int ncol);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
template <typename ViewT, typename GasConcsT>
void read_atmos(const std::string& input_file, ViewT &p_lay, ViewT &t_lay, ViewT &p_lev, ViewT &t_lev,
                GasConcsT &gas_concs, ViewT &col_dry, int ncol)
{
  using DeviceT = typename ViewT::device_type;
  using LayoutT = typename ViewT::array_layout;
  using RealT   = typename ViewT::non_const_value_type;
  using MDRP    = typename conv::MDRP<LayoutT>;

  conv::SimpleNetCDF io;
  io.open(input_file, conv::NETCDF_MODE_READ);

  int nlay = io.getDimSize("lay");
  int nlev = io.getDimSize("lev");

  p_lay = ViewT("p_lay",ncol,nlay);
  t_lay = ViewT("t_lay",ncol,nlay);
  p_lev = ViewT("p_lev",ncol,nlev);
  t_lev = ViewT("t_lev",ncol,nlev);

  ViewT tmp2d;
  // p_lay
  io.read(tmp2d,"p_lay");
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay, icol, ilay,
    p_lay(icol,ilay) = tmp2d(0,ilay);
  ));
  // t_lay
  io.read(tmp2d,"t_lay");
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay, icol, ilay,
    t_lay(icol,ilay) = tmp2d(0,ilay);
  ));
  // p_lev
  io.read(tmp2d,"p_lev");
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlev, icol, ilev,
    p_lev(icol,ilev) = tmp2d(0,ilev);
  ));
  // t_lev
  io.read(tmp2d,"t_lev");
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlev, icol, ilev,
    t_lev(icol,ilev) = tmp2d(0,ilev);
  ));

  std::vector<std::string> gas_names = {
    "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"
  };

  // Initialize GasConcs object with an "ncol" given from the calling program
  gas_concs.init(gas_names,ncol,nlay);

  tmp2d = ViewT();     // Reset the tmp2d variable
  for (auto& gas_name : gas_names) {
    std::string vmr_name = "vmr_"+gas_name;
    if ( ! io.varExists(vmr_name) ) { stoprun("ERROR: gas does not exist in input file"); }
    // Read in 2-D varaible
    io.read(tmp2d,vmr_name);
    // Create 1-D variable with just the first column
    Kokkos::View<RealT*, LayoutT, DeviceT> tmp1d("tmp1d",nlay);
    TIMED_KERNEL(Kokkos::parallel_for( nlay, KOKKOS_LAMBDA (int i) {
      tmp1d(i) = tmp2d(0,i);
    }));
    // Call set_vmr with only the first column from the data file copied among all of the model columns
    gas_concs.set_vmr( gas_name , tmp1d );
  }

  if ( io.varExists("col_dry") ) {
    col_dry = ViewT("col_dry",ncol,nlay);
    tmp2d = ViewT();     // Reset the tmp2d variable
    io.read(tmp2d,"col_dry");
    TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay, icol, ilay,
      col_dry(icol,ilay) = tmp2d(0,ilay);
    ));
  }

  io.close();
}

template <typename FluxupT, typename FluxdnT, typename FluxdirT>
void write_sw_fluxes(const std::string& fileName, FluxupT const &flux_up, FluxdnT const &flux_dn, FluxdirT const &flux_dir, int ncol)
{
  conv::SimpleNetCDF io;
  io.open(fileName , conv::NETCDF_MODE_WRITE);
  io.write(flux_up  , "sw_flux_up_result"  , {"col_new","lev"});
  io.write(flux_dn  , "sw_flux_dn_result"  , {"col_new","lev"});
  io.write(flux_dir , "sw_flux_dir_result" , {"col_new","lev"});
  io.close();
}

template <typename FluxupT, typename FluxdnT>
void write_lw_fluxes(const std::string& fileName, FluxupT const &flux_up, FluxdnT const &flux_dn, int ncol)
{
  conv::SimpleNetCDF io;
  io.open(fileName , conv::NETCDF_MODE_WRITE);
  io.write(flux_up , "lw_flux_up_result" , {"col_new","lev"});
  io.write(flux_dn , "lw_flux_dn_result" , {"col_new","lev"});
  io.close();
}
#endif
