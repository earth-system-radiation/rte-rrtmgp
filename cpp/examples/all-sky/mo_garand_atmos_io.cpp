
#include "mo_garand_atmos_io.h"
#include "rrtmgp_conversion.h"

void read_atmos(std::string input_file, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, real2d &col_dry, int ncol) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  yakl::SimpleNetCDF io;
  io.open(input_file , yakl::NETCDF_MODE_READ);

  int nlay = io.getDimSize("lay");
  int nlev = io.getDimSize("lev");

  p_lay = real2d("p_lay",ncol,nlay);
  t_lay = real2d("t_lay",ncol,nlay);
  p_lev = real2d("p_lev",ncol,nlev);
  t_lev = real2d("t_lev",ncol,nlev);

  real2d tmp2d;
  // p_lay
  io.read(tmp2d,"p_lay");
  // for (int ilay=1 ; ilay <= nlay ; ilay++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
    p_lay(icol,ilay) = tmp2d(1,ilay);
  });
  // t_lay
  io.read(tmp2d,"t_lay");
  // for (int ilay=1 ; ilay <= nlay ; ilay++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA ( int ilay, int icol) {
    t_lay(icol,ilay) = tmp2d(1,ilay);
  });
  // p_lev
  tmp2d = real2d();  // Reset tmp2d to avoid warnings about reallocating during file read
  io.read(tmp2d,"p_lev");
  // for (int ilev=1 ; ilev <= nlev ; ilev++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlev,ncol) , YAKL_LAMBDA ( int ilev, int icol) {
    p_lev(icol,ilev) = tmp2d(1,ilev);
  });
  // t_lev
  io.read(tmp2d,"t_lev");
  // for (int ilev=1 ; ilev <= nlev ; ilev++) {
  //   for (int icol=1 ; icol <= ncol ; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlev,ncol) , YAKL_LAMBDA( int ilev, int icol) {
    t_lev(icol,ilev) = tmp2d(1,ilev);
  });

  int ngas = 8;
  std::vector<std::string> gas_names = {
    "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"
  };

  // Initialize GasConcs object with an "ncol" given from the calling program
  gas_concs.init(gas_names,ncol,nlay);

  tmp2d = real2d();     // Reset the tmp2d variable
  for (auto& gas_name : gas_names) {
    std::string vmr_name = "vmr_"+gas_name;
    if ( ! io.varExists(vmr_name) ) { stoprun("ERROR: gas does not exist in input file"); }
    // Read in 2-D varaible
    io.read(tmp2d,vmr_name);
    // Create 1-D variable with just the first column
    real1d tmp1d("tmp1d",nlay);
    // for (int i=1 ; i <= nlay ; i++) {
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>(nlay) , YAKL_LAMBDA (int i) {
      tmp1d(i) = tmp2d(1,i);
    });
    // Call set_vmr with only the first column from the data file copied among all of the model columns
    gas_concs.set_vmr( gas_name , tmp1d );
  }

  if ( io.varExists("col_dry") ) {
    col_dry = real2d("col_dry",ncol,nlay);
    tmp2d = real2d();     // Reset the tmp2d variable
    io.read(tmp2d,"col_dry");
    // for (int ilay=1 ; ilay <= nlay ; ilay++) {
    //   for (int icol=1 ; icol <= ncol ; icol++) {
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA( int ilay, int icol) {
      col_dry(icol,ilay) = tmp2d(1,ilay);
    });
  }

  io.close();
}

#ifdef RRTMGP_ENABLE_KOKKOS
void read_atmos(std::string input_file, real2dk &p_lay, real2dk &t_lay, real2dk &p_lev, real2dk &t_lev,
                GasConcsK &gas_concs, real2dk &col_dry, int ncol) {
  conv::SimpleNetCDF io;
  io.open(input_file , conv::NETCDF_MODE_READ);

  int nlay = io.getDimSize("lay");
  int nlev = io.getDimSize("lev");

  p_lay = real2dk("p_lay",ncol,nlay);
  t_lay = real2dk("t_lay",ncol,nlay);
  p_lev = real2dk("p_lev",ncol,nlev);
  t_lev = real2dk("t_lev",ncol,nlev);

  real2dk tmp2d;
  // p_lay
  io.read(tmp2d,"p_lay");
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}), KOKKOS_LAMBDA (int ilay, int icol) {
    p_lay(icol,ilay) = tmp2d(0,ilay);
  });
  // t_lay
  io.read(tmp2d,"t_lay");
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}), KOKKOS_LAMBDA (int ilay, int icol) {
    t_lay(icol,ilay) = tmp2d(0,ilay);
  });
  // p_lev
  io.read(tmp2d,"p_lev");
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlev,ncol}), KOKKOS_LAMBDA (int ilev, int icol) {
    p_lev(icol,ilev) = tmp2d(0,ilev);
  });
  // t_lev
  io.read(tmp2d,"t_lev");
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlev,ncol}), KOKKOS_LAMBDA (int ilev, int icol) {
    t_lev(icol,ilev) = tmp2d(0,ilev);
  });

  int ngas = 8;
  std::vector<std::string> gas_names = {
    "h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"
  };

  // Initialize GasConcs object with an "ncol" given from the calling program
  gas_concs.init(gas_names,ncol,nlay);

  tmp2d = real2dk();     // Reset the tmp2d variable
  for (auto& gas_name : gas_names) {
    std::string vmr_name = "vmr_"+gas_name;
    if ( ! io.varExists(vmr_name) ) { stoprun("ERROR: gas does not exist in input file"); }
    // Read in 2-D varaible
    io.read(tmp2d,vmr_name);
    // Create 1-D variable with just the first column
    real1dk tmp1d("tmp1d",nlay);
    Kokkos::parallel_for( nlay, KOKKOS_LAMBDA (int i) {
      tmp1d(i) = tmp2d(0,i);
    });
    // Call set_vmr with only the first column from the data file copied among all of the model columns
    gas_concs.set_vmr( gas_name , tmp1d );
  }

  if ( io.varExists("col_dry") ) {
    col_dry = real2dk("col_dry",ncol,nlay);
    tmp2d = real2dk();     // Reset the tmp2d variable
    io.read(tmp2d,"col_dry");
    Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}), KOKKOS_LAMBDA (int ilay, int icol) {
      col_dry(icol,ilay) = tmp2d(0,ilay);
    });
  }

  io.close();
}
#endif

void write_sw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, real2d const &flux_dir, int ncol) {
  yakl::SimpleNetCDF io;
  io.open(fileName , yakl::NETCDF_MODE_WRITE);
  io.write(flux_up  , "sw_flux_up_result"  , {"col_new","lev"});
  io.write(flux_dn  , "sw_flux_dn_result"  , {"col_new","lev"});
  io.write(flux_dir , "sw_flux_dir_result" , {"col_new","lev"});
  io.close();
}

void write_lw_fluxes(std::string fileName, real2d const &flux_up, real2d const &flux_dn, int ncol) {
  yakl::SimpleNetCDF io;
  io.open(fileName , yakl::NETCDF_MODE_WRITE);
  io.write(flux_up , "lw_flux_up_result" , {"col_new","lev"});
  io.write(flux_dn , "lw_flux_dn_result" , {"col_new","lev"});
  io.close();
}

#ifdef RRTMGP_ENABLE_KOKKOS
void write_sw_fluxes(std::string fileName, real2dk const &flux_up, real2dk const &flux_dn, real2dk const &flux_dir, int ncol) {
  conv::SimpleNetCDF io;
  io.open(fileName , conv::NETCDF_MODE_WRITE);
  io.write(flux_up  , "sw_flux_up_result"  , {"col_new","lev"});
  io.write(flux_dn  , "sw_flux_dn_result"  , {"col_new","lev"});
  io.write(flux_dir , "sw_flux_dir_result" , {"col_new","lev"});
  io.close();
}

void write_lw_fluxes(std::string fileName, real2dk const &flux_up, real2dk const &flux_dn, int ncol) {
  conv::SimpleNetCDF io;
  io.open(fileName , conv::NETCDF_MODE_WRITE);
  io.write(flux_up , "lw_flux_up_result" , {"col_new","lev"});
  io.write(flux_dn , "lw_flux_dn_result" , {"col_new","lev"});
  io.close();
}
#endif
