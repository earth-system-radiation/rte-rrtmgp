
#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props.h"
#include "mo_cloud_optics.h"

#ifdef RRTMGP_ENABLE_YAKL
void load_cld_lutcoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);

void load_cld_padecoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// read cloud optical property LUT coefficients from NetCDF file
template <typename RealT, typename LayoutT, typename DeviceT>
void load_cld_lutcoeff(CloudOpticsK<RealT, LayoutT, DeviceT> &cloud_spec,
                       const std::string& cld_coeff_file) {
  using real2dk = Kokkos::View<RealT**,   LayoutT, DeviceT>;
  using real3dk = Kokkos::View<RealT***,  LayoutT, DeviceT>;

  conv::SimpleNetCDF io;
  // Open cloud optical property coefficient file
  io.open(cld_coeff_file , conv::NETCDF_MODE_READ);

  // Read LUT coefficient dimensions
  int nband     = io.getDimSize("nband");
  int nrghice   = io.getDimSize("nrghice");
  int nsize_liq = io.getDimSize("nsize_liq");
  int nsize_ice = io.getDimSize("nsize_ice");

  real2dk band_lims_wvn("band_lims_wvn",2,nband);
  io.read(band_lims_wvn,"bnd_limits_wavenumber");

  // Read LUT constants
  RealT radliq_lwr, radliq_upr, radliq_fac, radice_lwr, radice_upr, radice_fac;
  io.read(radliq_lwr , "radliq_lwr");
  io.read(radliq_upr , "radliq_upr");
  io.read(radliq_fac , "radliq_fac");
  io.read(radice_lwr , "radice_lwr");
  io.read(radice_upr , "radice_upr");
  io.read(radice_fac , "radice_fac");

  // Allocate cloud property lookup table input arrays
  real2dk lut_extliq("lut_extliq",nsize_liq, nband);
  real2dk lut_ssaliq("lut_ssaliq",nsize_liq, nband);
  real2dk lut_asyliq("lut_asyliq",nsize_liq, nband);
  real3dk lut_extice("lut_extice",nsize_ice, nband, nrghice);
  real3dk lut_ssaice("lut_ssaice",nsize_ice, nband, nrghice);
  real3dk lut_asyice("lut_asyice",nsize_ice, nband, nrghice);
  // Read LUT coefficients
  io.read(lut_extliq , "lut_extliq");
  io.read(lut_ssaliq , "lut_ssaliq");
  io.read(lut_asyliq , "lut_asyliq");
  io.read(lut_extice , "lut_extice");
  io.read(lut_ssaice , "lut_ssaice");
  io.read(lut_asyice , "lut_asyice");

  io.close();

  cloud_spec.load(band_lims_wvn, radliq_lwr, radliq_upr, radliq_fac, radice_lwr, radice_upr, radice_fac,
                                 lut_extliq, lut_ssaliq, lut_asyliq, lut_extice, lut_ssaice, lut_asyice);
}

// read cloud optical property Pade coefficients from NetCDF file
template <typename RealT, typename LayoutT, typename DeviceT>
void load_cld_padecoeff(CloudOpticsK<RealT, LayoutT, DeviceT> &cloud_spec,
                        const std::string& cld_coeff_file) {
  using real1dk = Kokkos::View<RealT*,     LayoutT, DeviceT>;
  using real2dk = Kokkos::View<RealT**,    LayoutT, DeviceT>;
  using real3dk = Kokkos::View<RealT***,   LayoutT, DeviceT>;
  using real4dk = Kokkos::View<RealT****,  LayoutT, DeviceT>;

  conv::SimpleNetCDF io;
  // Open cloud optical property coefficient file
  io.open(cld_coeff_file , conv::NETCDF_MODE_READ);

  // Read Pade coefficient dimensions
  int nband        = io.getDimSize("nband");
  int nrghice      = io.getDimSize("nrghice");
  int nsizereg     = io.getDimSize("nsizereg");
  int ncoeff_ext   = io.getDimSize("ncoeff_ext");
  int ncoeff_ssa_g = io.getDimSize("ncoeff_ssa_g");
  int nbound       = io.getDimSize("nbound");

  real2dk band_lims_wvn("band_lims_wvn",2,nband);
  io.read(band_lims_wvn, "bnd_limits_wavenumber");

  // Allocate cloud property Pade coefficient input arrays
  real3dk pade_extliq("pade_extliq",nband, nsizereg, ncoeff_ext);
  real3dk pade_ssaliq("pade_ssaliq",nband, nsizereg, ncoeff_ssa_g);
  real3dk pade_asyliq("pade_asyliq",nband, nsizereg, ncoeff_ssa_g);
  real4dk pade_extice("pade_extice",nband, nsizereg, ncoeff_ext,   nrghice);
  real4dk pade_ssaice("pade_ssaice",nband, nsizereg, ncoeff_ssa_g, nrghice);
  real4dk pade_asyice("pade_asyice",nband, nsizereg, ncoeff_ssa_g, nrghice);
  io.read(pade_extliq, "pade_extliq");
  io.read(pade_ssaliq, "pade_ssaliq");
  io.read(pade_asyliq, "pade_asyliq");
  io.read(pade_extice, "pade_extice");
  io.read(pade_ssaice, "pade_ssaice");
  io.read(pade_asyice, "pade_asyice");

  // Allocate cloud property Pade coefficient particle size boundary input arrays
  real1dk pade_sizreg_extliq("pade_sizreg_extliq",nbound);
  real1dk pade_sizreg_ssaliq("pade_sizreg_ssaliq",nbound);
  real1dk pade_sizreg_asyliq("pade_sizreg_asyliq",nbound);
  real1dk pade_sizreg_extice("pade_sizreg_extice",nbound);
  real1dk pade_sizreg_ssaice("pade_sizreg_ssaice",nbound);
  real1dk pade_sizreg_asyice("pade_sizreg_asyice",nbound);

  io.read(pade_sizreg_extliq, "pade_sizreg_extliq");
  io.read(pade_sizreg_ssaliq, "pade_sizreg_ssaliq");
  io.read(pade_sizreg_asyliq, "pade_sizreg_asyliq");
  io.read(pade_sizreg_extice, "pade_sizreg_extice");
  io.read(pade_sizreg_ssaice, "pade_sizreg_ssaice");
  io.read(pade_sizreg_asyice, "pade_sizreg_asyice");

  io.close();

  cloud_spec.load(band_lims_wvn, pade_extliq, pade_ssaliq, pade_asyliq,
                                 pade_extice, pade_ssaice, pade_asyice,
                                 pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq,
                                 pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice);
}
#endif
