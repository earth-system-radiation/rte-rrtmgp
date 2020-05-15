
#include "mo_rte_sw.h"

void rte_sw(OpticalProps2str const &atmos, bool top_at_1, real1d const &mu0, real2d const &inc_flux,
            real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, FluxesBroadband &fluxes, real2d const &inc_flux_dif) {
  real3d gpt_flux_up;
  real3d gpt_flux_dn;
  real3d gpt_flux_dir;
  real2d sfc_alb_dir_gpt;
  real2d sfc_alb_dif_gpt;
  int ncol  = atmos.get_ncol();
  int nlay  = atmos.get_nlay();
  int ngpt  = atmos.get_ngpt();
  int nband = atmos.get_nband();

  // Error checking -- consistency of sizes and validity of values
  if (! fluxes.are_desired()) { stoprun("rte_sw: no space allocated for fluxes"); }

  if (size(mu0,1) != ncol) { stoprun("rte_sw: mu0 inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (anyLT(mu0,0._wp) || anyGT(mu0,1._wp)) { stoprun("rte_sw: one or more mu0 <= 0 or > 1"); }
  #endif

  if (size(inc_flux,1) != ncol || size(inc_flux,2) != ngpt) { stoprun("rte_sw: inc_flux inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (anyLT(inc_flux,0._wp)) { stoprun("rte_sw: one or more inc_flux < 0"); }
  #endif
  if (allocated(inc_flux_dif)) {
    if (size(inc_flux_dif,1) != ncol || size(inc_flux_dif,2) != ngpt) { stoprun("rte_sw: inc_flux_dif inconsistently sized"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (anyLT(inc_flux_dif,0._wp)) { stoprun("rte_sw: one or more inc_flux_dif < 0"); }
    #endif
  }

  if (size(sfc_alb_dir,1) != nband || size(sfc_alb_dir,2) != ncol) { stoprun("rte_sw: sfc_alb_dir inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (anyLT(sfc_alb_dir,0._wp) || anyGT(sfc_alb_dir,1._wp)) { stoprun("rte_sw: sfc_alb_dir out of bounds [0,1]"); }
  #endif
  if (size(sfc_alb_dif,1) != nband || size(sfc_alb_dif,2) != ncol) { stoprun("rte_sw: sfc_alb_dif inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (anyLT(sfc_alb_dif,0._wp) || anyGT(sfc_alb_dif,1._wp)) { stoprun("rte_sw: sfc_alb_dif out of bounds [0,1]"); }
  #endif

  gpt_flux_up  = real3d("gpt_flux_up" ,ncol, nlay+1, ngpt);
  gpt_flux_dn  = real3d("gpt_flux_dn" ,ncol, nlay+1, ngpt);
  gpt_flux_dir = real3d("gpt_flux_dir",ncol, nlay+1, ngpt);
  sfc_alb_dir_gpt = real2d("sfc_alb_dir_gpt",ncol, ngpt);
  sfc_alb_dif_gpt = real2d("sfc_alb_dif_gpt",ncol, ngpt);
  // Lower boundary condition -- expand surface albedos by band to gpoints
  //   and switch dimension ordering
  expand_and_transpose(atmos, sfc_alb_dir, sfc_alb_dir_gpt);
  expand_and_transpose(atmos, sfc_alb_dif, sfc_alb_dif_gpt);

  // Compute the radiative transfer...
  // Apply boundary conditions
  //   On input flux_dn is the diffuse component; the last action in each solver is to add
  //   direct and diffuse to represent the total, consistent with the LW
  apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, mu0, gpt_flux_dir);
  if (allocated(inc_flux_dif)) {
    apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif,  gpt_flux_dn );
  } else {
    apply_BC(ncol, nlay, ngpt, top_at_1,                gpt_flux_dn );
  }

  // two-stream calculation with scattering
  atmos.validate();
  sw_solver_2stream(ncol, nlay, ngpt, top_at_1, 
                    atmos.tau, atmos.ssa, atmos.g, mu0,      
                    sfc_alb_dir_gpt, sfc_alb_dif_gpt,        
                    gpt_flux_up, gpt_flux_dn, gpt_flux_dir);

  // ...and reduce spectral fluxes to desired output quantities
  fluxes.reduce(gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir);
}


