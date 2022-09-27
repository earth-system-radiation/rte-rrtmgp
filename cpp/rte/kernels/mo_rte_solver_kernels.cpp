
#include "mo_rte_solver_kernels.h"



void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &flux_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  //   Upper boundary condition
  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      flux_dn(icol,      1, igpt)  = 0;
    });
  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      flux_dn(icol, nlay+1, igpt)  = 0;
    });
  }
}



void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &inc_flux, real1d const &factor, real3d const &flux_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      flux_dn(icol,      1, igpt)  = inc_flux(icol,igpt) * factor(icol);
    });
  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      flux_dn(icol, nlay+1, igpt)  = inc_flux(icol,igpt) * factor(icol);
    });
  }
}



// Upper boundary condition
void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &inc_flux, real3d const &flux_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  //   Upper boundary condition
  if (top_at_1) {
    //$acc  parallel loop collapse(2)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      flux_dn(icol,      1, igpt)  = inc_flux(icol,igpt);
    });
  } else {
    //$acc  parallel loop collapse(2)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      flux_dn(icol, nlay+1, igpt)  = inc_flux(icol,igpt);
    });
  }
}



// Transport of diffuse radiation through a vertically layered atmosphere.
//   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
//   This routine is shared by longwave and shortwave
void adding(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &albedo_sfc, real3d const &rdif, real3d const &tdif,
            real3d const &src_dn, real3d const &src_up, real2d const &src_sfc, real3d const &flux_up, real3d const &flux_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real3d albedo("albedo",ncol,nlay+1,ngpt);
  real3d src   ("src   ",ncol,nlay+1,ngpt);
  real3d denom ("denom ",ncol,nlay  ,ngpt);

  // Indexing into arrays for upward and downward propagation depends on the vertical
  //   orientation of the arrays (whether the domain top is at the first or last index)
  // We write the loops out explicitly so compilers will have no trouble optimizing them.
  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      int ilev = nlay + 1;
      // Albedo of lowest level is the surface albedo...
      albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
      // ... and source of diffuse radiation is surface emission
      src(icol,ilev,igpt) = src_sfc(icol,igpt);

      // From bottom to top of atmosphere --
      //   compute albedo and source of upward radiation
      for (ilev=nlay; ilev>=1; ilev--) {
        denom(icol,ilev,igpt) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(icol,ilev+1,igpt));    // Eq 10
        albedo(icol,ilev,igpt) = rdif(icol,ilev,igpt) + 
                                 tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev+1,igpt) * denom(icol,ilev,igpt); // Equation 9
        // Equation 11 -- source is emitted upward radiation at top of layer plus
        //   radiation emitted at bottom of layer,
        //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
        src(icol,ilev,igpt) =  src_up(icol, ilev, igpt) + 
                               tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *       
                               (src(icol,ilev+1,igpt) + albedo(icol,ilev+1,igpt)*src_dn(icol,ilev,igpt));
      }

      // Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = 1;
      flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
                                src(icol,ilev,igpt);                                  // emission from below

      // From the top of the atmosphere downward -- compute fluxes
      for (ilev = 2; ilev <= nlay+1; ilev++) {
        flux_dn(icol,ilev,igpt) = (tdif(icol,ilev-1,igpt)*flux_dn(icol,ilev-1,igpt) +   // Equation 13
                                  rdif(icol,ilev-1,igpt)*src(icol,ilev,igpt) +       
                                  src_dn(icol,ilev-1,igpt)) * denom(icol,ilev-1,igpt);
        flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                  src(icol,ilev,igpt);
      }
    });

  } else {

    #ifdef RRTMGP_CPU_KERNELS
      #ifdef YAKL_AUTO_PROFILE
        auto timername = std::string(YAKL_AUTO_LABEL());
        yakl::timer_start(timername.c_str());
      #endif
      #ifdef YAKL_ARCH_OPENMP
        #pragma omp parallel for
      #endif
      for (int igpt = 1; igpt <= ngpt; igpt++) {
        int ilev = 1;
        for (int icol = 1; icol <= ncol; icol++) {
          // Albedo of lowest level is the surface albedo...
          albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
          // ... and source of diffuse radiation is surface emission
          src(icol,ilev,igpt) = src_sfc(icol,igpt);
        }

          // From bottom to top of atmosphere --
          //   compute albedo and source of upward radiation
        for (ilev = 1; ilev <= nlay; ilev++) {
          for (int icol = 1; icol <= ncol; icol++) {
            denom (icol,ilev  ,igpt) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(icol,ilev,igpt));                // Eq 10
            albedo(icol,ilev+1,igpt) = rdif(icol,ilev,igpt) + 
                                       tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev,igpt) * denom(icol,ilev,igpt); // Equation 9
            // Equation 11 -- source is emitted upward radiation at top of layer plus
            //   radiation emitted at bottom of layer,
            //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
            src(icol,ilev+1,igpt) =  src_up(icol, ilev, igpt) +  
                                     tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *       
                                     (src(icol,ilev,igpt) + albedo(icol,ilev,igpt)*src_dn(icol,ilev,igpt));
          }
        }

        // Eq 12, at the top of the domain upwelling diffuse is due to ...
        ilev = nlay+1;
        for (int icol = 1; icol <= ncol; icol++) {
          flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
                                    src(icol,ilev,igpt);                          // scattering by the direct beam below
        }

        // From the top of the atmosphere downward -- compute fluxes
        for (ilev=nlay; ilev >= 1; ilev--) {
          for (int icol = 1; icol <= ncol; icol++) {
            flux_dn(icol,ilev,igpt) = (tdif(icol,ilev,igpt)*flux_dn(icol,ilev+1,igpt) +   // Equation 13
                                      rdif(icol,ilev,igpt)*src(icol,ilev,igpt) + 
                                      src_dn(icol, ilev, igpt)) * denom(icol,ilev,igpt);
            flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                      src(icol,ilev,igpt);

          }
        }
      }
      #ifdef YAKL_AUTO_PROFILE
        yakl::timer_stop(timername.c_str());
      #endif
    #else
      // do igpt = 1, ngpt
      //   do icol = 1, ncol
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
        int ilev = 1;
        // Albedo of lowest level is the surface albedo...
        albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
        // ... and source of diffuse radiation is surface emission
        src(icol,ilev,igpt) = src_sfc(icol,igpt);

        // From bottom to top of atmosphere --
        //   compute albedo and source of upward radiation
        for (ilev = 1; ilev <= nlay; ilev++) {
          denom (icol,ilev  ,igpt) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(icol,ilev,igpt));                // Eq 10
          albedo(icol,ilev+1,igpt) = rdif(icol,ilev,igpt) + 
                                     tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev,igpt) * denom(icol,ilev,igpt); // Equation 9
          // Equation 11 -- source is emitted upward radiation at top of layer plus
          //   radiation emitted at bottom of layer,
          //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
          src(icol,ilev+1,igpt) =  src_up(icol, ilev, igpt) +  
                                   tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *       
                                   (src(icol,ilev,igpt) + albedo(icol,ilev,igpt)*src_dn(icol,ilev,igpt));
        }

        // Eq 12, at the top of the domain upwelling diffuse is due to ...
        ilev = nlay+1;
        flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
                                  src(icol,ilev,igpt);                          // scattering by the direct beam below

        // From the top of the atmosphere downward -- compute fluxes
        for (ilev=nlay; ilev >= 1; ilev--) {
          flux_dn(icol,ilev,igpt) = (tdif(icol,ilev,igpt)*flux_dn(icol,ilev+1,igpt) +   // Equation 13
                                    rdif(icol,ilev,igpt)*src(icol,ilev,igpt) + 
                                    src_dn(icol, ilev, igpt)) * denom(icol,ilev,igpt);
          flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                    src(icol,ilev,igpt);

        }
      });
    #endif
  }
}



void sw_solver_2stream(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real3d const &ssa, real3d const &g,
                       real1d const &mu0, real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, real3d const &flux_up,
                       real3d const &flux_dn, real3d const &flux_dir) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real3d Rdif      ("Rdif      ",ncol,nlay,ngpt);         
  real3d Tdif      ("Tdif      ",ncol,nlay,ngpt);         
  real3d Rdir      ("Rdir      ",ncol,nlay,ngpt);         
  real3d Tdir      ("Tdir      ",ncol,nlay,ngpt);         
  real3d Tnoscat   ("Tnoscat   ",ncol,nlay,ngpt);         
  real3d source_up ("source_up ",ncol,nlay,ngpt);         
  real3d source_dn ("source_dn ",ncol,nlay,ngpt);         
  real2d source_srf("source_srf",ncol     ,ngpt);         

  // Cell properties: transmittance and reflectance for direct and diffuse radiation
  sw_two_stream(ncol, nlay, ngpt, mu0, 
                tau , ssa , g   ,      
                Rdif, Tdif, Rdir, Tdir, Tnoscat);

  sw_source_2str(ncol, nlay, ngpt, top_at_1,       
                 Rdir, Tdir, Tnoscat, sfc_alb_dir, 
                 source_up, source_dn, source_srf, flux_dir);

  adding(ncol, nlay, ngpt, top_at_1,   
         sfc_alb_dif, Rdif, Tdif,      
         source_dn, source_up, source_srf, flux_up, flux_dn);

  // adding computes only diffuse flux; flux_dn is total
  //
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay+1
  //     do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay+1,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    flux_dn(icol,ilay,igpt) = flux_dn(icol,ilay,igpt) + flux_dir(icol,ilay,igpt);
  });
}



// Top-level longwave kernels
//
// LW fluxes, no scattering, mu (cosine of integration angle) specified by column
//   Does radiation calculation at user-supplied angles; converts radiances to flux
//   using user-supplied weights
void lw_solver_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &D, real1d const &weights, int weight_ind, real3d const &tau,
                      real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                      real2d const &sfc_emis, real2d const &sfc_src, real3d const &radn_up, real3d const &radn_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real3d tau_loc   ("tau_loc   ",ncol,nlay,ngpt);             
  real3d trans     ("trans     ",ncol,nlay,ngpt);             
  real2d source_sfc("source_sfc",ncol,     ngpt);             
  real2d sfc_albedo("sfc_albedo",ncol,     ngpt);             

  real tau_thresh = sqrt( std::numeric_limits<real>::epsilon() );

  real constexpr pi = M_PI;

  // Which way is up?
  // Level Planck sources for upward and downward radiation
  // When top_at_1, lev_source_up => lev_source_dec
  //                lev_source_dn => lev_source_inc, and vice-versa
  int top_level;
  real3d lev_source_up;
  real3d lev_source_dn;
  if (top_at_1) {
    top_level = 1;
    // Recall below that equating two arrays is like assigning pointers in Fortran. No data is copied.
    // The LHS just uses the same data pointer as the RHS so that changing one's data changes the other's as well.
    lev_source_up = lev_source_dec;
    lev_source_dn = lev_source_inc;
  } else {
    top_level = nlay+1;
    // Recall below that equating two arrays is like assigning pointers in Fortran. No data is copied.
    // The LHS just uses the same data pointer as the RHS so that changing one's data changes the other's as well.
    lev_source_up = lev_source_inc;
    lev_source_dn = lev_source_dec;
  }

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
    // Transport is for intensity
    //   convert flux at top of domain to intensity assuming azimuthal isotropy
    radn_dn(icol,top_level,igpt) = radn_dn(icol,top_level,igpt)/(2._wp * pi * weights(weight_ind));
    
    // Surface albedo, surface source function
    sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt);
    source_sfc(icol,igpt) = sfc_emis(icol,igpt) * sfc_src(icol,igpt);
  });

  real3d source_dn ("source_dn ",ncol,nlay,ngpt);             
  real3d source_up ("source_up ",ncol,nlay,ngpt);             
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    // Optical path and transmission, used in source function and transport calculations
    tau_loc(icol,ilay,igpt) = tau(icol,ilay,igpt)*D(icol,igpt);
    trans  (icol,ilay,igpt) = exp(-tau_loc(icol,ilay,igpt));

    lw_source_noscat_stencil(ncol, nlay, ngpt, icol, ilay, igpt,        
                             lay_source, lev_source_up, lev_source_dn,  
                             tau_loc, trans,                            
                             source_dn, source_up, tau_thresh);
  });

  // Transport
  lw_transport_noscat(ncol, nlay, ngpt, top_at_1,  
                      tau_loc, trans, sfc_albedo, source_dn, source_up, source_sfc, 
                      radn_up, radn_dn);

  // Convert intensity to flux assuming azimuthal isotropy and quadrature weight
  // do igpt = 1, ngpt
  //   do ilev = 1, nlay+1
  //     do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay+1,ncol) , YAKL_LAMBDA (int igpt, int ilev, int icol) {
    radn_dn(icol,ilev,igpt) = 2._wp * pi * weights(weight_ind) * radn_dn(icol,ilev,igpt);
    radn_up(icol,ilev,igpt) = 2._wp * pi * weights(weight_ind) * radn_up(icol,ilev,igpt);
  });
}



// LW transport, no scattering, multi-angle quadrature
//   Users provide a set of weights and quadrature angles
//   Routine sums over single-angle solutions for each sets of angles/weights
void lw_solver_noscat_GaussQuad(int ncol, int nlay, int ngpt, bool top_at_1, int nmus, real1d const &Ds, real1d const &weights, 
                                real3d const &tau, real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                                real2d const &sfc_emis, real2d const &sfc_src, real3d const &flux_up, real3d const &flux_dn) {
  using yakl::intrinsics::merge;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  // Local variables
  real3d radn_dn ("radn_dn ",ncol,nlay+1,ngpt);  
  real3d radn_up ("radn_up ",ncol,nlay+1,ngpt);  
  real2d Ds_ncol ("Ds_ncol ",ncol,       ngpt);  
  real2d flux_top("flux_top",ncol,       ngpt);  

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
    Ds_ncol(icol, igpt) = Ds(1);
  });

  lw_solver_noscat(ncol, nlay, ngpt, 
                   top_at_1, Ds_ncol, weights, 1, tau, 
                   lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, 
                   flux_up, flux_dn);
  //
  // For more than one angle use local arrays
  int top_level = merge(1, nlay+1, top_at_1);

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
    flux_top(icol,igpt) = flux_dn(icol,top_level,igpt);
  });

  apply_BC(ncol, nlay, ngpt, top_at_1, flux_top, radn_dn);

  for (int imu=2; imu<=nmus; imu++) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      Ds_ncol(icol, igpt) = Ds(imu);
    });

    lw_solver_noscat(ncol, nlay, ngpt, 
                     top_at_1, Ds_ncol, weights, imu, tau, 
                     lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, 
                     radn_up, radn_dn);

    // do igpt = 1, ngpt
    //   do ilev = 1, nlay+1
    //     do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay+1,ncol) , YAKL_LAMBDA (int igpt, int ilev, int icol) {
      flux_up(icol,ilev,ngpt) = flux_up(icol,ilev,ngpt) + radn_up(icol,ilev,ngpt);
      flux_dn(icol,ilev,ngpt) = flux_dn(icol,ilev,ngpt) + radn_dn(icol,ilev,ngpt);
    });

  } // imu
}



// Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
//   This version straight from ECRAD
//   Source is provided as W/m2-str; factor of pi converts to flux units
void lw_source_2str(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &sfc_emis, real2d const &sfc_src,
                    real3d const &lay_source, real3d const &lev_source, real3d const &gamma1, real3d const &gamma2,
                    real3d const &rdif, real3d const &tdif, real3d const &tau, real3d const &source_dn, real3d const &source_up,
                    real2d const &source_sfc) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real constexpr pi = 3.14159265358979323846;

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    if ( tau(icol,ilay,ngpt) > 1.0e-8_wp ) {
      real lev_source_top, lev_source_bot;
      if (top_at_1) {
        lev_source_top = lev_source(icol,ilay  ,ngpt);
        lev_source_bot = lev_source(icol,ilay+1,ngpt);
      } else {
        lev_source_top = lev_source(icol,ilay+1,ngpt);
        lev_source_bot = lev_source(icol,ilay  ,ngpt);
      }
      // Toon et al. (JGR 1989) Eqs 26-27
      real Z = (lev_source_bot-lev_source_top) / (tau(icol,ilay,igpt)*(gamma1(icol,ilay,igpt)+gamma2(icol,ilay,igpt)));
      real Zup_top        =  Z + lev_source_top;
      real Zup_bottom     =  Z + lev_source_bot;
      real Zdn_top        = -Z + lev_source_top;
      real Zdn_bottom     = -Z + lev_source_bot;
      source_up(icol,ilay,igpt) = pi * (Zup_top    - rdif(icol,ilay,igpt) * Zdn_top    - tdif(icol,ilay,igpt) * Zup_bottom);
      source_dn(icol,ilay,igpt) = pi * (Zdn_bottom - rdif(icol,ilay,igpt) * Zup_bottom - tdif(icol,ilay,igpt) * Zdn_top);
    } else {
      source_up(icol,ilay,igpt) = 0._wp;
      source_dn(icol,ilay,igpt) = 0._wp;
    }
    if(ilay == 1) {
      source_sfc(icol,igpt) = pi * sfc_emis(icol,igpt) * sfc_src(icol,igpt);
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// Source function combination
// RRTMGP provides two source functions at each level
//   using the spectral mapping from each of the adjascent layers.
//   Need to combine these for use in two-stream calculation.
void lw_combine_sources(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &lev_src_inc, real3d const &lev_src_dec,
                        real3d const &lev_source) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay+1
  //     do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay+1,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    if (ilay == 1) {
      lev_source(icol, ilay, igpt) =      lev_src_dec(icol, ilay,   igpt);
    } else if (ilay == nlay+1) {
      lev_source(icol, ilay, igpt) =      lev_src_inc(icol, ilay-1, igpt);
    } else {
      lev_source(icol, ilay, igpt) = sqrt(lev_src_dec(icol, ilay, igpt) * 
                                          lev_src_inc(icol, ilay-1, igpt));
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
//    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
// Equations are developed in Meador and Weaver, 1980,
//    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
void lw_two_stream(int ncol, int nlay, int ngpt, real3d const &tau, real3d const &w0, real3d const &g, real3d const &gamma1,
                   real3d const &gamma2, real3d const &Rdif, real3d const &Tdif) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real constexpr LW_diff_sec = 1.66;  // 1./cos(diffusivity angle)

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    // Coefficients differ from SW implementation because the phase function is more isotropic
    //   Here we follow Fu et al. 1997, doi:10.1175/1520-0469(1997)054<2799:MSPITI>2.0.CO;2
    //   and use a diffusivity sec of 1.66
    gamma1(icol,ilay,igpt)= LW_diff_sec * (1._wp - 0.5_wp * w0(icol,ilay,igpt) * (1._wp + g(icol,ilay,igpt))); // Fu et al. Eq 2.9
    gamma2(icol,ilay,igpt)= LW_diff_sec *          0.5_wp * w0(icol,ilay,igpt) * (1._wp - g(icol,ilay,igpt));  // Fu et al. Eq 2.10

    // Written to encourage vectorization of exponential, square root
    // Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
    //   k = 0 for isotropic, conservative scattering; this lower limit on k
    //   gives relative error with respect to conservative solution
    //   of < 0.1% in Rdif down to tau = 10^-9
    real k = sqrt(std::max((gamma1(icol,ilay,igpt) - gamma2(icol,ilay,igpt)) * 
                           (gamma1(icol,ilay,igpt) + gamma2(icol,ilay,igpt)) , 1.e-12_wp));
    real exp_minusktau = exp(-tau(icol,ilay,igpt)*k);

    // Diffuse reflection and transmission
    real exp_minus2ktau = exp_minusktau * exp_minusktau;

    // Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    real RT_term = 1._wp / (k * (1._wp + exp_minus2ktau)  +  gamma1(icol,ilay,igpt) * (1._wp - exp_minus2ktau) );

    // Equation 25
    Rdif(icol,ilay,igpt) = RT_term * gamma2(icol,ilay,igpt) * (1._wp - exp_minus2ktau);

    // Equation 26
    Tdif(icol,ilay,igpt) = RT_term * 2._wp * k * exp_minusktau;
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



//   Top-level shortwave kernels
// -------------------------------------------------------------------------------------------------
//   Extinction-only i.e. solar direct beam
void sw_solver_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real1d const &mu0, real3d const &flux_dir) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real1d mu0_inv("mu0_inv",ncol);

  parallel_for( YAKL_AUTO_LABEL() , ncol , YAKL_LAMBDA (int icol) {
    mu0_inv(icol) = 1._wp/mu0(icol);
  });

  // Indexing into arrays for upward and downward propagation depends on the vertical
  //   orientation of the arrays (whether the domain top is at the first or last index)
  // We write the loops out explicitly so compilers will have no trouble optimizing them.
  // Downward propagation
  if (top_at_1) {
    // For the flux at this level, what was the previous level, and which layer has the
    //   radiation just passed through?
    // layer index = level index - 1
    // previous level is up (-1)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      for (int ilev=2; ilev<=nlay+1; ilev++) {
        flux_dir(icol,ilev,igpt) = flux_dir(icol,ilev-1,igpt) * exp(-tau(icol,ilev,igpt)*mu0_inv(icol));
      }
    });
  } else {
    // layer index = level index
    // previous level is up (+1)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      for (int ilev=nlay; ilev>=1; ilev--) {
        flux_dir(icol,ilev,igpt) = flux_dir(icol,ilev+1,igpt) * exp(-tau(icol,ilev,igpt)*mu0_inv(icol));
      }
    });
  }
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// Longwave two-stream calculation:
//   combine RRTMGP-specific sources at levels
//   compute layer reflectance, transmittance
//   compute total source function at levels using linear-in-tau
//   transport
void lw_solver_2stream(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real3d const &ssa, real3d const &g,
                       real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                       real2d const &sfc_emis, real2d const &sfc_src, real3d const &flux_up, real3d const &flux_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real3d Rdif      ("Rdif      ",ncol,nlay  ,ngpt);    
  real3d Tdif      ("Tdif      ",ncol,nlay  ,ngpt);    
  real3d gamma1    ("gamma1    ",ncol,nlay  ,ngpt);      
  real3d gamma2    ("gamma2    ",ncol,nlay  ,ngpt);      
  real2d sfc_albedo("sfc_albedo",ncol       ,ngpt);          
  real3d lev_source("lev_source",ncol,nlay+1,ngpt);          
  real3d source_dn ("source_dn ",ncol,nlay  ,ngpt);         
  real3d source_up ("source_up ",ncol,nlay  ,ngpt);         
  real2d source_sfc("source_sfc",ncol       ,ngpt);          

  // RRTMGP provides source functions at each level using the spectral mapping
  //   of each adjacent layer. Combine these for two-stream calculations
  lw_combine_sources(ncol, nlay, ngpt, top_at_1, 
                     lev_source_inc, lev_source_dec, 
                     lev_source);

  // Cell properties: reflection, transmission for diffuse radiation
  //   Coupling coefficients needed for source function
  lw_two_stream(ncol, nlay, ngpt, 
                tau , ssa, g,     
                gamma1, gamma2, Rdif, Tdif);

  // Source function for diffuse radiation
  lw_source_2str(ncol, nlay, ngpt, top_at_1, 
                 sfc_emis, sfc_src, 
                 lay_source, lev_source, 
                 gamma1, gamma2, Rdif, Tdif, tau, 
                 source_dn, source_up, source_sfc);

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
    sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt);
  });

  // Transport
  adding(ncol, nlay, ngpt, top_at_1,        
         sfc_albedo,                        
         Rdif, Tdif,                        
         source_dn, source_up, source_sfc,  
         flux_up, flux_dn);

  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}





