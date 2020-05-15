
#include "mo_optical_props_kernels.h"



// increment 2stream by 2stream
void inc_2stream_by_2stream_bybnd(int ncol, int nlay, int ngpt,
                                  real3d       &tau1, real3d       &ssa1, real3d       &g1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2,
                                  int nbnd, int2d const &gpt_lims) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1 , ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1,ibnd) && igpt <= gpt_lims(2,ibnd) ) {
        // t=tau1 + tau2
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        // w=(tau1*ssa1 + tau2*ssa2) / t
        real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                         tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd);
        g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(icol,ilay,igpt) + 
                              tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd) * g2(icol,ilay,ibnd)) / max(eps,tauscat12);
        ssa1(icol,ilay,igpt) = tauscat12 / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
      }
    }
  });
}



// Incrementing when the second set of optical properties is defined at lower spectral resolution
//   (e.g. by band instead of by gpoint)
void inc_1scalar_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2,
                                  int nbnd, int2d const &gpt_lims) {
  // do igpt = 1 , ngpt
  //   do ilay = 1 , nlay
  //     do icol = 1 , ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1,ibnd) && igpt <= gpt_lims(2,ibnd) ) {
        tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
      }
    }
  });
}



// Delta-scale
//   f = g*g
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, real3d &tau, real3d &ssa, real3d &g) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    if (tau(icol,ilay,igpt) > eps) {
      real f  = g  (icol,ilay,igpt) * g  (icol,ilay,igpt);
      real wf = ssa(icol,ilay,igpt) * f;
      tau(icol,ilay,igpt) = (1._wp - wf) * tau(icol,ilay,igpt);
      ssa(icol,ilay,igpt) = (ssa(icol,ilay,igpt) - wf) / (1.0_wp - wf);
      g  (icol,ilay,igpt) = (g  (icol,ilay,igpt) -  f) / (1.0_wp -  f);
    }
  });
}


// Delta-scaling, provided only for two-stream properties at present
// -------------------------------------------------------------------------------------------------
// Delta-scale
//   user-provided value of f (forward scattering)
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, real3d &tau, real3d &ssa, real3d &g, real3d const &f) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    if (tau(icol,ilay,igpt) > eps) {
      real wf = ssa(icol,ilay,igpt) * f(icol,ilay,igpt);
      tau(icol,ilay,igpt) = (1._wp - wf) * tau(icol,ilay,igpt);
      ssa(icol,ilay,igpt) = (ssa(icol,ilay,igpt) - wf) /  (1.0_wp - wf);
      g  (icol,ilay,igpt) = (g  (icol,ilay,igpt) - f(icol,ilay,igpt)) / (1._wp - f(icol,ilay,igpt));
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// Addition of optical properties: the first set are incremented by the second set.
//
//   There are three possible representations of optical properties (scalar = optical depth only;
//   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
//   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
//   each choice of representation on the left hand side times three representations of the
//   optical properties to be added.
//
//   There are two sets of these nine routines. In the first the two sets of optical
//   properties are defined at the same spectral resolution. There is also a set of routines
//   to add properties defined at lower spectral resolution to a set defined at higher spectral
//   resolution (adding properties defined by band to those defined by g-point)
void increment_1scalar_by_1scalar(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 1scalar by 2stream
void increment_1scalar_by_2stream(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2, real3d const &ssa2) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt) * (1._wp - ssa2(icol,ilay,igpt));
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 1scalar by nstream
void increment_1scalar_by_nstream(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2, real3d const &ssa2) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt) * (1._wp - ssa2(icol,ilay,igpt));
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 2stream by 1scalar
void increment_2stream_by_1scalar(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d const &tau2) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    if (tau12 > eps) {
      ssa1(icol,ilay,igpt) = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) / tau12;
      tau1(icol,ilay,igpt) = tau12;
      // g is unchanged
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 2stream by 2stream
void increment_2stream_by_2stream(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d &g1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    // t=tau1 + tau2
    real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    // w=(tau1*ssa1 + tau2*ssa2) / t
    real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                     tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt);
    if (tauscat12 > eps) {
      g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(icol,ilay,igpt) + 
                            tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt) * g2(icol,ilay,igpt)) 
                             / tauscat12;
      ssa1(icol,ilay,igpt) = tauscat12 / tau12;
      tau1(icol,ilay,igpt) = tau12;
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}


// increment 2stream by nstream
void increment_2stream_by_nstream(int ncol, int nlay, int ngpt, int nmom2, real3d &tau1, real3d &ssa1,
                                  real3d &g1, real3d const &tau2, real3d const &ssa2, real4d const &p2) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    // t=tau1 + tau2
    real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    // w=(tau1*ssa1 + tau2*ssa2) / t
    real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                     tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt);
    if (tauscat12 > eps) {
      g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(   icol,ilay,igpt)+ 
                            tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt) * p2(1, icol,ilay,igpt)) / tauscat12;
      ssa1(icol,ilay,igpt) = tauscat12 / tau12;
      tau1(icol,ilay,igpt) = tau12;
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment nstream by 1scalar
void increment_nstream_by_1scalar(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d const &tau2) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    if (tau12 > eps) {
      ssa1(icol,ilay,igpt) = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) / tau12;
      tau1(icol,ilay,igpt) = tau12;
      // p is unchanged
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment nstream by 2stream
void increment_nstream_by_2stream(int ncol, int nlay, int ngpt, int nmom1, real3d &tau1, real3d &ssa1,
                                  real4d &p1, real3d const &tau2, real3d const &ssa2, real3d const &g2) {
  real1d temp_moms("temp_moms",nmom1);
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                     tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt);
    // Here assume Henyey-Greenstein
    if (tauscat12 > eps) {
      temp_moms(1) = g2(icol,ilay,igpt);
      for (int imom=2; imom<=nmom1; imom++) {
        temp_moms(imom) = temp_moms(imom-1) * g2(icol,ilay,igpt);
      }
      for (int imom=1; imom<=nmom1; imom++) {
        p1(imom, icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * p1(imom, icol,ilay,igpt) + 
                                    tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt) * temp_moms(imom)  ) / tauscat12;
      }
      ssa1(icol,ilay,igpt) = tauscat12 / tau12;
      tau1(icol,ilay,igpt) = tau12;
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment nstream by nstream
void increment_nstream_by_nstream(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real3d &tau1,
                                  real3d &ssa1, real4d &p1, real3d const &tau2, real3d const &ssa2, real4d const &p2) {
  real eps = 3*std::numeric_limits<real>::min();
  int mom_lim = min(nmom1, nmom2);

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                     tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt);
    if (tauscat12 > eps) {
      // If op2 has more moments than op1 these are ignored;
      //   if it has fewer moments the higher orders are assumed to be 0
      for (int imom=1; imom<=mom_lim; imom++) {
        p1(imom, icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * p1(imom, icol,ilay,igpt) + 
                                    tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt) * p2(imom, icol,ilay,igpt)) / max(eps,tauscat12);
      }
      ssa1(icol,ilay,igpt) = tauscat12 / max(eps,tau12);
      tau1(icol,ilay,igpt) = tau12;
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 1scalar by 2stream
void inc_1scalar_by_2stream_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2,
                                  real3d const &ssa2, int nbnd, int2d const &gpt_lims) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd) * (1._wp - ssa2(icol,ilay,ibnd));
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 1scalar by nstream
void inc_1scalar_by_nstream_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d const &tau2,
                                  real3d const &ssa2, int nbnd, int2d const &gpt_lims) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd) * (1._wp - ssa2(icol,ilay,ibnd));
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment 2stream by 1scalar
void inc_2stream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1,
                                  real3d const &tau2, int nbnd, int2d const &gpt_lims) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        ssa1(icol,ilay,igpt) = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
        // g is unchanged
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}


// increment 2stream by nstream
void inc_2stream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom2, real3d &tau1, real3d &ssa1,
                                  real3d &g1, real3d const &tau2, real3d const &ssa2, real4d const &p2,
                                  int nbnd, int2d const &gpt_lims) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        // t=tau1 + tau2
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        // w=(tau1*ssa1 + tau2*ssa2) / t
        real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                         tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd);
        g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(   icol,ilay,igpt)+ 
                              tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd) * p2(1, icol,ilay,ibnd)) / max(eps,tauscat12);
        ssa1(icol,ilay,igpt) = tauscat12 / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment nstream by 1scalar
void inc_nstream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d &tau1, real3d &ssa1, real3d const &tau2,
                                  int nbnd, int2d const &gpt_lims) {
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        ssa1(icol,ilay,igpt) = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
        // p is unchanged
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment nstream by 2stream
void inc_nstream_by_2stream_bybnd(int ncol, int nlay, int ngpt, int nmom1, real3d &tau1, real3d &ssa1, real4d &p1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2, int nbnd,
                                  int2d const &gpt_lims) {
  real1d temp_moms("temp_moms",nmom1);
  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) +
                         tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd);
        // Here assume Henyey-Greenstein
        temp_moms(1) = g2(icol,ilay,ibnd);
        for (int imom=2; imom<=nmom1; imom++) {
          temp_moms(imom) = temp_moms(imom-1) * g2(icol,ilay,ibnd);
        }
        for (int imom=1; imom<=nmom1; imom++) {
          p1(imom, icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * p1(imom, icol,ilay,igpt) + 
                                      tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd) * temp_moms(imom)  ) / max(eps,tauscat12);
        }
        ssa1(icol,ilay,igpt) = tauscat12 / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// increment nstream by nstream
void inc_nstream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real3d &tau1, real3d &ssa1, real4d &p1,
                                  real3d const &tau2, real3d const &ssa2, real4d const &p2, int nbnd, int2d const &gpt_lims) {
  real eps = 3*std::numeric_limits<real>::min();
  int mom_lim = min(nmom1, nmom2);

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=0; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1, ibnd) && igpt <= gpt_lims(2, ibnd) ) {
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                         tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd);
        // If op2 has more moments than op1 these are ignored;
        //   if it has fewer moments the higher orders are assumed to be 0
        for (int imom=1; imom<=mom_lim; imom++) {
          p1(imom, icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * p1(imom, icol,ilay,igpt) + 
                                      tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd) * p2(imom, icol,ilay,ibnd)) / max(eps,tauscat12);
        }
        ssa1(icol,ilay,igpt) = tauscat12 / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
      }
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// Subsetting, meaning extracting some portion of the 3D domain
void extract_subset_dim1_3d(int ncol, int nlay, int ngpt, real3d const &array_in, int colS, int colE, real3d &array_out) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = colS, colE
  parallel_for( Bounds<3>(ngpt,nlay,{colS,colE}) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    array_out(icol-colS+1, ilay, igpt) = array_in(icol, ilay, igpt);
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



void extract_subset_dim2_4d(int nmom, int ncol, int nlay, int ngpt, real3d const &array_in, int colS, int colE, real3d &array_out) {
  //do igpt = 1, ngpt
  //  do ilay = 1, nlay
  //    do icol = colS, colE
  //      do imom = 1, nmom
  parallel_for( Bounds<4>(ngpt,nlay,{colS,colE},nmom) , YAKL_LAMBDA (int igpt, int ilay, int icol, int imom) {
    array_out(imom, icol-colS+1, ilay, igpt) = array_in(imom, icol, ilay, igpt);
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}


// Extract the absorption optical thickness which requires mulitplying by 1 - ssa
void extract_subset_absorption_tau(int ncol, int nlay, int ngpt, real3d const &tau_in, real3d const &ssa_in, int colS, int colE,
                                   real3d &tau_out) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = colS, colE
  parallel_for( Bounds<3>(ngpt,nlay,{colS,colE}) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    tau_out(icol-colS+1, ilay, igpt) = tau_in(icol, ilay, igpt) * (1._wp - ssa_in(icol, ilay, igpt));
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}




