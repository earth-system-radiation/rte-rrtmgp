
#pragma once

#include "rrtmgp_const.h"

template <typename RealT>
struct rrtmgp_constants
{

// -----------------------------------------
// Physical constants, 2018 SI defintion of metric system
//   doi:10.1088/1681-7575/aa950a (see also https://www.nist.gov/si-redefinition/meet-constants)
// Boltzmann constant [J/K] = [(kg m^2)/(K s^2)]
static RealT constexpr k_boltz = 1.380649e-23;

//  molecular weight of water [kg/mol]
static RealT constexpr m_h2o =  0.018016;

// Avogadro's number [molec/mol]
static RealT constexpr avogad = 6.02214076e23;

// Universal gas constant [J/(mol K)]
static RealT constexpr R_univ_gconst = avogad * k_boltz;

static inline RealT m_dry = 0.028964;
static inline RealT grav = 9.80665;
static inline RealT cp_dry = 1004.64;

void init_constants(RealT gravity, RealT mol_weight_dry_air, RealT heat_capacity_dry_air)
{
  grav   = gravity;
  m_dry  = mol_weight_dry_air;
  cp_dry = heat_capacity_dry_air;
}

};
