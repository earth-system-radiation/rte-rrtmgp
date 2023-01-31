
#pragma once

#include "rrtmgp_const.h"

// -----------------------------------------
// Physical constants, 2018 SI defintion of metric system
//   doi:10.1088/1681-7575/aa950a (see also https://www.nist.gov/si-redefinition/meet-constants)
// Boltzmann constant [J/K] = [(kg m^2)/(K s^2)]
real constexpr k_boltz = 1.380649e-23;

//  molecular weight of water [kg/mol]
real constexpr m_h2o =  0.018016;

// Avogadro's number [molec/mol]
real constexpr avogad = 6.02214076e23;

// Universal gas constant [J/(mol K)]
real constexpr R_univ_gconst = avogad * k_boltz;


extern real m_dry;
extern real grav;
extern real cp_dry;


void init_constants(real gravity, real mol_weight_dry_air, real heat_capacity_dry_air);


