
#include "mo_rrtmgp_constants.h"


real m_dry = 0.028964_wp;
real grav = 9.80665_wp;
real cp_dry = 1004.64_wp;


void init_constants(real gravity, real mol_weight_dry_air, real heat_capacity_dry_air) {
  grav   = gravity;
  m_dry  = mol_weight_dry_air;
  cp_dry = heat_capacity_dry_air;
}



