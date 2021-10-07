
#pragma once

#include "rrtmgp_const.h"
#include "mo_rrtmgp_util_reorder_kernels.h"


// (x,y,z) -> (z,x,y)
void reorder123x312(int d1, int d2, int d3, real3d const &array, real3d &array_out);


// (x,y,z) -> (z,y,x)
void reorder123x321(int d1, int d2, int d3, real3d const &array, real3d &array_out);

