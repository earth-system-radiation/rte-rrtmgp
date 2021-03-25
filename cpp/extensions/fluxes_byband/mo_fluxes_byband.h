#pragma once
#include "mo_fluxes.h"
#include "mo_fluxes_byband_kernels.h"
#include "mo_optical_props.h"
#include "const.h"
//#include "mo_fluxes_broadband_kernels.h"
#include <iomanip>


class FluxesByband : public FluxesBroadband {
    public:
        real3d bnd_flux_up;
        real3d bnd_flux_dn;
        real3d bnd_flux_dn_dir;
        real3d bnd_flux_net;

    void reduce(real3d const &gpt_flux_up, const real3d &gpt_flux_dn, OpticalProps const &spectral_disc,
                bool top_at_1, real3d const &gpt_flux_dn_dir=real3d()) {
        int ncol = size(gpt_flux_up,1);
        int nlev = size(gpt_flux_up,2);
        int ngpt = size(gpt_flux_up,3);
        int nbnd = size(this->bnd_flux_up,3);

        // Base clase reduce
        FluxesBroadband::reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir);
        // Reduce byband fluxes
        auto band2gpt = spectral_disc.band2gpt;
        if (allocated(this->bnd_flux_up    )) { sum_byband(ncol, nlev, ngpt, nbnd, band2gpt, gpt_flux_up,     this->bnd_flux_up    ); }
        if (allocated(this->bnd_flux_dn    )) { sum_byband(ncol, nlev, ngpt, nbnd, band2gpt, gpt_flux_dn,     this->bnd_flux_dn    ); }
        if (allocated(this->bnd_flux_dn_dir)) { sum_byband(ncol, nlev, ngpt, nbnd, band2gpt, gpt_flux_dn_dir, this->bnd_flux_dn_dir); }
        // Compute net fluxes
        if (allocated(this->bnd_flux_net)) {
            if (allocated(this->bnd_flux_dn) && allocated(this->bnd_flux_up)) {
                net_byband(ncol, nlev, nbnd, this->bnd_flux_dn, this->bnd_flux_up, this->bnd_flux_net);
            } else {
                stoprun("reduce: bnd_flux_net requested but bnd_flux_dn or bnd_flux_up not allocated.");
            }
        }
    }
    
};
