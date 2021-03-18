import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sb
import colorcet as cc

def mae(diff, col_dim):
    #
    # Mean absolute error
    #
    return(xr.ufuncs.fabs(diff).mean(dim=col_dim))
def rms(diff, col_dim):
    #
    # Root mean square error
    #
    return(xr.ufuncs.sqrt(xr.ufuncs.square(diff).mean(dim=col_dim)))

def make_comparison_plot(variants, labels, reference, vscale, col_dim="site", lay_dim="layer"):
    #
    # Make a plot comparing differences with respect to reference
    #
    if type(variants) is not list:
        make_comparison_plot([variants], [labels], reference, vscale)
    else:
        for i in np.arange(0, len(variants)):
            delta = variants[i] - reference
            plt.plot(mae(delta, col_dim),
                     vscale, '-',\
                     color=cols[i], label=labels[i])
            plt.plot(rms(delta, col_dim),
                     vscale, '--',\
                     color=cols[i])
        # Reverse vertical ordering
        plt.legend()
        # Reverse vertical ordering
        plt.ylim(vscale.max(), vscale.min())

def construct_lbl_esgf_name(var, esgf_node="dkrz"):
    #
    # For a given variable name, provide the OpenDAP URL for the LBLRM line-by-line results
    #
    prefix = "http://esgf3.dkrz.de/thredds/dodsC/cmip6/RFMIP/AER/LBLRTM-12-8/rad-irf/r1i1p1f1/Efx/"
    if(esgf_node is "llnl"): prefix = "http://esgf-data1.llnl.gov/thredds/dodsC/css03_data/CMIP6/RFMIP/AER/LBLRTM-12-8/rad-irf/r1i1p1f1/Efx/"
    return(prefix + var + "/gn/v20190514/" + var + "_Efx_LBLRTM-12-8_rad-irf_r1i1p1f1_gn.nc")

########################################################################
if __name__ == '__main__':
    gp = xr.open_dataset("test_atmospheres.nc")
    #
    # Does the flux plus the Jacobian equal a calculation with perturbed surface temperature?
    #
    gp['lw_flux_up_from_deriv'] = gp.lw_flux_up  + gp.lw_jaco_up
    gp.lw_flux_up_from_deriv.attrs = {"description":"LW flux up, surface T+1K, computed from Jacobian"}
    lbl = xr.open_mfdataset([construct_lbl_esgf_name(v) for v in ["rsd", "rsu", "rld", "rlu"]],
                            combine = "by_coords").sel(expt=0)
    ########################################################################
    #
    # The RFMIP cases are on an irregular pressure grid so we can't compute errors
    #   as a function of pressure
    # Interpolate onto a uniform pressure level - every 10 hPa
    #   plev varies from site to site - could the interpolation be made more concise?
    #
    plevs = np.arange(0, lbl.plev.max(), 1000)
    plevs[0] = lbl.plev.min().values
    plev = xr.DataArray(plevs, dims=['plev'], coords={'plev': plevs})
    lbli = xr.concat([lbl.sel(site=i).reset_coords().swap_dims({"level":"plev"}).interp(plev=plev) for i in np.arange(0, lbl.site.size)], dim = 'site')
    gpi  = xr.concat([ gp.sel(site=i).rename({"pres_level":"plev"}).\
                                      reset_coords().swap_dims({"level":"plev"}).interp(plev=plev) for i in np.arange(0,  gp.site.size)], dim = 'site')

    cols = cc.glasbey_dark
    plev.load()
    gpi.load()
    lbli.load()
    with PdfPages('validation-figures.pdf') as pdf:
        ########################################################################
        # Longwave
        ########################################################################
        #
        # Accuracy - 3-angle and single-angle
        #
        variants = [[gpi.lw_flux_dn, gpi.lw_flux_dn_optang, gpi.lw_flux_dn_3ang, gpi.lw_flux_dn_2str],
                    [gpi.lw_flux_up, gpi.lw_flux_up_optang, gpi.lw_flux_up_3ang, gpi.lw_flux_up_2str],
                    [gpi.lw_flux_net,
                    gpi.lw_flux_dn_optang - gpi.lw_flux_up_optang,
                    gpi.lw_flux_dn_3ang   - gpi.lw_flux_up_3ang,
                    gpi.lw_flux_dn_2str   - gpi.lw_flux_up_2str]]
        refs = [lbli.rld,  lbli.rlu,  lbli.rld -  lbli.rlu]
        titles = ["Accuracy wrt LBLRTM: LW down", "Accuracy wrt LBLRTM: LW up", "Accuracy: LW net"]
        for v, r, t in zip(variants, refs, titles):
            make_comparison_plot(v, \
                                 labels = ["default", "optimal-angle", "3-angle", "2-stream"], \
                                 reference = r,            \
                                 vscale = plev/100.)
            plt.ylabel("Pressure (Pa)")
            plt.xlabel("Error (W/m2), solid=mean, dash=RMS")
            plt.title(t)
            pdf.savefig()
            plt.close()

        #
        # Variants on the default - should be near but not precisely 0
        #
        # Level temperatures not provided
        # Using scattering representations of optical properties
        #   to do a no-scattering calculation
        #
        variants = [[gpi.lw_flux_dn_notlev,  gpi.lw_flux_dn_2str],
                    [gpi.lw_flux_up_notlev,  gpi.lw_flux_up_2str]]
        refs =     [gpi.lw_flux_dn,         gpi.lw_flux_up]
        titles =   ["Variants: LW down",    "Variants: LW up"]
        for v, r, t in zip(variants, refs, titles):
            make_comparison_plot(v, \
                                 labels = ["no-tlev", "2str"], \
                                 reference = r,            \
                                 vscale = plev/100.)
            plt.ylabel("Pressure (Pa)")
            plt.xlabel("Difference (W/m2), solid=mean, dash=RMS")
            plt.title(t)
            pdf.savefig()
            plt.close()
        ########################################################################
        # Shortwave
        #   These are in W/m2 so profiles with high run count for more.
        #   Would they be better scaled to incoming solar flux?
        ########################################################################
        #
        # Accuracy comparison
        #
        variants = [gpi.sw_flux_dn, gpi.sw_flux_up]
        refs =     [lbli.rsd,       lbli.rsu]
        titles = ["Accuracy: SW down", "Accuracy: SW up"]
        for v, r, t in zip(variants, refs, titles):
            make_comparison_plot(v, \
                                 labels = "default", \
                                 reference = r, \
                                 vscale = plev/100.)
            plt.ylabel("Pressure (Pa)")
            plt.xlabel("Error (W/m2), solid=mean, dash=RMS")
            plt.title(t)
            pdf.savefig()
            plt.close()
