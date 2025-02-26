---
layout: page
title: Migrating from RRTMG
---

- **Band ordering in RRTMGP_SW** - The band ordering in the SW has changed, so users should be careful when revising RRTMG inputs like surface albedo for use in RRTMGP. The band ordering in RRTMG_SW began where the bands in RRTMG_LW ended, with the near-infrared bands first and the ultraviolet bands last. The one remaining band, which covered the infrared region (820-2600 cm-1), was placed after the UV bands. This idiosyncratic ordering has not been used in RRTMGP_SW. Instead, the bands are ordered in wavenumber, with the infrared band occurring first. Please consider this reordering when providing external input (surface albedo, aerosol and cloud optical properties) to the SW code.

- **RRTMGP's band boundaries** - A few of the band boundaries have shifted from RRTMG's boundaries, which should be kept in mind when interpreting the spectral fluxes and heating rates. The first two LW bands in RRTMGP are 10-250 cm-1 and 250-500 cm-1 (in RRTMG, the middle boundary was 350 cm-1). Also, the band boundaries at 2380 and 2600 cm-1 in RRTMG (relevant bands - 2250-2380 cm-1, 2380-2600 cm-1, 2600-3250 cm-1) have been shifted in RRTMGP to 2390 and 2680 cm-1, respectively (relevant bands - 2250-2390 cm-1, 2390-2680 cm-1, 2680-3250 cm-1).

- **Extraterrestial solar irradiance** - The solar constant in RRTMGP is 1360.9 W/m2 and is based on the NRLSSI2 model. (The default solar constant in RRTMG was ~1368 W/m2.) Both the total solar irradiance and its spectral distribution can be adjusted via calls to ty_gas_optics_rrtmgp%set_solar_variablity()
