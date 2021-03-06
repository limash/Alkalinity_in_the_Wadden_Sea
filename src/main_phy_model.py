import numpy as np
import xarray as xr
import brom_functions as bf

ds = xr.open_dataset('wadden_sea_out.nc')
df = ds.to_dataframe()
levelcntr = df.groupby('levelcntr').get_group(0.625)
levelface = levelcntr.groupby('levelface').get_group(0)

par = levelface['par'].values.astype(np.float64)
temperature = levelface['temperature'].values.astype(np.float64)
nh4_data = levelface['ammonium'].values.astype(np.float64)
no3_data = levelface['nitrate'].values.astype(np.float64)
po4_data = levelface['phosphate'].values.astype(np.float64)
si_data = levelface['silicate'].values.astype(np.float64)
o2_data = levelface['oxygen'].values.astype(np.float64)
chl_a_data = levelface['chl_a'].values.astype(np.float64)
om_flux = levelface['om_flux'].values[0:365].astype(np.float64)

# initial data
# some common variables
depth, k, latitude, days = (0.625, 0, 54.88, np.arange(0, 364, 1))
# nutrients
nh4 = np.zeros(365)
nh4[0] = 2
no2 = np.zeros(365)
no2[0] = 1.5
no3 = np.zeros(365)
no3[0] = no3_data[0]
si = np.zeros(365)
si[0] = si_data[0]
po4 = np.zeros(365)
po4[0] = po4_data[0]
o2 = np.zeros(365)
o2[0] = o2_data[0]
# phy
phy = np.zeros(365)
phy[0] = 90
# daily irradiance, convertion microM per second to M per day
irradiance = par*86400/1000000
# het
het = np.zeros(365)
het[0] = 90
# om
pom = np.zeros(365)
pom[0] = 100
dom = np.zeros(365)
dom[0] = 50

# initial parameters values
# horizontal advection
k_mix = 145
# phy
knh4_lim, knox_lim, ksi_lim, kpo4_lim = (0.1, 0.1, 0.1, 0.1)
pbm, alpha, kexc, kmortality = (8, 0.04, 0.015, 0.0001)
# het
k_het_phy_gro, k_het_phy_lim = (0.2, 0.4)
k_het_pom_gro = k_het_phy_gro
k_het_pom_lim = k_het_phy_lim
k_het_res, k_het_mort, uz, hz = (0.015, 0.1, 0.5, 0.5)
# nitrification
k_nfix, k_nitrif1, k_nitrif2, o2s_nf = (0.4, 0.1, 0.1, 5)
k_anammox, o2s_dn = (0.8, 10)
# om respiration
k_pom_dom, k_omox_o2, tref = (0.15, 1, 0)
k_dom_ox, k_pom_ox = (0.01, 0.002)

chl_a, phy_dgrate, rations = bf.calculate(
   depth, k, k_mix, latitude,
   days, temperature,
   nh4, no2, no3, si, po4, o2,
   nh4_data, no3_data, si_data, po4_data, o2_data,
   phy, par, irradiance,
   knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmortality,
   het, k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim,
   k_het_res, k_het_mort, uz, hz,
   k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
   pom, dom, om_flux,
   k_pom_dom, k_omox_o2, tref, k_dom_ox, k_pom_ox)
