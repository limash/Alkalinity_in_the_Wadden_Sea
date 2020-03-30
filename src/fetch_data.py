import numpy as np
import xarray as xr


def get_data():
    ds = xr.open_dataset('data/wadden_sea_out.nc')
    df = ds.to_dataframe()

    levelcntr = df.groupby('levelcntr').get_group(0.625)
    levelface = levelcntr.groupby('levelface').get_group(0)
    levelface.describe()

    par = levelface['par'].values[0:365].astype(np.float64)
    temperature = levelface['temperature'].values[0:365].astype(np.float64)
    no3 = levelface['nitrate'].values[0:365].astype(np.float64)
    ammonium = levelface['ammonium'].values[0:365].astype(np.float64)
    po4 = levelface['phosphate'].values[0:365].astype(np.float64)
    si = levelface['silicate'].values[0:365].astype(np.float64)
    irradiance = par*86400/1000000  # convertion microM per second to M per day

    return par, temperature, no3, ammonium, po4, si, irradiance
