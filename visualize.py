# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import numpy as np # imports a fast numerical programming library
import matplotlib.pyplot as plt #sets up plotting under plt
from matplotlib import rc
import matplotlib.dates as mdates
import pandas as pd #lets us handle data as dataframes
from datetime import datetime
import seaborn as sns #sets up styles and gives us more plotting options

#sets up pandas table display
pd.set_option('display.width', 500)
pd.set_option('display.max_columns', 100)
pd.set_option('display.notebook_repr_html', True)
pd.options.mode.chained_assignment = None
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
rc('font', **font)
sns.set()


def addseason(datestring):
    item = datetime.strptime(datestring,"%Y-%m-%d %H:%M:%S")
    if 2 < item.month < 6:
        return 'spring' # spring
    elif 5 < item.month < 9:
        return 'summer' # summer
    elif 8 < item.month < 12:
        return 'autumn' # autumn
    else:
        return 'winter' # winter

def addmonth(datestring):
    item = datetime.strptime(datestring,"%Y-%m-%d %H:%M:%S")
    year = {1:'january', 2:'february', 3:'march', 4:'april',
            5:'may', 6:'june', 7:'july', 8:'august', 9:'september',
            10: 'october', 11: 'november', 12: 'december'}
    return year[item.month]

def addseconds(datestring):
    t = datetime.strptime(datestring,"%Y-%m-%d %H:%M:%S")
    return (t - datetime(2017,1,1)).total_seconds()

def addsecondstolvl(datestring):
    datestring = datestring[0:-4]
    t = datetime.strptime(datestring,"%d/%m/%Y %H:%M:%S")
    return (t - datetime(2017,1,1)).total_seconds()

def addband(longitude):
    if 6.95 < longitude < 7:
        return 1
    else:
        return 0

def addphase(seconds):
    period = (12 * 60 * 60) + (25.2 * 60) 
    half_period = period / 2
    startphase = (half_period / 12) * 8
    modulus = (seconds - startphase) % period
    if modulus < half_period:
        return 'low'
    else:
        return 'high'

def addphase_2(slev):
    if slev > 0:
        return 'high'
    else:
        return 'low'

def commafix(string):
    return float(string.replace(',','.'))

def calculateTA(method,t,s):
    if method == 'Bellerby':
        if s >= 34.65:
            return 66.96 * s - 36.803 # Bellerby & Canoba
        else:
            return 3887 - 46.25 * s # Borges & Frankignoulle & Canoba
    elif method == 'Millero':
        if t < 20:
            return s / 35 * (2291 - 2.69 * (t - 20) - 0.046 * np.square(t - 20)) # Millero et al, MarChem, 1998
        else:
            return 520.1 + 51.24 * s # Millero et al, MarChem, 1998

def treatlvl(sealvldata):
    sealvldata['Seconds_since_start_of_the_year'] = sealvldata.TIME.map(addsecondstolvl)
    try:
        sealvldata['SLEV'] = sealvldata.SLEV.map(commafix)
    except AttributeError:
        pass
    sealvldata = sealvldata[sealvldata.SLEV.values[:] != -999]
    sealvldata['Phase'] = sealvldata.SLEV.map(addphase_2)
    return sealvldata

def treatbiogeodata(biogeodata):
    biogeodata['Season'] = biogeodata.Datetime.map(addseason)
    biogeodata['Month'] = biogeodata.Datetime.map(addmonth)
    biogeodata['Seconds_since_start_of_the_year'] = biogeodata.Datetime.map(addseconds)
    biogeodata['TAfromS'] = [calculateTA('Millero',t, s) for t, s in zip(biogeodata.Temperature.values,
                                                                          biogeodata.Salinity.values)]
    # biogeodata['Band'] = biogeodata.Longitude.map(addband)
    return biogeodata

def addlvlphase(biogeodata, sealvldata):
    """biogeodata and sealvldata for the current month"""
    biogeodata['SLEV'] = [np.interp(x, sealvldata.Seconds_since_start_of_the_year.values, sealvldata.SLEV.values)
                          for x in biogeodata.Seconds_since_start_of_the_year.values]
    biogeodata['Phase'] = biogeodata.SLEV.map(addphase_2)
    return biogeodata

def plotTA(biogeodata):
    fig, ax = plt.subplots(figsize=(12,5))
    Time = biogeodata.Seconds_since_start_of_the_year.values / 86400
    TA = biogeodata.TA.values
    TAfromS = biogeodata.TAfromS.values
    size = 14
    ax.scatter(Time, TA, label= 'Total Alkalinity Measured',s=size)
    ax.scatter(Time, TAfromS,
               label='\nTotal Alkalinity, \nCalculated from salinity',s=size)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    ax.legend(fontsize = 14,loc = 'upper left')
    #ax.set_title('Figure 1. Surface water TA measured and calculated from the alkalinity-salinity relation (Millero 1997),\n FerryBox measurements (54.19$^\circ$ N, 6.99$^\circ$ E, 2017)', fontsize=16, fontweight='bold')
    plt.xlabel("Day in a year")
    plt.ylabel('Total Alkalinity, $ \mu M$')

def plot_intro():
    north7 = pd.read_csv("HafniaDataNorth7Shamil.csv")
    north7 = treatbiogeodata(north7)
    fig = plotTA(north7)

def plot_alkalinity_flux_low():
    import xarray as xr
    
    ds1 = xr.open_dataset('data/low_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2 = xr.open_dataset('data/low_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3 = xr.open_dataset('data/low_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4 = xr.open_dataset('data/low_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    ds5 = xr.open_dataset('data/low_sulfate_reduction_rate/6_po75-25_di15e-9/water.nc')
    ds6 = xr.open_dataset('data/low_sulfate_reduction_rate/7_po75-25_di20e-9/water.nc')
    ds7 = xr.open_dataset('data/low_sulfate_reduction_rate/8_po75-25_di25e-9/water.nc')
    ds8 = xr.open_dataset('data/low_sulfate_reduction_rate/9_po75-25_di30e-9/water.nc')
    ds9 = xr.open_dataset('data/low_sulfate_reduction_rate/10_po75-25_di35e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9):#, ds10, ds11, ds12, ds13):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        alk_year.append(alk.loc['2011-01-01':'2011-12-31'])
        alkflux_bottom_year.append(alkflux_bottom.loc['2011-01-01':'2011-12-31'])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()
        i += 1
    plt.style.use('classic')
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'$1e-9$')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'$2e-9$')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'$5e-9$')
    ax.plot(alkflux_bottom_year[3]['time'], alkflux_bottom_year[3]['B_C_Alk   _flux'], linewidth=2, label=r'$10e-9$')
    ax.plot(alkflux_bottom_year[4]['time'], alkflux_bottom_year[4]['B_C_Alk   _flux'], linewidth=2, label=r'$15e-9$')
    ax.plot(alkflux_bottom_year[5]['time'], alkflux_bottom_year[5]['B_C_Alk   _flux'], linewidth=2, label=r'$20e-9$')
    ax.plot(alkflux_bottom_year[6]['time'], alkflux_bottom_year[6]['B_C_Alk   _flux'], linewidth=2, label=r'$25e-9$')
    ax.plot(alkflux_bottom_year[7]['time'], alkflux_bottom_year[7]['B_C_Alk   _flux'], linewidth=2, label=r'$30e-9$')
    ax.plot(alkflux_bottom_year[8]['time'], alkflux_bottom_year[8]['B_C_Alk   _flux'], linewidth=2, label=r'$35e-9$')
    
    # --- add title and axis labels
    ax.set_title('The lower bound of sulfate reduction rates')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='lower left', title='$kz_{dispersion}$')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

def plot_alkalinity_low():
    import xarray as xr
    
    ds1 = xr.open_dataset('data/low_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2 = xr.open_dataset('data/low_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3 = xr.open_dataset('data/low_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4 = xr.open_dataset('data/low_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    ds5 = xr.open_dataset('data/low_sulfate_reduction_rate/6_po75-25_di15e-9/water.nc')
    ds6 = xr.open_dataset('data/low_sulfate_reduction_rate/7_po75-25_di20e-9/water.nc')
    ds7 = xr.open_dataset('data/low_sulfate_reduction_rate/8_po75-25_di25e-9/water.nc')
    ds8 = xr.open_dataset('data/low_sulfate_reduction_rate/9_po75-25_di30e-9/water.nc')
    ds9 = xr.open_dataset('data/low_sulfate_reduction_rate/10_po75-25_di35e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9):#, ds10, ds11, ds12, ds13):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        alk_year.append(alk.loc['2011-01-01':'2011-12-31'])
        alkflux_bottom_year.append(alkflux_bottom.loc['2011-01-01':'2011-12-31'])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()
        i += 1
    plt.style.use('classic')
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'$1e-9$')
    ax.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'$2e-9$')
    ax.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'$5e-9$')
    ax.plot(alk_year[3]['time'], alk_year[3]['B_C_Alk'], linewidth=2, label=r'$10e-9$')
    ax.plot(alk_year[4]['time'], alk_year[4]['B_C_Alk'], linewidth=2, label=r'$15e-9$')
    ax.plot(alk_year[5]['time'], alk_year[5]['B_C_Alk'], linewidth=2, label=r'$20e-9$')
    ax.plot(alk_year[6]['time'], alk_year[6]['B_C_Alk'], linewidth=2, label=r'$25e-9$')
    ax.plot(alk_year[7]['time'], alk_year[7]['B_C_Alk'], linewidth=2, label=r'$30e-9$')
    ax.plot(alk_year[8]['time'], alk_year[8]['B_C_Alk'], linewidth=2, label=r'$35e-9$')
    #ax.plot(alk_year[9]['time'], alk_year[9]['B_C_Alk'], linewidth=2, label=r'$20e-9$')
    #ax.plot(alk_year[10]['time'], alk_year[10]['B_C_Alk'], linewidth=2, label=r'$24e-9$')
    #ax.plot(alk_year[11]['time'], alk_year[11]['B_C_Alk'], linewidth=2, label=r'$28e-9$')
    #ax.plot(alk_year[12]['time'], alk_year[12]['B_C_Alk'], linewidth=2, label=r'$35e-9$')
    # --- add title and axis labels
    ax.set_title('The lower bound of sulfate reduction rates')
    ax.set_ylabel('TA increment, mmol m$^{-3}$', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='lower left', title='$kz_{dispersion}$')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

def plot_alkalinity_flux_high():
    import xarray as xr
    
    ds1 = xr.open_dataset('data/high_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2 = xr.open_dataset('data/high_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3 = xr.open_dataset('data/high_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4 = xr.open_dataset('data/high_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    ds5 = xr.open_dataset('data/high_sulfate_reduction_rate/6_po75-25_di15e-9/water.nc')
    ds6 = xr.open_dataset('data/high_sulfate_reduction_rate/7_po75-25_di20e-9/water.nc')
    ds7 = xr.open_dataset('data/high_sulfate_reduction_rate/8_po75-25_di25e-9/water.nc')
    ds8 = xr.open_dataset('data/high_sulfate_reduction_rate/9_po75-25_di30e-9/water.nc')
    ds9 = xr.open_dataset('data/high_sulfate_reduction_rate/10_po75-25_di35e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9):#, ds10, ds11, ds12, ds13):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        alk_year.append(alk.loc['2011-01-01':'2011-12-31'])
        alkflux_bottom_year.append(alkflux_bottom.loc['2011-01-01':'2011-12-31'])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()
        i += 1
    plt.style.use('classic')
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'$1e-9$')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'$2e-9$')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'$5e-9$')
    ax.plot(alkflux_bottom_year[3]['time'], alkflux_bottom_year[3]['B_C_Alk   _flux'], linewidth=2, label=r'$10e-9$')
    ax.plot(alkflux_bottom_year[4]['time'], alkflux_bottom_year[4]['B_C_Alk   _flux'], linewidth=2, label=r'$15e-9$')
    ax.plot(alkflux_bottom_year[5]['time'], alkflux_bottom_year[5]['B_C_Alk   _flux'], linewidth=2, label=r'$20e-9$')
    ax.plot(alkflux_bottom_year[6]['time'], alkflux_bottom_year[6]['B_C_Alk   _flux'], linewidth=2, label=r'$25e-9$')
    ax.plot(alkflux_bottom_year[7]['time'], alkflux_bottom_year[7]['B_C_Alk   _flux'], linewidth=2, label=r'$30e-9$')
    ax.plot(alkflux_bottom_year[8]['time'], alkflux_bottom_year[8]['B_C_Alk   _flux'], linewidth=2, label=r'$35e-9$')
    # --- add title and axis labels
    ax.set_title('The upper bound of sulfate reduction rates')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='lower left', title='$kz_{dispersion}$')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

def plot_alkalinity_high():
    import xarray as xr
    
    ds1 = xr.open_dataset('data/high_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2 = xr.open_dataset('data/high_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3 = xr.open_dataset('data/high_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4 = xr.open_dataset('data/high_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    ds5 = xr.open_dataset('data/high_sulfate_reduction_rate/6_po75-25_di15e-9/water.nc')
    ds6 = xr.open_dataset('data/high_sulfate_reduction_rate/7_po75-25_di20e-9/water.nc')
    ds7 = xr.open_dataset('data/high_sulfate_reduction_rate/8_po75-25_di25e-9/water.nc')
    ds8 = xr.open_dataset('data/high_sulfate_reduction_rate/9_po75-25_di30e-9/water.nc')
    ds9 = xr.open_dataset('data/high_sulfate_reduction_rate/10_po75-25_di35e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9):#, ds10, ds11, ds12, ds13):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        alk_year.append(alk.loc['2011-01-01':'2011-12-31'])
        alkflux_bottom_year.append(alkflux_bottom.loc['2011-01-01':'2011-12-31'])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()
        i += 1
    plt.style.use('classic')
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'$1e-9$')
    ax.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'$2e-9$')
    ax.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'$5e-9$')
    ax.plot(alk_year[3]['time'], alk_year[3]['B_C_Alk'], linewidth=2, label=r'$10e-9$')
    ax.plot(alk_year[4]['time'], alk_year[4]['B_C_Alk'], linewidth=2, label=r'$15e-9$')
    ax.plot(alk_year[5]['time'], alk_year[5]['B_C_Alk'], linewidth=2, label=r'$20e-9$')
    ax.plot(alk_year[6]['time'], alk_year[6]['B_C_Alk'], linewidth=2, label=r'$25e-9$')
    ax.plot(alk_year[7]['time'], alk_year[7]['B_C_Alk'], linewidth=2, label=r'$30e-9$')
    ax.plot(alk_year[8]['time'], alk_year[8]['B_C_Alk'], linewidth=2, label=r'$35e-9$')
    #ax.plot(alk_year[9]['time'], alk_year[9]['B_C_Alk'], linewidth=2, label=r'$20e-9$')
    #ax.plot(alk_year[10]['time'], alk_year[10]['B_C_Alk'], linewidth=2, label=r'$24e-9$')
    #ax.plot(alk_year[11]['time'], alk_year[11]['B_C_Alk'], linewidth=2, label=r'$28e-9$')
    #ax.plot(alk_year[12]['time'], alk_year[12]['B_C_Alk'], linewidth=2, label=r'$35e-9$')
    # --- add title and axis labels
    ax.set_title('The upper bound of sulfate reduction rates')
    ax.set_ylabel('TA increment, mmol m$^{-3}$', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='lower left', title='$kz_{dispersion}$')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

def plot_alkalinity_porosities():
    import xarray as xr
  
    ds0 = xr.open_dataset('data/different_porosities/0_po55-25_di10e-9/water.nc')
    ds1 = xr.open_dataset('data/different_porosities/1_po65-25_di10e-9/water.nc')
    ds2 = xr.open_dataset('data/different_porosities/2_po75-25_di10e-9/water.nc')
    ds3 = xr.open_dataset('data/different_porosities/3_po85-25_di10e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds0, ds1, ds2, ds3):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        alk_year.append(alk.loc['2011-01-01':'2011-12-31'])
        alkflux_bottom_year.append(alkflux_bottom.loc['2011-01-01':'2011-12-31'])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()
        i += 1
    plt.style.use('classic')
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'55-25')
    ax.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'65-25')
    ax.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'75-25')
    ax.plot(alk_year[3]['time'], alk_year[3]['B_C_Alk'], linewidth=2, label=r'85-25')
    # --- add title and axis labels
    ax.set_title('Alkalinity')
    ax.set_ylabel('Delta', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='best')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

def plot_alkalinity_flux_porosities():
    import xarray as xr
    
  
    ds0 = xr.open_dataset('data/different_porosities/0_po55-25_di10e-9/water.nc')
    ds1 = xr.open_dataset('data/different_porosities/1_po65-25_di10e-9/water.nc')
    ds2 = xr.open_dataset('data/different_porosities/2_po75-25_di10e-9/water.nc')
    ds3 = xr.open_dataset('data/different_porosities/3_po85-25_di10e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds0, ds1, ds2, ds3):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        alk_year.append(alk.loc['2011-01-01':'2011-12-31'])
        alkflux_bottom_year.append(alkflux_bottom.loc['2011-01-01':'2011-12-31'])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()
        i += 1
    plt.style.use('classic')
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'55-25')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'65-25')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'75-25')
    ax.plot(alkflux_bottom_year[3]['time'], alkflux_bottom_year[3]['B_C_Alk   _flux'], linewidth=2, label=r'85-25')
    
    # --- add title and axis labels
    ax.set_title('Alkalinity_fluxes')
    ax.set_ylabel('Flux', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='best')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)


