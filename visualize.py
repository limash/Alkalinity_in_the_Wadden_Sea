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
sns.set()

#sets up pandas table display
pd.set_option('display.width', 500)
pd.set_option('display.max_columns', 100)
pd.set_option('display.notebook_repr_html', True)
pd.options.mode.chained_assignment = None
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
#rc('font', **font)



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

def returndate(datestring):
    return datetime.strptime(datestring,"%Y-%m-%d %H:%M:%S")

def plotTA(biogeodata):
    fig, ax = plt.subplots(figsize=(12,5))
    Time = biogeodata.Datetime.map(returndate).values
    TA = biogeodata.TA.values
    TAfromS = biogeodata.TAfromS.values
    size = 14
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.scatter(Time, TA, label= 'Total Alkalinity, measured', s=size)
    ax.scatter(Time, TAfromS,
               label='\nTotal Alkalinity, \n calculated from salinity', s=size)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    ax.legend(fontsize = 14,loc = 'upper left')
    plt.ylabel('Total Alkalinity, $ \mu M$')
    plt.show()
def plot_intro():
    north7 = pd.read_csv("HafniaDataNorth7Shamil.csv")
    north7 = treatbiogeodata(north7)
    fig = plotTA(north7)

def plot_alkalinity_flux_low_high():
    import xarray as xr
    base_path = 'data/low_sulfate_reduction_rate'
    ds1 = xr.open_dataset('{}/2_po75-25_di1e-9/water.nc'.format(base_path))
    ds2 = xr.open_dataset('{}/3_po75-25_di2e-9/water.nc'.format(base_path))
    ds3 = xr.open_dataset('{}/4_po75-25_di5e-9/water.nc'.format(base_path))
    ds4 = xr.open_dataset('{}/5_po75-25_di10e-9/water.nc'.format(base_path))
    ds5 = xr.open_dataset('{}/6_po75-25_di15e-9/water.nc'.format(base_path))
    ds6 = xr.open_dataset('{}/7_po75-25_di20e-9/water.nc'.format(base_path))
    ds7 = xr.open_dataset('{}/8_po75-25_di25e-9/water.nc'.format(base_path))
    ds8 = xr.open_dataset('{}/9_po75-25_di30e-9/water.nc'.format(base_path))
    ds9 = xr.open_dataset('{}/10_po75-25_di35e-9/water.nc'.format(base_path))

    ds1_high = xr.open_dataset('data/high_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2_high = xr.open_dataset('data/high_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3_high = xr.open_dataset('data/high_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4_high = xr.open_dataset('data/high_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    ds5_high = xr.open_dataset('data/high_sulfate_reduction_rate/6_po75-25_di15e-9/water.nc')
    ds6_high = xr.open_dataset('data/high_sulfate_reduction_rate/7_po75-25_di20e-9/water.nc')
    ds7_high = xr.open_dataset('data/high_sulfate_reduction_rate/8_po75-25_di25e-9/water.nc')
    ds8_high = xr.open_dataset('data/high_sulfate_reduction_rate/9_po75-25_di30e-9/water.nc')
    ds9_high = xr.open_dataset('data/high_sulfate_reduction_rate/10_po75-25_di35e-9/water.nc')

    strt ='2011-01-01'
    stp = '2011-12-31'
    alk_year = []
    alkflux_bottom_year = []

    for i,ds in enumerate((ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9),start = 0):
        alk_df = ds['B_C_Alk'].to_dataframe()
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625).reset_index('z',drop = True)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5).reset_index('z_faces',drop = True)

        alk_year.append(alk[strt:stp])
        alkflux_bottom_year.append(alkflux_bottom[strt:stp])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i]['B_C_Alk'] = alk_year[i]['B_C_Alk']-alk_year[i]['B_C_Alk'].min()

   
    alk_year_high, alkflux_bottom_year_high = [],[]
    for i,ds_high in enumerate((ds1_high, ds2_high, ds3_high, ds4_high, ds5_high, ds6_high, ds7_high, ds8_high, ds9_high),start = 0):
        alk_df = ds_high['B_C_Alk'].to_dataframe()
        alkflux_df = ds_high['B_C_Alk   _flux'].to_dataframe()
        alk_high = alk_df.groupby('z').get_group(0.625).reset_index('z',drop=True)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5).reset_index('z_faces',drop=True)
        alk_year_high.append(alk_high[strt:stp])
        alkflux_bottom_year_high.append(alkflux_bottom[strt:stp])
        alk_year_high[i] = alk_year_high[i].reset_index()
        alkflux_bottom_year_high[i] = alkflux_bottom_year_high[i].reset_index()
        alk_year_high[i]['B_C_Alk'] = alk_year_high[i]['B_C_Alk']-alk_year_high[i]['B_C_Alk'].min()


    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(2, 2, 1) # row-col-num
    ax1 = fig.add_subplot(2, 2, 2) # row-col-num
    ax_2 = fig.add_subplot(2, 2, 3) # row-col-num
    ax_3 = fig.add_subplot(2, 2, 4) # row-col-num

    lnw = 1.5
    fntsz = 16
    labels = [r'$1e-9$',r'$2e-9$',r'$5e-9$',r'$10e-9$',r'$15e-9$',r'$20e-9$',r'$25e-9$',r'$30e-9$',r'$35e-9$']

    for n in range(0,9):
        ax.plot(alkflux_bottom_year[n]['time'], alkflux_bottom_year[n]['B_C_Alk   _flux'], linewidth=lnw, label=labels[n])
        ax1.plot(alkflux_bottom_year_high[n]['time'], alkflux_bottom_year_high[n]['B_C_Alk   _flux'], linewidth=lnw, label=labels[n])

        ax_2.plot(alk_year[n]['time'], alk_year[n]['B_C_Alk'], linewidth=2, label=labels[n])
        ax_3.plot(alk_year_high[n]['time'], alk_year_high[n]['B_C_Alk'], linewidth=2, label=labels[n])

    # --- add title and axis labels
    ax.set_title('Low sulfate reduction', fontsize=fntsz)
    ax1.set_title('High sulfate reduction', fontsize=fntsz)


    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)
    ax1.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)

    ax_2.set_ylabel('Relative TA, mmol m$^{-3}$', fontsize=fntsz)
    ax_3.set_ylabel('Relative TA, mmol m$^{-3}$', fontsize=fntsz)

    ax.legend(loc='best', title='$kz_{dispersion}$, m$^2$ s$^{-1}$')

    # --- improve the layout   
    for axis in (ax,ax1,ax_2,ax_3):
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    fig.tight_layout(pad=1)
    plt.show()

plot_alkalinity_flux_low_high()

def plot_alkalinity_flux_sulfur_oxidation():
    import xarray as xr
  
    ds0 = xr.open_dataset('data/different_sulfur_oxidation/high/water.nc')
    ds1 = xr.open_dataset('data/different_sulfur_oxidation/low/water.nc')
    ds2 = xr.open_dataset('data/different_sulfur_oxidation/regular/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds0, ds1, ds2):
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
    ax = fig.add_subplot(1, 2, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'high')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'low')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'regular')
    
    # --- add title and axis labels
    #ax.set_title('Alkalinity fluxes')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    #ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='upper left', title='Sulfur compounds oxidation rates')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    ax1 = fig.add_subplot(1, 2, 2) # row-col-num
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax1.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'high')
    ax1.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'low')
    ax1.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'regular')
    # --- add title and axis labels
    #ax1.set_title('Alkalinity')
    ax1.set_ylabel('TA increment, mmol m$^{-3}$', fontsize=16)
    #ax1.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax1.legend(loc='upper left', title='Sulfur compounds oxidation rates')
    # --- add grid – not in default classic style
    ax1.grid(True)
    fig.tight_layout(pad=1)

def plot_alkalinity_flux_porosities():
    import xarray as xr
    
  
    ds0 = xr.open_dataset('data/different_porosities/4_po45-25_di10e-9/water.nc')
    ds1 = xr.open_dataset('data/different_porosities/0_po55-25_di10e-9/water.nc')
    ds2 = xr.open_dataset('data/different_porosities/1_po65-25_di10e-9/water.nc')
    ds3 = xr.open_dataset('data/different_porosities/2_po75-25_di10e-9/water.nc')
    ds4 = xr.open_dataset('data/different_porosities/3_po85-25_di10e-9/water.nc')
    
    alk_year = []
    alkflux_bottom_year = []
    
    i = 0
    for ds in (ds0, ds1, ds2, ds3, ds4):
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
    ax = fig.add_subplot(1, 2, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'45-25')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'55-25')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'65-25')
    ax.plot(alkflux_bottom_year[3]['time'], alkflux_bottom_year[3]['B_C_Alk   _flux'], linewidth=2, label=r'75-25')
    ax.plot(alkflux_bottom_year[4]['time'], alkflux_bottom_year[4]['B_C_Alk   _flux'], linewidth=2, label=r'85-25')
    
    # --- add title and axis labels
    #ax.set_title('Alkalinity fluxes')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    #ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='upper left', title='Porosities')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    ax1 = fig.add_subplot(1, 2, 2) # row-col-num
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax1.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'45-25')
    ax1.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'55-25')
    ax1.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'65-25')
    ax1.plot(alk_year[3]['time'], alk_year[3]['B_C_Alk'], linewidth=2, label=r'75-25')
    ax1.plot(alk_year[4]['time'], alk_year[4]['B_C_Alk'], linewidth=2, label=r'85-25')
    # --- add title and axis labels
    #ax1.set_title('Alkalinity')
    ax1.set_ylabel('TA increment, mmol m$^{-3}$', fontsize=16)
    #ax1.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax1.legend(loc='upper left', title='Porosities')
    # --- add grid – not in default classic style
    ax1.grid(True)
    fig.tight_layout(pad=1)

def plot_alkalinity_flux_porosities_2():
    import xarray as xr
    
  
    ds0 = xr.open_dataset('data/different_porosities_2/0_po75-05_di10e-9/water.nc')
    ds1 = xr.open_dataset('data/different_porosities_2/1_po75-15_di10e-9/water.nc')
    ds2 = xr.open_dataset('data/different_porosities_2/2_po75-25_di10e-9/water.nc')
    ds3 = xr.open_dataset('data/different_porosities_2/3_po75-35_di10e-9/water.nc')
    
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
    ax = fig.add_subplot(1, 2, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'75-05')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'75-15')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'75-25')
    ax.plot(alkflux_bottom_year[3]['time'], alkflux_bottom_year[3]['B_C_Alk   _flux'], linewidth=2, label=r'75-35')
    
    # --- add title and axis labels
    #ax.set_title('Alkalinity fluxes')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    #ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='upper left', title='Porosities')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    ax1 = fig.add_subplot(1, 2, 2) # row-col-num
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax1.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'75-05')
    ax1.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'75-15')
    ax1.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'75-25')
    ax1.plot(alk_year[3]['time'], alk_year[3]['B_C_Alk'], linewidth=2, label=r'75-35')
    # --- add title and axis labels
    #ax1.set_title('Alkalinity')
    ax1.set_ylabel('TA increment, mmol m$^{-3}$', fontsize=16)
    #ax1.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax1.legend(loc='upper left', title='Porosities')
    # --- add grid – not in default classic style
    ax1.grid(True)
    fig.tight_layout(pad=1)

def plot_alkalinity_flux_porosities_3():
    import xarray as xr
    
    ds0 = xr.open_dataset('data/different_porosities_3/0_po63-32_di10e-9/water.nc')
    ds1 = xr.open_dataset('data/different_porosities_3/1_po70-28_di10e-9/water.nc')
    ds2 = xr.open_dataset('data/different_porosities_3/2_po75-25_di10e-9/water.nc')
    ds3 = xr.open_dataset('data/different_porosities_3/3_po82-21_di10e-9/water.nc')
    
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
    ax = fig.add_subplot(1, 2, 1) # row-col-num
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.plot(alkflux_bottom_year[0]['time'], alkflux_bottom_year[0]['B_C_Alk   _flux'], linewidth=2, label=r'63-32')
    ax.plot(alkflux_bottom_year[1]['time'], alkflux_bottom_year[1]['B_C_Alk   _flux'], linewidth=2, label=r'70-28')
    ax.plot(alkflux_bottom_year[2]['time'], alkflux_bottom_year[2]['B_C_Alk   _flux'], linewidth=2, label=r'75-25')
    ax.plot(alkflux_bottom_year[3]['time'], alkflux_bottom_year[3]['B_C_Alk   _flux'], linewidth=2, label=r'82-21')
    
    # --- add title and axis labels
    #ax.set_title('Alkalinity fluxes')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    #ax.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='upper left', title='Porosities')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    ax1 = fig.add_subplot(1, 2, 2) # row-col-num
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax1.plot(alk_year[0]['time'], alk_year[0]['B_C_Alk'], linewidth=2, label=r'63-32')
    ax1.plot(alk_year[1]['time'], alk_year[1]['B_C_Alk'], linewidth=2, label=r'70-28')
    ax1.plot(alk_year[2]['time'], alk_year[2]['B_C_Alk'], linewidth=2, label=r'75-25')
    ax1.plot(alk_year[3]['time'], alk_year[3]['B_C_Alk'], linewidth=2, label=r'82-21')
    # --- add title and axis labels
    #ax1.set_title('Alkalinity')
    ax1.set_ylabel('TA increment, mmol m$^{-3}$', fontsize=16)
    #ax1.set_xlabel('Month', fontsize=16)
    # --- plot a legend in the best location
    ax1.legend(loc='upper left', title='Porosities')
    # --- add grid – not in default classic style
    ax1.grid(True)
    fig.tight_layout(pad=1)

def plot_alk_sulfur_fluxes():
    import xarray as xr

    ds1 = xr.open_dataset('data/low_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2 = xr.open_dataset('data/low_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3 = xr.open_dataset('data/low_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4 = xr.open_dataset('data/low_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')

    alkflux_bottom_july = []
    nh4flux_bottom_july = []
    no2flux_bottom_july = []
    no3flux_bottom_july = []
    po4flux_bottom_july = []
    so4flux_bottom_july = []

    alk_mean = []
    nh4_mean = []
    no2_mean = []
    no3_mean = []
    po4_mean = []
    so4_mean = []

    i = 0
    for ds in (ds1, ds2, ds3, ds4):
        alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
        nh4flux_df = ds['B_NUT_NH4 _flux'].to_dataframe()
        no2flux_df = ds['B_NUT_NO2 _flux'].to_dataframe()
        no3flux_df = ds['B_NUT_NO3 _flux'].to_dataframe()
        po4flux_df = ds['B_NUT_PO4 _flux'].to_dataframe()
        so4flux_df = ds['B_S_SO4   _flux'].to_dataframe()

        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5)
        nh4flux_bottom = nh4flux_df.groupby('z_faces').get_group(2.5)
        no2flux_bottom = no2flux_df.groupby('z_faces').get_group(2.5)
        no3flux_bottom = no3flux_df.groupby('z_faces').get_group(2.5)
        po4flux_bottom = po4flux_df.groupby('z_faces').get_group(2.5)
        so4flux_bottom = so4flux_df.groupby('z_faces').get_group(2.5)

        alkflux_bottom_july.append(alkflux_bottom.loc['2011-07-01':'2011-08-01'])
        nh4flux_bottom_july.append(nh4flux_bottom.loc['2011-07-01':'2011-08-01'])
        no2flux_bottom_july.append(no2flux_bottom.loc['2011-07-01':'2011-08-01'])
        no3flux_bottom_july.append(no3flux_bottom.loc['2011-07-01':'2011-08-01'])
        po4flux_bottom_july.append(po4flux_bottom.loc['2011-07-01':'2011-08-01'])
        so4flux_bottom_july.append(so4flux_bottom.loc['2011-07-01':'2011-08-01'])

        alkflux_bottom_july[i] = alkflux_bottom_july[i].reset_index()
        nh4flux_bottom_july[i] = nh4flux_bottom_july[i].reset_index()
        no2flux_bottom_july[i] = no2flux_bottom_july[i].reset_index()
        no3flux_bottom_july[i] = no3flux_bottom_july[i].reset_index()
        po4flux_bottom_july[i] = po4flux_bottom_july[i].reset_index()
        so4flux_bottom_july[i] = so4flux_bottom_july[i].reset_index()

        alk_mean.append(alkflux_bottom_july[i]['B_C_Alk   _flux'].mean())
        nh4_mean.append(nh4flux_bottom_july[i]['B_NUT_NH4 _flux'].mean())
        no2_mean.append(no2flux_bottom_july[i]['B_NUT_NO2 _flux'].mean())
        no3_mean.append(no3flux_bottom_july[i]['B_NUT_NO3 _flux'].mean())
        po4_mean.append(po4flux_bottom_july[i]['B_NUT_PO4 _flux'].mean())
        so4_mean.append(so4flux_bottom_july[i]['B_S_SO4   _flux'].mean())

        i += 1

    alk = np.array(alk_mean)
    nh4 = np.array(nh4_mean)
    no2 = np.array(no2_mean)
    no3 = np.array(no3_mean)
    po4 = np.array(po4_mean)
    so4 = np.array(so4_mean)

    h2s__flux_bottom_july = []
    s0___flux_bottom_july = []
    s2o3_flux_bottom_july = []

    h2s__mean = []
    s0___mean = []
    s2o3_mean = []

    i = 0
    for ds in (ds1, ds2, ds3, ds4):
        h2s__flux_df = ds['B_S_S0    _flux'].to_dataframe()
        s0___flux_df = ds['B_S_H2S   _flux'].to_dataframe()
        s2o3_flux_df = ds['B_S_S2O3  _flux'].to_dataframe()

        h2s__flux_bottom = h2s__flux_df.groupby('z_faces').get_group(2.5)
        s0___flux_bottom = s0___flux_df.groupby('z_faces').get_group(2.5)
        s2o3_flux_bottom = s2o3_flux_df.groupby('z_faces').get_group(2.5)

        h2s__flux_bottom_july.append(h2s__flux_bottom.loc['2011-07-01':'2011-08-01'])
        s0___flux_bottom_july.append(s0___flux_bottom.loc['2011-07-01':'2011-08-01'])
        s2o3_flux_bottom_july.append(s2o3_flux_bottom.loc['2011-07-01':'2011-08-01'])

        h2s__flux_bottom_july[i] = h2s__flux_bottom_july[i].reset_index()
        s0___flux_bottom_july[i] = s0___flux_bottom_july[i].reset_index()
        s2o3_flux_bottom_july[i] = s2o3_flux_bottom_july[i].reset_index()

        h2s__mean.append(h2s__flux_bottom_july[i]['B_S_S0    _flux'].mean())
        s0___mean.append(s0___flux_bottom_july[i]['B_S_H2S   _flux'].mean())
        s2o3_mean.append(s2o3_flux_bottom_july[i]['B_S_S2O3  _flux'].mean())

        i += 1

    h2s = np.array(h2s__mean)
    s0 = np.array(s0___mean)
    s2o3 = np.array(s2o3_mean)

    alk_calc = nh4-no2-no3-po4-2*so4
    s_total = h2s + s0 + s2o3
    x = np.array([1e-9, 2e-9, 5e-9, 10e-9])

    plt.style.use('classic')
    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.plot(x, alk, linewidth=2, label=r'alkalinity flux')
    ax.plot(x, s_total, linewidth=2, label=r'sulfur flux')
    # --- add title and axis labels
    #ax.set_title('Fluxes')
    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=16)
    ax.set_xlabel('$kz_{dispersion}$, m$^2$ s$^{-1}$', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='upper left', title='Fluxes')
    # --- add grid – not in default classic style
    ax.grid(True)


