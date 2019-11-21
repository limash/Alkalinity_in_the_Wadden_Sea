import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.dates as mdates
import pandas as pd
from datetime import datetime
import matplotlib as mpl
import seaborn as sns
sns.set()
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

#sets up pandas table display
pd.set_option('display.width', 500)
pd.set_option('display.max_columns', 100)
pd.set_option('display.notebook_repr_html', True)
pd.options.mode.chained_assignment = None
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
#rc('font', **font)

lnw = 1.5
fntsz = 10
lgndsz = 8
mpl.rcParams['xtick.labelsize'] = lgndsz
mpl.rcParams['ytick.labelsize'] = lgndsz
strt ='2011-01-01'
stp = '2011-12-31'
x_text = 0.97
y_text = 0.98

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

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def plotTA(biogeodata):
    fig, ax = plt.subplots(figsize=cm2inch(18, 5), constrained_layout=True)
    Time = biogeodata.Datetime.map(returndate).values
    TA = biogeodata.TA.values
    TAfromS = biogeodata.TAfromS.values
    size = fntsz
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.scatter(Time, TA, label= 'Total Alkalinity, measured', s=size)
    ax.scatter(Time, TAfromS, label='Total Alkalinity, calculated from salinity', s=size)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fntsz)
    ax.legend(loc = 'upper left',fontsize = lgndsz)
    plt.ylabel('Total Alkalinity, $\mu M$')
    plt.savefig('Figure1.png', dpi = 300)
    #plt.show()

def plot_intro():
    north7 = pd.read_csv("HafniaDataNorth7Shamil.csv")
    north7 = treatbiogeodata(north7)
    plotTA(north7)

alk_var = 'B_C_Alk'
alkflux_var = 'B_C_Alk   _flux'
def get_data_time(dtsts):

    alk_year,alkflux_bottom_year  = [], []

    for i,ds in enumerate((dtsts),start = 0):
        alk_df = ds[alk_var].to_dataframe()
        alkflux_df = ds[alkflux_var].to_dataframe()
        alk = alk_df.groupby('z').get_group(0.625).reset_index('z',drop = True)
        alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5).reset_index('z_faces',drop = True)

        alk_year.append(alk[strt:stp])
        alkflux_bottom_year.append(alkflux_bottom[strt:stp])
        alk_year[i] = alk_year[i].reset_index()
        alkflux_bottom_year[i] = alkflux_bottom_year[i].reset_index()
        alk_year[i][alk_var] = alk_year[i][alk_var]-alk_year[i][alk_var].min()

    return alk_year,alkflux_bottom_year

def plot_alkalinity_flux_low_high():
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

    #ds1_high = xr.open_dataset('data/high_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    #ds2_high = xr.open_dataset('data/high_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    #ds3_high = xr.open_dataset('data/high_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    #ds4_high = xr.open_dataset('data/high_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    #ds5_high = xr.open_dataset('data/high_sulfate_reduction_rate/6_po75-25_di15e-9/water.nc')
    #ds6_high = xr.open_dataset('data/high_sulfate_reduction_rate/7_po75-25_di20e-9/water.nc')
    #ds7_high = xr.open_dataset('data/high_sulfate_reduction_rate/8_po75-25_di25e-9/water.nc')
    #ds8_high = xr.open_dataset('data/high_sulfate_reduction_rate/9_po75-25_di30e-9/water.nc')
    #ds9_high = xr.open_dataset('data/high_sulfate_reduction_rate/10_po75-25_di35e-9/water.nc')

    alk_year,alkflux_bottom_year = get_data_time([ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9])
    #alk_year_high, alkflux_bottom_year_high = get_data_time([ds1_high, ds2_high, ds3_high, ds4_high,
    #                                        ds5_high, ds6_high, ds7_high, ds8_high, ds9_high])

    fig = plt.figure(figsize=cm2inch(18, 7))
    ax = fig.add_subplot(1, 2, 1) # row-col-num
    ax1 = fig.add_subplot(1, 2, 2) # row-col-num
    #ax_2 = fig.add_subplot(2, 2, 3) # row-col-num
    #ax_3 = fig.add_subplot(2, 2, 4) # row-col-num

    labels = [r'$1e-9$',r'$2e-9$',r'$5e-9$',r'$10e-9$',r'$15e-9$',r'$20e-9$',r'$25e-9$',r'$30e-9$',r'$35e-9$']

    for n in range(0,9):
        ax.plot(alkflux_bottom_year[n]['time'], alkflux_bottom_year[n][alkflux_var], linewidth=lnw, label=labels[n])
        #ax1.plot(alkflux_bottom_year_high[n]['time'], alkflux_bottom_year_high[n][alkflux_var], linewidth=lnw, label=labels[n])

        ax1.plot(alk_year[n]['time'], alk_year[n][alk_var], linewidth=lnw, label=labels[n])
        #ax_3.plot(alk_year_high[n]['time'], alk_year_high[n][alk_var], linewidth=lnw, label=labels[n])

    # --- add title and axis labels
    #ax.set_title('Low sulfate reduction', fontsize=fntsz)
    #ax1.set_title('High sulfate reduction', fontsize=fntsz)

    ax.set_ylabel('TA fluxes, mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)
    #ax1.set_ylabel('TA fluxes, mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)

    ax1.set_ylabel('Relative TA, mmol m$^{-3}$', fontsize=fntsz)
    #ax_3.set_ylabel('Relative TA, mmol m$^{-3}$', fontsize=fntsz)

    ax.legend(loc='best', title='$kz_{dispersion}$, m$^2$ s$^{-1}$',fontsize = lgndsz, title_fontsize = lgndsz)

    labels = ('(A) ','(B)')#,'(C) ','(D)')
    # --- improve the layout
    for i,axis in enumerate((ax,ax1)):#,ax_2,ax_3)):
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=fntsz, fontweight='bold', va='top', ha='right')

    fig.tight_layout(pad=0.5)
    plt.savefig('Figure3-new.png', dpi = 300)
    #plt.show()

def plot_alkalinity_flux_sulfur_oxidation():

    ds0 = xr.open_dataset('data/different_sulfur_oxidation/high/water.nc')
    ds1 = xr.open_dataset('data/different_sulfur_oxidation/low/water.nc')
    ds2 = xr.open_dataset('data/different_sulfur_oxidation/regular/water.nc')

    alk_year,alkflux_bottom_year = get_data_time([ds0, ds1, ds2])

    fig = plt.figure(figsize=cm2inch(18, 7))
    x_text = 0.97
    y_text = 0.98
    ax = fig.add_subplot(1, 2, 1) # row-col-num
    ax1 = fig.add_subplot(1, 2, 2) # row-col-num

    labels = ['high','low','base']
    for n in range(0,3):
        ax.plot(alkflux_bottom_year[n]['time'], alkflux_bottom_year[n][alkflux_var], linewidth=lnw, label=labels[n])
        ax1.plot(alk_year[n]['time'], alk_year[n][alk_var], linewidth=lnw,label=labels[n])

    ax.set_ylabel('mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)
    ax.set_title('TA fluxes',fontsize = fntsz)
    ax1.set_ylabel('mmol m$^{-3}$', fontsize=fntsz)
    ax1.set_title('Relative Total Alkalinity', fontsize=fntsz)

    labels = ('(A) ','(B)')
    for i,axis in enumerate((ax,ax1)):
        axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=fntsz, fontweight='bold', va='top', ha='right')
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    ax.legend(loc='upper left', title='Sulfur compounds \noxidation rates',fontsize = lgndsz, title_fontsize = lgndsz)
    fig.tight_layout(pad=0.5)
    plt.savefig('Figure4.png', dpi = 300)
    #plt.show()

def plot_alkalinity_flux_porosities1_2_3():

    ds0 = xr.open_dataset('data/different_porosities/4_po45-25_di10e-9/water.nc')
    ds1 = xr.open_dataset('data/different_porosities/0_po55-25_di10e-9/water.nc')
    ds2 = xr.open_dataset('data/different_porosities/1_po65-25_di10e-9/water.nc')
    ds3 = xr.open_dataset('data/different_porosities/2_po75-25_di10e-9/water.nc')
    ds4 = xr.open_dataset('data/different_porosities/3_po85-25_di10e-9/water.nc')

    ds0_2 = xr.open_dataset('data/different_porosities_2/0_po75-05_di10e-9/water.nc')
    ds1_2 = xr.open_dataset('data/different_porosities_2/1_po75-15_di10e-9/water.nc')
    ds2_2 = xr.open_dataset('data/different_porosities_2/2_po75-25_di10e-9/water.nc')
    ds3_2 = xr.open_dataset('data/different_porosities_2/3_po75-35_di10e-9/water.nc')

    ds0_3 = xr.open_dataset('data/different_porosities_3/0_po63-32_di10e-9/water.nc')
    ds1_3 = xr.open_dataset('data/different_porosities_3/1_po70-28_di10e-9/water.nc')
    ds2_3 = xr.open_dataset('data/different_porosities_3/2_po75-25_di10e-9/water.nc')
    ds3_3 = xr.open_dataset('data/different_porosities_3/3_po82-21_di10e-9/water.nc')

    alk_year, alkflux_bottom_year = get_data_time([ds0, ds1, ds2, ds3, ds4])
    alk_year_2,alkflux_bottom_year_2 = get_data_time([ds0_2, ds1_2, ds2_2, ds3_2])
    alk_year_3,alkflux_bottom_year_3 = get_data_time([ds0_3, ds1_3, ds2_3, ds3_3])

    fig = plt.figure(figsize=cm2inch(18, 18))
    ax = fig.add_subplot(3, 2, 1) # row-col-num
    ax1 = fig.add_subplot(3, 2, 2) # row-col-num
    ax_2 = fig.add_subplot(3, 2, 3) # row-col-num
    ax1_2 = fig.add_subplot(3, 2, 4) # row-col-num
    ax1_3 = fig.add_subplot(3, 2, 6) # row-col-num
    ax_3 = fig.add_subplot(3, 2, 5) # row-col-num

    labels = ('0.45-0.25','0.55-0.25','0.65-0.25','0.75-0.25','0.85-0.25')
    labels_2 =('0.75-0.05','0.75-0.15','0.75-0.25','0.75-0.35')
    labels_3 = ('0.63-0.32','0.70-0.28','0.75-0.25','0.82-0.21')

    for n in range(0,5):
        ax.plot(alkflux_bottom_year[n]['time'], alkflux_bottom_year[n]['B_C_Alk   _flux'], linewidth=lnw,alpha = 1, label=labels[n])
        ax1.plot(alk_year[n]['time'], alk_year[n]['B_C_Alk'], linewidth=lnw, label=labels[n])

    for n in range(0,4):
        ax_2.plot(alkflux_bottom_year_2[n]['time'], alkflux_bottom_year_2[n]['B_C_Alk   _flux'], linewidth=lnw, label=labels_2[n])
        ax1_2.plot(alk_year_2[n]['time'], alk_year_2[n]['B_C_Alk'], linewidth=lnw, label= labels_2[n])

        ax_3.plot(alkflux_bottom_year_3[n]['time'], alkflux_bottom_year_3[n]['B_C_Alk   _flux'], linewidth=2, label=labels_3[n])
        ax1_3.plot(alk_year_3[n]['time'], alk_year_3[n]['B_C_Alk'], linewidth=lnw, label=labels_3[n])

    for axis in [ax,ax1,ax_2,ax1_2,ax_3,ax1_3]:
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    for axis in [ax,ax_2,ax_3]:
        axis.set_ylabel('mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)
        axis.legend(loc='upper left', title='Porosities', fontsize=lgndsz, title_fontsize = lgndsz)
        axis.set_ylim(0,27)

    for axis in [ax1,ax1_2,ax1_3]:
        axis.set_ylabel('mmol m$^{-3}$', fontsize=fntsz)
        axis.set_ylim(0,180)
    ax.set_title('TA fluxes', fontsize=fntsz)
    ax1.set_title('Relative Total Alkalinity', fontsize=fntsz)

    labels = ('(A) ','(B)','(C) ','(D)','(E)','(F)')
    for i,axis in enumerate((ax,ax1,ax_2,ax1_2,ax_3,ax1_3)):
        axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=fntsz, fontweight='bold', va='top', ha='right')

    fig.tight_layout(pad=0.5)
    plt.savefig('Figure5-new.png', dpi = 300)
    #plt.show()

def plot_alk_sulfur_fluxes():

    ds1 = xr.open_dataset('data/low_sulfate_reduction_rate/2_po75-25_di1e-9/water.nc')
    ds2 = xr.open_dataset('data/low_sulfate_reduction_rate/3_po75-25_di2e-9/water.nc')
    ds3 = xr.open_dataset('data/low_sulfate_reduction_rate/4_po75-25_di5e-9/water.nc')
    ds4 = xr.open_dataset('data/low_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')

    def get_var_data_time(dtsts,varname):
        varflux_bottom_july,var_mean = [],[]
        for i,ds in enumerate(dtsts,start = 0):
            varflux_df = ds[varname].to_dataframe()
            varflux_bottom = varflux_df.groupby('z_faces').get_group(2.5).reset_index('z_faces',drop = True)
            varflux_bottom_july.append(varflux_bottom['2011-07-01':'2011-08-01'])
            varflux_bottom_july[i] = varflux_bottom_july[i].reset_index()
            var_mean.append(varflux_bottom_july[i][varname].mean())
        return np.array(var_mean),varflux_bottom_july

    dtsts = [ds1, ds2, ds3, ds4]
    alk, alkflux_bottom_july = get_var_data_time(dtsts,'B_C_Alk   _flux')
    nh4, nh4flux_bottom_july = get_var_data_time(dtsts,'B_NUT_NH4 _flux')
    no2, no2flux_bottom_july = get_var_data_time(dtsts,'B_NUT_NO2 _flux')
    no3, no3flux_bottom_july = get_var_data_time(dtsts,'B_NUT_NO3 _flux')
    po4, po4flux_bottom_july = get_var_data_time(dtsts,'B_NUT_PO4 _flux')
    so4, so4flux_bottom_july = get_var_data_time(dtsts,'B_S_SO4   _flux')
    h2s, h2sflux_bottom_july = get_var_data_time(dtsts,'B_S_H2S   _flux')
    s0, s0flux_bottom_july = get_var_data_time(dtsts,'B_S_S0    _flux')
    s2o3, s2o3flux_july = get_var_data_time(dtsts,'B_S_S2O3  _flux')
    alk_calc = nh4-no2-no3-po4-2*so4
    s_total = h2s + s0 + 2*s2o3
    x = np.array([1e-9, 2e-9, 5e-9, 10e-9])

    fig = plt.figure(figsize=cm2inch(8.5, 6), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    ax.plot(x, alk, linewidth=lnw, label=r'alkalinity flux')
    ax.plot(x, s_total, linewidth=lnw, label=r'sulfur flux')

    ax.set_ylabel('Flux, mmol m$^{-2}$ d$^{-1}$', fontsize=fntsz)
    ax.set_xlabel('$kz_{dispersion}$, m$^2$ s$^{-1}$', fontsize=fntsz)

    ax.legend(loc='upper left', title='Fluxes',fontsize = lgndsz, title_fontsize = lgndsz)
    #fig.tight_layout(pad=0.5)
    plt.savefig('Figure6-new.png', dpi = 300)
    #plt.show()

def plot_caco3():
    import matplotlib.gridspec as gs

    ds = xr.open_dataset('data/low_sulfate_reduction_rate/5_po75-25_di10e-9/water.nc')
    #ds_ad = xr.open_dataset('data/advective_ta/1_10e-9/water.nc')

    alkflux_df = ds['B_C_Alk   _flux'].to_dataframe()
    biogrow_df = ds['B_BIO_GrowthPhy'].to_dataframe()
    omresp_df  = ds['B_BIO_DcPOM_O2'].to_dataframe()
    alk_df     = ds['B_C_Alk'].to_dataframe()

    #minus values mean alkalinity goes to the North Sea
    #nh4ta_df = ds_ad['TA_due_to_NH4'].to_dataframe()
    #no3ta_df = ds_ad['TA_due_to_NO3'].to_dataframe()
    #po4ta_df = ds_ad['TA_due_to_PO4'].to_dataframe()

    #nh4ta_df.loc[:'2010-04-01'] = 0; nh4ta_df.loc['2012-05-01':] = 0
    #nh4ta_df.loc['2010-05-01':'2011-04-01'] = 0; nh4ta_df.loc['2011-05-01':'2012-04-01'] = 0
    #no3ta_df.loc[:'2010-04-01'] = 0; no3ta_df.loc['2012-05-01':] = 0
    #no3ta_df.loc['2010-05-01':'2011-04-01'] = 0; no3ta_df.loc['2011-05-01':'2012-04-01'] = 0
    #po4ta_df.loc[:'2010-04-01'] = 0; po4ta_df.loc['2012-05-01':] = 0
    #po4ta_df.loc['2010-05-01':'2011-04-01'] = 0; po4ta_df.loc['2011-05-01':'2012-04-01'] = 0

    #nh4ta_rolled = nh4ta_df.rolling(window=30).sum()
    #no3ta_rolled = no3ta_df.rolling(window=30).sum()
    #po4ta_rolled = po4ta_df.rolling(window=30).sum()

    #nh4ta_year = nh4ta_rolled.loc['2011-01-01':'2011-12-31']
    #no3ta_year = no3ta_rolled.loc['2011-01-01':'2011-12-31']
    #po4ta_year = po4ta_rolled.loc['2011-01-01':'2011-12-31']

    #nh4ta = np.array(nh4ta_year.TA_due_to_NH4.values)
    #no3ta = np.array(no3ta_year.TA_due_to_NO3.values)
    #po4ta = np.array(po4ta_year.TA_due_to_PO4.values)
    #total = nh4ta+no3ta+po4ta
    #total = total * (-1/2.5) # depth is 2.5, minus since we return nutrients influence back

    alkflux_bottom = alkflux_df.groupby('z_faces').get_group(2.5).reset_index('z_faces',drop = True)
    omresp_bottom  = omresp_df.groupby('z').get_group(2.4749999046325684).reset_index('z',drop = True)
    biogrow_surfac = biogrow_df.groupby('z').get_group(0.625).reset_index('z',drop = True)
    alk_surface    = alk_df.groupby('z').get_group(0.625).reset_index('z',drop = True)
    alk_surface_year = alk_surface[strt:stp].reset_index()

    year = (('2011-01-01','2011-01-31'), ('2011-02-01','2011-02-28'), ('2011-03-01','2011-03-31'), ('2011-04-01','2011-04-30'),
        ('2011-05-01','2011-05-31'), ('2011-06-01','2011-06-30'), ('2011-07-01','2011-07-31'), ('2011-08-01','2011-08-31'),
        ('2011-09-01','2011-09-30'), ('2011-10-01','2011-10-31'), ('2011-11-01','2011-11-30'), ('2011-12-01','2011-12-31'))

    year_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    #year_acc_days = np.cumsum(year_days)

    alk_year_delta = []
    alk_year = []
    bio_year = []
    res_year = []
    for month in year:
        alk_delta_month = alk_surface[month[0]:month[1]]
        alk_month = alkflux_bottom[month[0]:month[1]]
        bio_month = biogrow_surfac[month[0]:month[1]]
        res_month =  omresp_bottom[month[0]:month[1]]
        alk_year_delta.append(alk_delta_month['B_C_Alk'][0])
        alk_year.append(alk_month['B_C_Alk   _flux'].mean())
        bio_year.append(bio_month['B_BIO_GrowthPhy'].mean())
        res_year.append(res_month['B_BIO_DcPOM_O2'].mean())

    bio_year_quotas = np.array(bio_year)/sum(bio_year)
    alk_year_quotas = np.array(alk_year)/sum(alk_year)
    res_year_quotas = np.array(res_year)/sum(res_year)
    caco3_precipitation = bio_year_quotas*1000/year_days #mM
    caco3_dissolution = res_year_quotas*1000/year_days
    ca_flux = caco3_dissolution - caco3_precipitation
    ca_array = np.array(ca_flux)/2.5*2

    alk_array = np.array(alk_surface_year['B_C_Alk'])
    alkflux_bottom_year = alkflux_bottom[strt:stp].reset_index()

    calpart = np.zeros(365)
    day = 0
    last_entry = 0
    for month, increment in zip(year_days, ca_array):
        temp = np.linspace(last_entry+increment, last_entry+increment*month, num=month)
        calpart[day:day+month] = temp
        last_entry = temp[-1]
        day += month

    result_array = alk_array + calpart

    caco3_dis = np.zeros(365)
    day = 0
    for month, increment in zip(year_days, caco3_dissolution):
        caco3_dis[day:day+month] = increment
        day += month

    caco3_pre = np.zeros(365)
    day = 0
    for month, increment in zip(year_days, caco3_precipitation):
        caco3_pre[day:day+month] = increment
        day += month

    #gs = gs.GridSpec(2, 2)
    #fig = plt.figure(figsize=cm2inch(18, 14))
    #ax1 = fig.add_subplot(gs[0,0])
    #ax2 = fig.add_subplot(gs[0,1])
    #ax3 = fig.add_subplot(gs[1,0:])
    
    fig = plt.figure(figsize=cm2inch(18, 7))
    ax1 = fig.add_subplot(1, 2, 1) # row-col-num
    ax2 = fig.add_subplot(1, 2, 2) # row-col-num

    ax1.xaxis_date()
    ax2.xaxis_date()
    #ax3.xaxis_date()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    #ax3.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    ax1.plot(alk_surface_year['time'], caco3_dis*2, label = 'CaCO$_3$ dissolution')
    ax1.plot(alk_surface_year['time'], caco3_pre*2, label = 'CaCO$_3$ precipitation')
    ax1.plot(alk_surface_year['time'], alkflux_bottom_year['B_C_Alk   _flux'], label = 'Modelled TA flux at the SWI')
    ax1.plot(alk_surface_year['time'], caco3_dis*2+alkflux_bottom_year['B_C_Alk   _flux'], linewidth=2, label = r'CaCO$_3$ dissolution + TA flux at the SWI')
    ax1.set_ylabel('TA fluxes, mmol m$^{-2}$ d$^{-1}$', fontsize = fntsz)
    ax1.legend(fontsize = lgndsz, title_fontsize = lgndsz, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", borderaxespad=0)

    ax2.plot(alk_surface_year['time'], calpart-calpart.min(), linewidth=2, label=r'Due to CaCO$_3$ dissolution/precipitation')
    ax2.plot(alk_surface_year['time'], alk_array-alk_array.min(), linewidth=2, label=r'From the model calculations')
    ax2.plot(alk_surface_year['time'], result_array - result_array.min(), linewidth=2, label=r'CaCO$_3$ + model calculations')
    ax2.set_ylabel('Relative TA, mmol m$^{-3}$', fontsize=fntsz)
    ax2.legend(fontsize = lgndsz, title_fontsize = lgndsz, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", borderaxespad=0)

    #ax3.plot(alk_surface_year['time'], alk_array-alk_array.min(), linewidth=2, label=r'From the model calculations')
    #ax3.plot(alk_surface_year['time'], result_array - result_array.min(), linewidth=2, label=r'CaCO$_3$ + model calculations')
    #topicture = result_array + total
    #topicture = topicture - topicture.min()
    #ax3.plot(alk_surface_year['time'], topicture, linewidth=2, label=r'CaCO$_3$ + model calculations + nutrients consumption')
    #ax3.set_ylabel('Relative TA, mmol m$^{-3}$', fontsize=fntsz)
    #ax3.legend(fontsize = lgndsz, title_fontsize = lgndsz, loc="best", borderaxespad=0)

    labels = ('(A) ','(B)')#,'(C) ')
    # --- improve the layout
    for i,axis in enumerate((ax1,ax2)):#,ax3)):
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=fntsz, fontweight='bold', va='top', ha='right')
    fig.tight_layout(pad=0.5)

    plt.savefig('Figure7-new.png', dpi = 300)
    #plt.show()

if __name__ == "__main__":
    #plot_intro()
    plot_alkalinity_flux_low_high()
    #plot_alkalinity_flux_sulfur_oxidation()
    plot_alkalinity_flux_porosities1_2_3()
    plot_alk_sulfur_fluxes()
    plot_caco3()
