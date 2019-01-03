import numpy as np # imports a fast numerical programming library
import matplotlib.pyplot as plt #sets up plotting under plt
import pandas as pd #lets us handle data as dataframes
from datetime import datetime
from matplotlib import rc
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
    ax.set_title('Figure 1. Surface water TA measured and calculated from the alkalinity-salinity relation (Millero 1997),\n FerryBox measurements (54.19$^\circ$ N, 6.99$^\circ$ E, 2017)', fontsize=16, fontweight='bold')
    plt.xlabel("Day in a year")
    plt.ylabel('Total Alkalinity, $ \mu M$')

def plot_intro():
    north7 = pd.read_csv("HafniaDataNorth7Shamil.csv")
    north7 = treatbiogeodata(north7)
    fig = plotTA(north7)