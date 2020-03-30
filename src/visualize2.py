'''
Created on 28. jun. 2017
original @author: ELP
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import datetime
from dateutil.relativedelta import relativedelta
from matplotlib import gridspec
from pandas.plotting import register_matplotlib_converters
from netCDF4 import Dataset, num2date, date2index

sns.set()
register_matplotlib_converters()
font = {'size': 11}


def plot_fig(start_year, stop_year, name, title, path):

    plt.clf()
    directory = os.getcwd()
    water_fname = os.path.join(directory, path+'/water.nc')
    sediments_fname = os.path.join(directory, path+'/sediments.nc')

    figure = plt.figure(figsize=(8, 6), dpi=100,
                        facecolor='None', edgecolor='None')

    fh_water = Dataset(water_fname)
    fh_sediments = Dataset(sediments_fname)

    depth_water = np.array(fh_water.variables['z_faces'][:])
    depth_sed = fh_sediments.variables['z_faces'][:]
    min_water = np.amin(depth_water)
    max_water = np.amax(depth_water)
    min_sed = np.amin(depth_sed)
    max_sed = np.amax(depth_sed)

    time = fh_water.variables['time']
    time2 = fh_water.variables['time'][:]
    time_units = fh_water.variables['time'].units
    format_time = num2date(time2, units=time_units, calendar='standard')

    to_start = datetime.datetime(start_year, 1, 1, 12, 0)
    to_stop = datetime.datetime(stop_year, 1, 1, 12, 0)

    start = date2index(to_start, time, calendar=None, select='nearest')
    stop = date2index(to_stop, time, calendar=None, select='nearest')

    data_units = fh_water.variables[name].units

    start_f = num2date(time2[start], units=time_units)

    X_water, Y_water = np.meshgrid(time2[start:stop], depth_water)
    X_water = num2date(X_water, units=time_units)

    X_sed, Y_sed = np.meshgrid(time2[start:stop], depth_sed)
    X_sed = num2date(X_sed, units=time_units)

    fh_water.close()
    fh_sediments.close()

    gs = gridspec.GridSpec(2, 1)
    gs.update(left=0.09, right=1, top=0.95,
              bottom=0.06, wspace=0.1, hspace=0.05)

    ax1 = figure.add_subplot(gs[0])
    ax2 = figure.add_subplot(gs[1])

    ax1.set_title(title + ' [' + str(data_units) + ']')

    import matplotlib.ticker as ticker

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    ax2.axhline(max_water, color='w', linestyle='--', linewidth=1)
    ax2.annotate('  Sediment Water Interface', xy=(start_f, max_water),
                 xytext=(start_f, max_water-0.01), color='w')

    if (stop-start) >= 367:
        dt = int((stop - start)/365)
        time_ticks = []
        for n in range(0, dt+1):
            time_ticks.append(
                format_time[start]+relativedelta(years=n))

    def add_colorbar(CS, axis, ma1):
        if ma1 > 10000 or ma1 < 0.001:
            cb = plt.colorbar(CS, ax=axis, pad=0.02, aspect=7,
                              shrink=0.9, format=ticker.FuncFormatter(fmt))
        else:
            cb = plt.colorbar(CS, ax=axis, pad=0.01,
                              aspect=7, shrink=0.9)
        return cb

    labels = ["Depth m", "Depth m"]

    n = 0
    for axis in (ax1, ax2):
        axis.set_xticks(time_ticks)
        axis.yaxis.set_label_coords(-0.08, 0.5)
        axis.set_ylabel(labels[n], fontsize=11)
        n = n+1

    ax1.set_ylim(max_water, min_water)
    ax2.set_ylim(max_sed, min_sed)
    ax1.set_xticklabels([])
    if (stop-start) > 367 and (stop-start) < 365*6:
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter("%b '%y"))
    elif (stop-start) >= 365*6:
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter('%Y'))
    else:
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter('%b'))
    plt.show()
