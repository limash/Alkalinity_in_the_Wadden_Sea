'''
Created on 28. jun. 2017
@author: ELP
'''
import sys,os
from netCDF4 import Dataset,num2date,date2num,date2index
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import matplotlib.dates as mdates
import datetime 
font = {'size' : 14}
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas)
from dateutil.relativedelta import relativedelta

def save_to_dir(dir_name):
    #dir_name = 'Results'
    script_dir = os.path.abspath(os.path.dirname(__file__))
    dir_to_save = os.path.abspath(os.path.join(script_dir,dir_name))
        
    if not os.path.isdir(dir_to_save):
        os.makedirs(dir_to_save)
    filename = '{}\ice_brom_{}.pdf'.format(dir_to_save,name)       
    plt.savefig(filename, format='pdf',transparent = True)
    #plt.savefig(filename, format='pdf', dpi=300,transparent = True)

def plot_3fig(start_year,stop_year,name,title,path): 

    plt.clf() 
    directory = os.getcwd()     
    water_fname = os.path.join(directory,path+'/water.nc')
    sediments_fname = os.path.join(directory,path+'/sediments.nc')

    figure = plt.figure(figsize=(8, 6), dpi=100,
                    facecolor='None',edgecolor='None')  

    fh_water =  Dataset(water_fname)
    fh_sediments =  Dataset(sediments_fname) 

    depth_water = np.array(fh_water.variables['z_faces'][:])
    depth_sed = fh_sediments.variables['z_faces'][:] 

    min_water = np.amin(depth_water)
    max_water = np.amax(depth_water)

    min_sed = np.amin(depth_sed)
    max_sed = np.amax(depth_sed)

    time = fh_water.variables['time']      
    time2 = fh_water.variables['time'][:]

    time_units = fh_water.variables['time'].units

    format_time = num2date(time2,units = time_units,calendar= 'standard')
                
    #assert start_year < stop_year
    to_start = datetime.datetime(start_year,1,1,12,0)
    to_stop= datetime.datetime(stop_year,1,1,12,0)

    start = date2index(to_start, time,#units = time_units,
                        calendar=None, select='nearest')
    stop = date2index(to_stop, time,#units = time_units,
                        calendar=None, select='nearest')

    var_water = np.array(
        fh_water.variables[name][:]).T 
    var_sediments = np.array(
        fh_sediments.variables[name][:]).T                
    data_units = fh_water.variables[name].units

    start_f = num2date(time2[start],units = time_units) 
    stop_f = num2date(time2[stop],units = time_units)    
                                            
    X_water,Y_water = np.meshgrid(time2[start:stop],depth_water)
    X_water = num2date(X_water,units = time_units)   

    X_sed, Y_sed = np.meshgrid(time2[start:stop],depth_sed)
    X_sed = num2date(X_sed,units = time_units)   

    fh_water.close()
    fh_sediments.close()          

    gs = gridspec.GridSpec(2, 1)
    gs.update(left=0.09, right= 1,top = 0.95,bottom = 0.06,
                        wspace=0.1,hspace=0.05)

    ax1 = figure.add_subplot(gs[0]) # water
    ax2 = figure.add_subplot(gs[1]) # sed

    #specify colormap
    #cmap = plt.cm.terrain #'plasma' #'terrain'        
    cmap = plt.get_cmap('viridis') 
    cmap_water = plt.get_cmap('CMRmap') 

    #min = ma.min(var_water)
    #max = ma.max(var_water)
    #var_levels = np.linspace(min,max,num = 20 )

    #plot 2d figures 
    #without interpolation 

    CS4 = ax1.pcolor(X_water,Y_water,var_water[:,start:stop],
                        cmap = cmap_water)
    CS7 = ax2.pcolor(X_sed,Y_sed,var_sediments[:,start:stop], cmap = cmap_water)  
                    
    ax1.set_title(title +' ['+ str(data_units)+']')

    import matplotlib.ticker as ticker

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    # interpolate data to plot 
    #CS1 = ax0.contourf(X,Y, var_ice[:,start:stop],
    #                   cmap = cmap,levels = var_levels)        
    #CS4 = ax1.contourf(X_water,Y_water, var_water[:,start:stop],
    #                   cmap = cmap) 
    #CS7 = ax2.contourf(X_sed,Y_sed, var_sed[:,start:stop],
    #                   cmap = cmap)    

    ax2.axhline(max_water, color='w', linestyle = '--',linewidth = 1 ) 
    ax2.annotate('  Sediment Water Interface',
                xy =(start_f,max_water),
                xytext=(start_f,max_water-0.01),color = 'w')
        
    #ax2.axhline(max_water, color='k', linestyle = '--',linewidth = 1 ) 
    #print (start_f,max_water)
    #ax2.annotate('Sediment Water Interface',
    #            xy =(start_f,max_water),
    #            xytext=(start_f,max_water-0.01),color = 'k')       
    ### Time ticks ### 

    if (stop-start)>= 367:
        dt =  int((stop - start)/365) #number of years
        time_ticks = []
        for n in range(0,dt+1):
            time_ticks.append(
                format_time[start]+relativedelta(years = n))

    def add_colorbar(CS,axis,ma1):
        if ma1 > 10000 or ma1 < 0.001:
            cb = plt.colorbar(CS,ax = axis,pad=0.02,
                    aspect = 7,shrink = 0.9,format=ticker.FuncFormatter(fmt)) 
        else: 
            cb = plt.colorbar(CS,ax = axis,pad=0.01,
                    aspect = 7,shrink = 0.9) 
        return cb

    ma1 = ma.max(var_water[:,start:stop])        
    cb1 = add_colorbar(CS4,ax1,ma1)
    cb2 = add_colorbar(CS7,ax2,ma1)

    #letters = ['(a)','(b)','(c)']
    labels = ["Depth m","Depth m" ]

    n = 0
    for axis in (ax1,ax2): 
        try:
            axis.set_xticks(time_ticks)
        except: NameError
        axis.yaxis.set_label_coords(-0.08, 0.5)
        #axis.text(-0.21, 0.9, letters[n], transform=axis.transAxes , 
        #        size= fontsize) #, weight='bold')
        axis.set_ylabel(labels[n], fontsize = 11)
            
        n=n+1
        #plt.tick_params(axis='both', which='major', labelsize= fontsize) 
        
    ax1.set_ylim(max_water,min_water)
    ax2.set_ylim(max_sed,min_sed)  #

    # hide horizontal axis labels 
    ax1.set_xticklabels([]) 

    #plt.yticks(fontsize=fontsize, rotation=90) 
    #ax.yaxis.label.set_size(16) 
    #plt.rcParams.update({'font.size': 14})
    if (stop-start) > 367 and (stop-start) < 365*6:
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter("%b '%y"))            
    elif (stop-start)>= 365*6:
        #ax5.xaxis.set_major_formatter(
        #     mdates.DateFormatter("%b '%y")) 
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter('%Y')) 
        #ticks = np.arange(time[start:stop],time[start:stop],50)

        #ax2.xaxis.set_major_formatter(
        #    mdates.DateFormatter('%m/%Y'))   
    else : 
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter('%b'))
    plt.show()