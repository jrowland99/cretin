#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 17:54:51 2023

@author: jeff
"""

"""
Created on Fri Jul 21 11:12:37 2023

@author: jeff

Create complete 9-panel temperature plot for neon cell
include cretin, cloudy TD, cloudy SS, Helios

"""
import numpy as np
import os
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/jeff/Research/Code_Projects/git/main/cretin/")
sys.path.append("/home/jeff/Research/Code_Projects/git/main/cloudy/")
sys.path.append("/home/jeff/Research/Code_Projects/git/main/prism/")
import matplotlib.patches as patches
import statistics
from NeonCellSS import *
from NeonCell import *
from readPS import *
from CretinData_old import *
plt.rcParams["figure.dpi"] = 650
plt.rcParams["lines.linewidth"] = 1
plt.rcParams.update({"font.size": 10})
plt.rcParams.update({"legend.fontsize": 11})
plt.rcParams["figure.figsize"] = (8, 3)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)

density_str_dict={30:'$n_i = 10^{18}cm^{-3}$',15:'$n_i =5 \cdot 10^{17}cm^{-3}$',7.5:'$n_i = 2.5 \cdot 10^{17}cm^{-3}$'}
gen_str_dict = {'gen1c':'Generation 1, Close \n','gen2c':'Generation 2, Close \n','gen2f':'Generation 2, Far \n'}
gens = ['gen1c','gen2c','gen2f']
press = ['7-5_torr','15_torr','30_torr']
pdict = {'7-5_torr':7.5,'15_torr':15,'30_torr':30}
BWA_flux = {'gen2c':1.009238e+12*1e7*0.80 ,'gen2f':4.520144e+11*1e7*0.80,'gen1c':7.443568e+11*1e7}
gplc_ws = [5.824228e+007,7.130530e+008,2.891976e+009,5.740544e+009,1.096900e+010,6.269020e+009,4.346340e+009,1.429054e+009,7.130530e+008,3.133853e+008,2.610261e+008]
def trim_value(times,values): #Averages the values array to match the size of the weights array
    time_range = np.linspace(start=95,stop=105,endpoint=True,num=11)
    value_df = pd.DataFrame({'time':times,'value':values})
    value_df = value_df[value_df['time'].ge(time_range[0])]
    value_df = value_df[value_df['time'].le(time_range[-1])]
    new_values = []
    for t in time_range:
        avg_array =[]
        for index,row in value_df.iterrows():
            if round(row['time']) == t:
                avg_array.append(row['value'])
        avg = np.average(avg_array)
        new_values.append(avg)
    return new_values
 

def backlighter_weight_avg(weights, values):
    assert len(weights)==len(values)
    return np.dot(weights,values)/np.sum(weights)
#--------------------------------------------------------------------------
#
#
#
#
#
#--------------------------------------------------------------------------
def plot_rectangle(x, y, xerr, yerrs, fig, ax,color):
    # x and y are currently the center of the box, xerr is symmetric, yerrs = [ylow,yhigh]
    a = x - xerr
    b = y -float(yerrs[0][0]) # y error is asymmetric
    width = 2*xerr
    height = float(yerrs[0][0])+float(yerrs[1][0])
    rect = patches.Rectangle((a, b), width, height, linewidth=1,linestyle='none', facecolor=color,alpha=0.1)
    ax.add_patch(rect)
def plot_exp(gen,pressure,fig,ax):
    xerr = 1.5
    color = 'tab:blue'
    plotOn = True
    offset = 100
    #these data are from D. Mayes' dissertation
    if gen =='gen1c':
        if pressure == 7.5:# 17.9 +4.7/-5.2 eV
            lower_error =  5.2
            upper_error =  4.7
            te = 17.9
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(100-offset,te,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$17.9^{+4.7}_{-5.2} eV$')
            
        elif pressure == 15:# 
            lower_error =  6.6
            upper_error =  6.5
            te = 23.6
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(100-offset,23.6,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$23.6^{+6.5}_{-6.6} eV$')
            
        elif pressure == 30:# 26 +-5 eV 
            asymmetric_error = [[4.2],[4.2]]
            te = 21.3
            if plotOn:
                ax.errorbar(100-offset,21.3,yerr=4.2,xerr=xerr,capsize=3, label='$21.3 \pm 4.2 eV$')
    if gen =='gen2c':
        if pressure == 7.5:# 25 +6./-6 eV
            lower_error =  7.1
            upper_error =  6.4
            te = 25
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(100-offset,25,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$25^{+6.4}_{-7.1} eV$')
            
        elif pressure == 15:# 24 +4/-5 eV
            lower_error =  8
            upper_error =  7.8
            te=26.9
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(100-offset,26.9,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$26.9^{+7.8}_{-8.0} eV$')
            
        elif pressure == 30:# 26 +-5 eV 
            asymmetric_error = [[4.7],[4.7]]
            te = 24
            if plotOn:
                ax.errorbar(100-offset,24,yerr=4.7,xerr=xerr,capsize=3, label='$24 \pm 4.7 eV$')
    if gen =='gen2f':
        if pressure == 7.5:# 24 +4/-6 eV
            lower_error =  3.8
            upper_error =  3.1
            te = 17.7
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(100-offset,17.7,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$17.7^{+3.1}_{-3.8} eV$')
            
        elif pressure == 15:# 24 +4/-5 eV
            lower_error =  6.8
            upper_error =  6.6
            te = 29.2
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(100-offset,29.2,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$29.2^{+6.6}_{-6.8} eV$')
            
        elif pressure == 30:# 26 +-5 eV
            te = 23.9
            asymmetric_error = [[10.2],[10.2]]
            if plotOn:
                ax.errorbar(100-offset,23.9,yerr=10.2,xerr=xerr,capsize=3, label='$23.9 \pm 10.2 eV$')

    plot_rectangle(100-offset,te, xerr, asymmetric_error, fig, ax,color)
def plot_2020_exp(gen,pressure,fig,ax):
    xerr = 1.5
    color = 'darkviolet'
    plotOn = True
    offset = 100
    if gen !='gen1c':
        return
    if pressure == 7.5:# 24 +4/-6 eV
        lower_error =  6
        upper_error =  4
        asymmetric_error = [[lower_error], [upper_error]]
        te = 24
        if plotOn:
            ax.errorbar(100-offset,24,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$24^{+4}_{-6} eV$')
        
    elif pressure == 15:# 24 +4/-5 eV
        lower_error =  5
        upper_error =  4
        asymmetric_error = [[lower_error], [upper_error]]
        te = 24
        if plotOn:
            ax.errorbar(100-offset,24,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$24^{+4}_{-5} eV$')
        
    elif pressure == 30:# 26 +-5 eV  
        lower_error =  5
        upper_error =  5
        asymmetric_error = [[lower_error], [upper_error]]
        te = 26
        if plotOn:
            ax.errorbar(100-offset,26,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$26 \pm 5 eV$')
        
    plot_rectangle(100-offset,te, xerr, asymmetric_error, fig, ax,color)
#--------------------------------------------------------------------------
#
#
#
#
#
#--------------------------------------------------------------------------
# #%% Fetch Cloudy Time-Dependent
# ctd_dir = '/home/jeff/Research/Code_Projects/Cloudy_Projects/Z_machine/Neon_cell/Time_dep/all_gen_all_SED/'
# ctdtimes = [80,90,92,94,96,98,100,102,104]
# pressures = [7.5,15,30]
# ctd_gen_dict = {}
# for gen in gens:
#     directory = os.path.join(ctd_dir,gen)
#     ctd_dict = {}
#     for pressure in pressures:# add in the name for each td avg
#         path = os.path.join(directory,str(pressure)+'torr/')
#         density = 1e18/(30/pressure)
#         td_dict={}
#         stats = {}
#         for time in ctdtimes:
#             td = NeonCell.from_file(path=path+'{}/Ne_{}ns'.format(time,time),name="Ne_{}ns".format(str(time)),density=density,pressure=pressure,cumul_cont_flag=False,read_spectrum_flag=False)
#             # Skip BWA for now until we have it clear for all generations
            
#             # trim_eden_avg =trim_value(td.times,td.eden)
#             # BWA = backlighter_weight_avg(gplc_ws,trim_eden_avg)
#             # ioniz_param = 4*np.pi*BWA_flux['gen1c']/BWA
#             # td.xi = ioniz_param
#             # print('P = {},t = {}, xi = {}'.format(pressure,time,round(ioniz_param,2)))
#             td_dict[time]       = td
#             stats[time]         = td.temps # do this so we can average the temperatures for each SED run

#         stats               = pd.DataFrame.from_dict(stats)
#         savedict            = {'time':td.times,'te':stats.mean(axis=1),'te_std':stats.std(axis=1)}
#         ctd_dict[pressure]  = pd.DataFrame(savedict)
#     ctd_gen_dict[gen]       = ctd_dict
#%% Plot Temperatures for all
gens = ['gen1c']
directory='/home/jeff/Research/Code_Projects/Cloudy_Projects/Z_machine/Neon_cell/Ne_cell_steady_state/all_gen/'
cretin_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/sw36_off/'
cretin_dir2 = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/steady_state_latest/sw36_off/'
helios_dir = '/home/jeff/Research/Code_Projects/Cloudy_Projects/Z_machine/Neon_cell/helios_data'
fig, axes = plt.subplots(1,3, sharey='row')
plotPeak = True
offset = 100
for i, gen in enumerate(gens):
    # ctdp  = ctd_gen_dict[gen] #cloudy time dependent pressure
    for j, p in enumerate(press):
        pressure = pdict[p]
        ax = axes[j]
        # #--------------------------------------------------------------------------
        # # Cloudy Steady State
        # #--------------------------------------------------------------------------
        # file = gen+'_'+p+'_thermal.csv'
        # cloudySSpath = os.path.join(directory,file)
        # cssdf = pd.read_csv(cloudySSpath)  
        # cssdf.columns = ['time','temp','htot','ctot']
        # ax.plot(cssdf['time']-offset,cssdf['temp'],label = 'Cloudy SS',color='k',linestyle='dotted',marker='x',markersize=2)
        # #--------------------------------------------------------------------------
        # # Cloudy Time-dependent
        # #--------------------------------------------------------------------------
        # ctd_df = ctdp[pressure]
        # ax.plot(ctd_df['time']-offset,ctd_df['te'],label = 'Cloudy TD',color='goldenrod')
        # ax.fill_between(ctd_df['time']-offset,ctd_df['te'],ctd_df['te']-ctd_df['te_std'],alpha=0.3,color='grey', edgecolor='none')
        # ax.fill_between(ctd_df['time']-offset,ctd_df['te'],ctd_df['te']+ctd_df['te_std'],alpha=0.3,color='grey', edgecolor='none')
        # #--------------------------------------------------------------------------
        # # Helios
        # #--------------------------------------------------------------------------
        
        # hfile = f"{gen}_{p}_te_avg.txt"
        # hpath = os.path.join(helios_dir,hfile)
        # hdf = pd.read_csv(hpath)
        # ax.plot(hdf['time']-offset,hdf['te'],label = 'HELIOS',color = 'r')
        #--------------------------------------------------------------------------    
        # Cretin
        #--------------------------------------------------------------------------
        cretin_plt = os.path.join(cretin_dir,gen,str(pressure)+'torr','ne_cell.plt')
        cretin_obj = CretinData(cretin_plt)
        cretin_obj.load_data()
        cretin_data = cretin_obj.get_data()
        ax.plot(cretin_data['time']-offset,cretin_data['te'],label = 'Cretin TD')
        #--------------------------------------------------------------------------
        cretin_plt = os.path.join(cretin_dir2,gen,str(pressure)+'torr','ne_cell.plt')
        cretin_obj = CretinData(cretin_plt)
        cretin_obj.load_data()
        cretin_data = cretin_obj.get_data()
        ax.plot(cretin_data['time']-offset,cretin_data['te'],label = 'Cretin SS')
        # plot_exp(gen,pressure,fig,ax)
        # plot_2020_exp(gen,pressure,fig,ax)
        if plotPeak:
            ax.set_xlim([-10,6])
        ax.set_xlabel('Time [ns]')
        if j==0:
            ax.set_ylabel('Temperature [eV]')
        ax.grid(linestyle=':', linewidth=0.5)
        ax.legend(loc='upper left')
        if j==1:
            titlestr = gen_str_dict[gen]+density_str_dict[pressure]+'('+str(pressure)+' torr)'
        else:
            titlestr = density_str_dict[pressure]+'('+str(pressure)+' torr)'
        ax.set_title(titlestr,fontsize=14)
        #--------------------------------------------------------------------------
        # ax_twin     = ax.twinx()
        # yticks      = ax.get_yticks()
        ylim        = [0,110]
        ax.set_ylim(ylim)
        # ax_twin.set_yticks( yticks)
        # ax_twin.set_ylim(ylim)
        if j == 2:
            # 
            ax_twin     = ax.twinx()
            yticks      = ax.get_yticks()
            yticklabels = ax.get_yticklabels()
            ylim        = ax.get_ylim()
            ax_twin.set_yticks( yticks)
            ax_twin.set_ylim(ylim)
            #ax_twin.set_yticklabels(yticklabels)
            #ax_twin.yaxis.tick_right()
            ax_twin.set_ylabel('Temperature [eV]')
            #ax_twin.tick_params(axis='y')
    plt.suptitle('Cretin w/o Continuum radiation transport')
plt.tight_layout()
# plt.savefig(os.path.join(directory,'all_gen_temperatures_all_time'+ext), transparent=False)
# plt.savefig(os.path.join(directory,'all_gen_temperatures'+ext), transparent=False)
if plotPeak:
    # plt.savefig(os.path.join(cretin_dir,'all_gen_Te_peak'+ext), transparent=False)
    plt.savefig(os.path.join(cretin_dir,'sw36_off_gen1c'+ext), transparent=False)
else:
    plt.savefig(os.path.join(cretin_dir,'all_gen_temperatures'+ext), transparent=False)