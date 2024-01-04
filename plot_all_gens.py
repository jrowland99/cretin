#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:56:57 2023

@author: jeff
"""
import numpy as np
import re
import os
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import statistics
import sys
from labellines import labelLines 
sys.path.append("/home/jeff/Research/Code_Projects/git/main/cloudy/")
from NeonCellSS import *
from NeonCell import *
from readPS import *


plt.rcParams["figure.dpi"] = 650
plt.rcParams["lines.linewidth"] = 1
plt.rcParams.update({"font.size": 10})
plt.rcParams.update({"legend.fontsize": 11})
plt.rcParams["figure.figsize"] = (12, 12)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)

density_str_dict={30:'$n_i = 10^{18}cm^{-3}$',15:'$n_i =5 \cdot 10^{17}cm^{-3}$',7.5:'$n_i = 2.5 \cdot 10^{17}cm^{-3}$'}
gen_str_dict = {'gen1c':'Generation 1, Close \n','gen2c':'Generation 2, Close \n','gen2f':'Generation 2, Far \n'}
BWA_flux = {'gen2c':1.009238e+12*1e7*0.80 ,'gen2f':4.520144e+11*1e7*0.80,'gen1c':7.443568e+11*1e7}
gplc_ws = [5.824228e+007,7.130530e+008,2.891976e+009,5.740544e+009,1.096900e+010,6.269020e+009,4.346340e+009,1.429054e+009,7.130530e+008,3.133853e+008,2.610261e+008]

def fetch_helios_temp(pressure):
    helios_time=[70,75,80,85,90,91,92,93,94,95,95.5,96,96.5,97,97.5,98,98.5,99,99.5,100,100.5,
                     101,101.5,102,102.5,103,103.5,104,104.5,105,106,107,108,109,110]
    if pressure == 7.5:   
        helios_ne_T = [2.646215,4.6113891,7.2893005,11.064016,15.448457,16.302372,17.124757,17.87929,18.60159,19.352553,19.808579,20.305306,20.860994,21.372897,21.773992,22.03143,22.170141,22.273682,22.37318,22.486271,22.608346,22.720299,22.826706,22.927493,23.009448,23.063882,23.128089,23.175455,23.215723,23.241685,23.331713,23.365972,23.42532,23.486416,23.543447]
        helios_ne_cool=[1.3137461,1.7730058,10.293328,68.026629,297.47736,375.27991,484.55859,641.33442,880.91361,1297.9348,1549.8419,1703.9586,1743.5785,1580.4132,1282.7768,970.21323,730.62485,652.35252,645.91101,657.0936,583.47138,502.24371,394.12879,280.01064,246.87938,218.08516,191.63859,166.88582,162.59843,159.17595,163.06919,173.34566,181.51599,186.57383,190.30638]
        helios_ne_heat=[38.617488,64.915095,129.45141,321.91172,566.61658,627.51151,721.95999,866.64786,1090.6413,1541.3488,1826.756,2017.4916,2061.4346,1857.7678,1482.2958,1082.2916,797.8278,711.96107,708.66463,721.66369,652.74571,568.43829,454.48742,332.16007,287.36111,250.71072,218.76507,190.53205,183.09227,178.10052,180.19988,189.55544,197.42372,202.31258,205.8285]
        helios_to_cloudy = 8.3675e13/(4*np.pi)
    elif pressure == 15:   
        helios_ne_T = [2.3766957,4.1598528,6.5820136,10.238093,15.020088,15.995314,16.904061,17.803864,18.690759,19.676253,20.27489,20.975236,21.722121,22.380288,22.928031,23.299639,23.495895,23.628329,23.776315,23.899768,24.011263,24.130996,24.23529,24.323765,24.380377,24.43559,24.490903,24.522264,24.589062,24.632186,24.762276,24.936785,25.114917,25.343334,25.58842]
        helios_ne_cool=[0.85481175,1.1921946,8.7183183,62.900471,296.6658,387.1497,519.85245,684.52775,911.84551,1287.2053,1532.8135,1619.8115,1604.3561,1462.3431,1241.7922,1038.1756,846.3176,774.76729,769.28444,778.08045,709.20011,629.33268,526.49607,419.20127,390.51664,365.32558,343.48906,324.51808,325.08789,326.27462,338.08705,352.98342,362.04355,364.8563,364.14134]
        helios_ne_heat=[30.850082,56.391139,104.15659,283.85743,584.4902,665.31147,782.53431,944.19647,1174.7022,1623.4723,1922.9339,2066.0725,2059.3286,1852.7946,1521.3993,1203.6467,945.45642,857.62507,849.50733,853.82978,782.12434,694.11907,582.09466,463.90396,425.39239,394.96014,370.98117,352.20536,354.33124,357.83345,374.91889,394.29856,404.90632,407.32306,404.96236]
        helios_to_cloudy = 1.674e14/(4*np.pi)
    elif pressure == 30:   
        helios_ne_T = [1.9003743,3.3420202,5.2529009,8.4383343,13.303585,14.467487,15.611836,16.774004,17.9219,19.231879,20.051885,20.971468,21.991124,23.085527,23.99607,24.720297,25.294825,25.696213,26.038329,26.423326,26.731916,27.069578,27.365693,27.600977,27.887887,28.109172,28.369847,28.62273,28.946274,29.244919,29.923936,30.590454,31.160992,31.707974,32.166729]
        helios_ne_cool=[0.54479415,0.74278274,5.7437449,41.633006,227.46188,316.3708,431.44551,614.77055,859.24846,1231.1612,1508.9884,1667.8502,1775.7088,1633.7349,1404.9319,1251.5054,1087.4656,1001.5815,980.72847,973.57186,901.28578,820.76547,716.47327,605.48675,572.52039,540.01546,511.00022,486.31451,477.69446,470.76172,462.46275,455.57667,442.37474,429.45329,418.54966]
        helios_ne_heat=[22.777038,38.928204,71.225292,203.87302,512.23297,619.29036,752.30756,949.96951,1208.7384,1694.5201,2066.3411,2314.0135,2452.0891,2245.7193,1882.3852,1570.2093,1280.4512,1139.5155,1096.8451,1074.5192,991.55718,900.57837,785.79556,667.4018,630.36755,596.97774,569.07787,547.10901,541.02714,535.93583,528.41763,517.54936,496.56781,476.16216,459.63413]
        helios_to_cloudy = 3.347e14/(4*np.pi) # This factor may need to change to a solid angle/4pi
        #helios_to_cloudy = 3.347e14
    helios_ne_T    = np.asarray(helios_ne_T)
    helios_ne_cool = np.asarray(helios_ne_cool)*helios_to_cloudy
    helios_ne_heat = np.asarray(helios_ne_heat)*helios_to_cloudy
    hdict = {'time':helios_time,'te':helios_ne_T,'heat':helios_ne_heat,'cool':helios_ne_cool}
    hdf = pd.DataFrame(hdict)
    return hdf

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
def plot_rectangle(x, y, xerr, yerrs, fig, ax,color):
    # x and y are currently the center of the box, xerr is symmetric, yerrs = [ylow,yhigh]
    a = x - xerr
    b = y -float(yerrs[0][0]) # y error is asymmetric
    width = 2*xerr
    height = float(yerrs[0][0])+float(yerrs[1][0])
    rect = patches.Rectangle((a, b), width, height, linewidth=1,linestyle='none', facecolor=color,alpha=0.2)
    ax.add_patch(rect)
def plot_exp(gen,pressure,fig,ax):
    xerr = 1.5
    color = 'tab:blue'
    plotOn = True
    offset = 0
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
    offset = 0
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

def get_cretin_temp(path):
    with open(path,'r') as f:
        chunks = f.read().split('#')
        chunks = chunks[1:]
        for chunk in chunks:
            lines = chunk.split('\n')
            lines = [line for line in lines if not line.startswith('$')]
            header = lines[0]
            lines = lines[1:]
            if header.strip() == 'TEV, TIV vs TIME':
                time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
                tev   = [float(line.split()[1]) for line in lines if not len(line)==0]
                tiv   = [float(line.split()[2]) for line in lines if not len(line)==0]
    
        results = {'time':time,'te':tev,'ti':tiv}
    return results
def get_cretin_eden(path):
    with open(path,'r') as f:
        chunks = f.read().split('#')
        chunks = chunks[1:]
        for chunk in chunks:
            lines = chunk.split('\n')
            lines = [line for line in lines if not line.startswith('$')]
            header = lines[0]
            lines = lines[1:]
            if header.strip() == 'NE vs TIME':
                time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
                ne   = [float(line.split()[1]) for line in lines if not len(line)==0]
    
        results = {'time':time,'ne':ne}
    return results
def get_cretin_erad(path):
    with open(path,'r') as f:
        chunks = f.read().split('#')
        chunks = chunks[1:]
        for chunk in chunks:
            lines = chunk.split('\n')
            lines = [line for line in lines if not line.startswith('$')]
            header = lines[0]
            lines = lines[1:]
            if header.strip() == 'ERAD vs TIME':
                time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
                erad  = [float(line.split()[1]) for line in lines if not len(line)==0]
    
        results = {'time':time,'erad':erad}
    return results
gens = ['gen1c','gen2c','gen2f']
#%% Fetch Cloudy Time-Dependent
ctd_dir = '/home/jeff/Research/Code_Projects/Cloudy_Projects/Z_machine/Neon_cell/Time_dep/td_weighted_avg/'
ctdtimes = [90,92,94,96,98,100]
pressures = [7.5,15,30]
ctd_gen_dict = {}
gens = ['gen1c','gen2c','gen2f']
for gen in gens:
    directory = os.path.join(ctd_dir,gen)
    ctd_dict = {}
    for pressure in pressures:# add in the name for each td avg
        path = os.path.join(directory,str(pressure)+'torr/')
        density = 1e18/(30/pressure)
        td_dict={}
        stats = {}
        for time in ctdtimes:
            td = NeonCell.from_file(path=path+'{}/Ne_{}ns'.format(time,time),name="Ne_{}ns".format(str(time)),density=density,pressure=pressure,cumul_cont_flag=False,read_spectrum_flag=False)
            # Skip BWA for now until we have it clear for all generations
            
            # trim_eden_avg =trim_value(td.times,td.eden)
            # BWA = backlighter_weight_avg(gplc_ws,trim_eden_avg)
            # ioniz_param = 4*np.pi*BWA_flux['gen1c']/BWA
            # td.xi = ioniz_param
            # print('P = {},t = {}, xi = {}'.format(pressure,time,round(ioniz_param,2)))
            td_dict[time]=td
            stats[time]=td.temps # do this so we can average the temperatures for each SED run
        stats= pd.DataFrame.from_dict(stats)
        savedict = {'time':td.times,'te':stats.mean(axis=1),'te_std':stats.std(axis=1)}
        ctd_dict[pressure] = pd.DataFrame(savedict)
    ctd_gen_dict[gen]= ctd_dict
#%% Get cretin temps
cretin_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/all_gen_v2/'
pressures = [7.5,15,30]
cretin_fname = 'ne_cell.plt'
cretin_gen_dict = {}
for gen in gens:
    directory = os.path.join(cretin_dir,gen)
    cretin_p_dict = {}
    for pressure in pressures:# add in the name for each td avg
        path = os.path.join(directory,str(pressure)+'torr/')
        print(os.path.join(path,cretin_fname))
        label = str(gen)+' '+str(pressure)+' torr'
        cretin_p_dict[pressure] = get_cretin_temp(os.path.join(path,cretin_fname))
        plt.plot(cretin_p_dict[pressure]['time'],cretin_p_dict[pressure]['te'], label = label)
        plt.grid(linestyle=':', linewidth=0.5)
        plt.legend()
    cretin_gen_dict[gen]=cretin_p_dict
#%%
# cretin_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/all_gen_v2/'
cretin_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/all_gen_1zone_v2/'
pressures = [7.5,15,30]
cretin_fname = 'ne_cell.plt'
cretin_gen_dict = {}
for gen in gens:
    directory = os.path.join(cretin_dir,gen)
    cretin_p_dict = {}
    fig,ax = plt.subplots()
    for pressure in pressures:# add in the name for each td avg
        path = os.path.join(directory,str(pressure)+'torr/')
        # print(os.path.join(path,cretin_fname))
        label = str(gen)+' '+str(pressure)+' torr'
        cretin_p_dict[pressure] = get_cretin_eden(os.path.join(path,cretin_fname))
        plt.plot(cretin_p_dict[pressure]['time'],cretin_p_dict[pressure]['ne'], label = label)
        plt.grid(linestyle=':', linewidth=0.5)
        plt.ylabel('Electron density $cm^{-3}$')
        plt.xlabel('Time [ns]')
        plt.legend()
    cretin_gen_dict[gen]=cretin_p_dict
    plt.savefig(os.path.join(cretin_dir,f'{gen}_eden'+ext), transparent=False)
#%%
cretin_gen_dict = {}
for gen in gens:
    directory = os.path.join(cretin_dir,gen)
    cretin_p_dict = {}
    fig,ax = plt.subplots()
    for pressure in pressures:# add in the name for each td avg
        path = os.path.join(directory,str(pressure)+'torr/')
        # print(os.path.join(path,cretin_fname))
        label = str(gen)+' '+str(pressure)+' torr'
        cretin_p_dict[pressure] = get_cretin_temp(os.path.join(path,cretin_fname))
        plt.plot(cretin_p_dict[pressure]['time'],cretin_p_dict[pressure]['te'], label = label)
        plt.grid(linestyle=':', linewidth=0.5)
        plt.ylabel('Temperature [eV]')
        plt.xlabel('Time [ns]')
        plt.legend()
    cretin_gen_dict[gen]=cretin_p_dict
    plt.savefig(os.path.join(cretin_dir,f'{gen}_temperature'+ext), transparent=False)
#%%
cretin_gen_dict = {}
for gen in gens:
    directory = os.path.join(cretin_dir,gen)
    cretin_p_dict = {}
    # fig,ax = plt.subplots()
    for pressure in pressures:# add in the name for each td avg
        path = os.path.join(directory,str(pressure)+'torr/')
        # print(os.path.join(path,cretin_fname))
        label = str(gen)+' '+str(pressure)+' torr'
        cretin_p_dict[pressure] = get_cretin_erad(os.path.join(path,cretin_fname))
        plt.plot(cretin_p_dict[pressure]['time'],cretin_p_dict[pressure]['erad'], label = label)
        plt.grid(linestyle=':', linewidth=0.5)
        plt.ylabel('Radiation Energy Density ergs/$cm^3$')
        plt.xlabel('Time [ns]')
        plt.yscale('log')
        plt.legend()
    cretin_gen_dict[gen]=cretin_p_dict
    plt.savefig(os.path.join(cretin_dir,f'{gen}_erad'+ext), transparent=False)    