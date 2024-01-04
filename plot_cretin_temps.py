#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 17:50:30 2023

@author: jeff
"""
import re
import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import statistics
import os
from labellines import labelLines 
sys.path.append("/home/jeff/Research/Code_Projects/git/main/cloudy/")
import NeonCell  as nc
import NeonCellSS as ncss
import warnings
from scipy.interpolate import interp1d
warnings.filterwarnings("ignore", category=FutureWarning)

plt.rcParams["figure.dpi"] = 250
plt.rcParams["lines.linewidth"] = 1
plt.rcParams.update({"font.size": 10})
plt.rcParams.update({"legend.fontsize": 10})
plt.rcParams["figure.figsize"] = (4, 4)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)
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
path = '/home/jeff/Research/Code_Projects/Cloudy_Projects/Z_machine/Neon_cell/Time_dep/5_17_23/30torr/100'
td_100ns = nc.NeonCell.from_file(path=os.path.join(path,'Ne_100ns'),name='Ne_100ns',density=1e18,pressure=30,cumul_cont_flag=False,read_spectrum_flag=False)
eden = dict(zip(td_100ns.times,td_100ns.eden))

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
    offset = -1
    #these data are from D. Mayes' dissertation
    if gen =='gen1c':
        if pressure == 7.5:# 17.9 +4.7/-5.2 eV
            lower_error =  5.2
            upper_error =  4.7
            te = 17.9
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(99-offset,te,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$17.9^{+4.7}_{-5.2} eV$')
            
        elif pressure == 15:# 
            lower_error =  6.6
            upper_error =  6.5
            te = 23.6
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(99-offset,23.6,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$23.6^{+6.5}_{-6.6} eV$')
            
        elif pressure == 30:# 26 +-5 eV 
            asymmetric_error = [[4.2],[4.2]]
            te = 21.3
            if plotOn:
                ax.errorbar(99-offset,21.3,yerr=4.2,xerr=xerr,capsize=3, label='$21.3 \pm 4.2 eV$')
    if gen =='gen2c':
        if pressure == 7.5:# 25 +6./-6 eV
            lower_error =  7.1
            upper_error =  6.4
            te = 25
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(99-offset,25,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$25^{+6.4}_{-7.1} eV$')
            
        elif pressure == 15:# 24 +4/-5 eV
            lower_error =  8
            upper_error =  7.8
            te=26.9
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(99-offset,26.9,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$26.9^{+7.8}_{-8.0} eV$')
            
        elif pressure == 30:# 26 +-5 eV 
            asymmetric_error = [[4.7],[4.7]]
            te = 24
            if plotOn:
                ax.errorbar(99-offset,24,yerr=4.7,xerr=xerr,capsize=3, label='$24 \pm 4.7 eV$')
    if gen =='gen2f':
        if pressure == 7.5:# 24 +4/-6 eV
            lower_error =  3.8
            upper_error =  3.1
            te = 17.7
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(99-offset,17.7,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$17.7^{+3.1}_{-3.8} eV$')
            
        elif pressure == 15:# 24 +4/-5 eV
            lower_error =  6.8
            upper_error =  6.6
            te = 29.2
            asymmetric_error = [[lower_error], [upper_error]]
            if plotOn:
                ax.errorbar(99-offset,29.2,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$29.2^{+6.6}_{-6.8} eV$')
            
        elif pressure == 30:# 26 +-5 eV
            te = 23.9
            asymmetric_error = [[10.2],[10.2]]
            if plotOn:
                ax.errorbar(99-offset,23.9,yerr=10.2,xerr=xerr,capsize=3, label='$23.9 \pm 10.2 eV$')

    plot_rectangle(99-offset,te, xerr, asymmetric_error, fig, ax,color)
def plot_2020_exp(gen,pressure,fig,ax):
    xerr = 1.5
    color = 'darkviolet'
    plotOn = True
    offset = -1
    if gen !='gen1c':
        return
    if pressure == 7.5:# 24 +4/-6 eV
        lower_error =  6
        upper_error =  4
        asymmetric_error = [[lower_error], [upper_error]]
        te = 24
        if plotOn:
            ax.errorbar(99-offset,24,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$24^{+4}_{-6} eV$')
        
    elif pressure == 15:# 24 +4/-5 eV
        lower_error =  5
        upper_error =  4
        asymmetric_error = [[lower_error], [upper_error]]
        te = 24
        if plotOn:
            ax.errorbar(99-offset,24,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$24^{+4}_{-5} eV$')
        
    elif pressure == 30:# 26 +-5 eV  
        lower_error =  5
        upper_error =  5
        asymmetric_error = [[lower_error], [upper_error]]
        te = 26
        if plotOn:
            ax.errorbar(99-offset,26,yerr=asymmetric_error,xerr=xerr,capsize=3, label='$26 \pm 5 eV$')
        
    plot_rectangle(99-offset,te, xerr, asymmetric_error, fig, ax,color)

#%%

version = 'v28'
path = '/home/jeff/Research/Export_Control/Cretin/Projects/Neon_Cell/visrad_drive/{}/ne_cell_visrad_xfile.plt'.format(version)
save_dir = '/home/jeff/Research/Export_Control/Cretin/Projects/Neon_Cell/visrad_drive/{}/'.format(version)
helios_dir = '/home/jeff/Research/Code_Projects/Cloudy_Projects/Z_machine/Neon_cell/helios_data'
gen = 'gen1c'
p = '30_torr'
with open(path,'r') as f:
    chunks = f.read().split('#')
    chunks = chunks[1:]
    for chunk in chunks:
        lines = chunk.split('\n')
        info1 = lines[1] # Lines[0] is # line
        info2 = lines[2]
        lines = [line for line in lines if not line.startswith('$')]
        header = lines[0]
        lines = lines[1:]
        # if header.strip().split(':')[-1]== "TEV, TIV, TRADV vs R":
            
        if header.strip() == 'TEV, TIV vs TIME':
            fig, ax = plt.subplots() 
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            tev   = [float(line.split()[1]) for line in lines if not len(line)==0]
            tiv   = [float(line.split()[2]) for line in lines if not len(line)==0]
            #--------------------------------------------------------------------------
            # Helios
            #--------------------------------------------------------------------------
            
            hfile = f"{gen}_{p}_te_avg.txt"
            hpath = os.path.join(helios_dir,hfile)
            hdf = pd.read_csv(hpath)
            
                  
            plt.plot(time,tev,label='CRETIN, $T_e$ ')
            plt.plot(time,tiv,label='CRETIN, $T_i$ ')
            ax.plot(hdf['time'],hdf['te'],label = 'HELIOS $T_e$',color = 'r')
            plt.plot(td_100ns.times,td_100ns.temps,label='Cloudy $T_e$')
            plot_exp('gen1c',30,fig,ax)
            plot_2020_exp('gen1c',30,fig,ax)
            plt.grid(linestyle=':', linewidth=0.5)
            # plt.ylim([0,40])
            # plt.xlim([70,106])
            plt.xlim([0,106])
            plt.xlabel('Time [ns]')
            plt.ylabel('$T_e$ [eV]')
            plt.title('Neon Gas Cell Temperature '+version)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_te_ti'+ext),transparent=False)
            plt.show() 
        elif header.strip() == 'NE vs TIME':
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            ne   = [float(line.split()[1]) for line in lines if not len(line)==0]
            fig, ax = plt.subplots()       
            plt.plot(time,ne,label='CRETIN, $n_e$ ')
            plt.grid(linestyle=':', linewidth=0.5)
            plt.xlabel('Time [ns]')
            plt.ylabel('Electron Density [$cm^{-3}$]')
            plt.title('Neon Gas Cell Electron Density '+version)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_ne'+ext),transparent=False)
            plt.show() 
        elif header.strip() == 'ZBAR vs TIME':
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            zbar  = [float(line.split()[1]) for line in lines if not len(line)==0]
            fig, ax = plt.subplots()       
            plt.plot(time,zbar,label='CRETIN, <Z> ')
            plt.grid(linestyle=':', linewidth=0.5)
            plt.xlabel('Time [ns]')
            plt.ylabel('<Z>')
            plt.title('Neon Gas Cell <Z> '+version)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_zbar'+ext),transparent=False)
            plt.show() 
        elif header.strip() == 'ERAD vs TIME':
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            erad  = [float(line.split()[1]) for line in lines if not len(line)==0]
            fig, ax = plt.subplots()       
            plt.plot(time,erad,label='CRETIN')
            plt.grid(linestyle=':', linewidth=0.5)
            plt.xlabel('Time [ns]')
            plt.ylabel('Radiation Energy Density [ergs/$cm^{3}$]')
            plt.title('Neon Gas Cell Radiation '+version)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_erad'+ext),transparent=False)
            plt.show() 

#%%
name = 'ne_rad'
path = f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_v1/no_dump/{name}.plt'
save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_v1/no_dump/'
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
            hdf = fetch_helios_temp(30)
            fig, ax = plt.subplots()       
            plt.plot(time,tev,label='CRETIN, $T_e$ ')
            plt.plot(time,tiv,label='CRETIN, $T_i$ ')
            ax.plot(hdf['time'],hdf['te'],label = 'HELIOS $T_e$',color = 'r')
            plt.plot(td_100ns.times,td_100ns.temps,label='Cloudy $T_e$')
            plot_exp('gen1c',30,fig,ax)
            plt.grid(linestyle=':', linewidth=0.5)
            # plt.ylim([0,40])
            plt.xlabel('Time [ns]')
            plt.ylabel('$T_e$ [eV]')
            plt.title('Neon Gas Cell Temperature 30 torr')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_te_ti'+ext),transparent=False)
            plt.show() 
        elif header.strip() == 'NE vs TIME':
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            ne   = [float(line.split()[1]) for line in lines if not len(line)==0]
            fig, ax = plt.subplots()       
            plt.plot(time,ne,label='CRETIN, $n_e$ ')
            plt.grid(linestyle=':', linewidth=0.5)
            plt.xlabel('Time [ns]')
            plt.ylabel('Electron Density [$cm^{-3}$]')
            plt.title('Neon Gas Cell Electron Density '+name)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_ne'+ext),transparent=False)
            plt.show() 
        elif header.strip() == 'ZBAR vs TIME':
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            zbar  = [float(line.split()[1]) for line in lines if not len(line)==0]
            fig, ax = plt.subplots()       
            plt.plot(time,zbar,label='CRETIN, <Z> ')
            plt.grid(linestyle=':', linewidth=0.5)
            plt.xlabel('Time [ns]')
            plt.ylabel('<Z>')
            plt.title('Neon Gas Cell <Z> '+name)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_zbar'+ext),transparent=False)
            plt.show() 
        elif header.strip() == 'ERAD vs TIME':
            time  = [float(line.split()[0])*1e9 for line in lines if not len(line)==0]
            erad  = [float(line.split()[1]) for line in lines if not len(line)==0]
            fig, ax = plt.subplots()       
            plt.plot(time,erad,label='CRETIN')
            plt.grid(linestyle=':', linewidth=0.5)
            plt.xlabel('Time [ns]')
            plt.ylabel('Radiation Energy Density [ergs/$cm^{3}$]')
            plt.title('Neon Gas Cell Radiation '+name)
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir,'Cretin_erad'+ext),transparent=False)
            plt.show()         
        

        
       
