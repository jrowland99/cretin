#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 14:10:50 2023

@author: jeff

Working file for CretinData
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
from CretinData import *

ne_ions  = ["Ne","Ne+",  "Ne+2", "Ne+3", "Ne+4", "Ne+5", "Ne+6", "Ne+7", "Ne+8", "Ne+9", "Ne+10"]
plt.rcParams["figure.dpi"] = 250
plt.rcParams["lines.linewidth"] = 2
plt.rcParams.update({"font.size": 15})
plt.rcParams.update({"legend.fontsize": 16})
plt.rcParams["figure.figsize"] = (8, 6)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)
ne_ion_dict = {"Ne+0":'Ne',"Ne+1":'F',  "Ne+2":'O', "Ne+3":'N', "Ne+4":'C', "Ne+5":'B', "Ne+6":'Be', "Ne+7":'Li', "Ne+8":'He', "Ne+9":'H', "Ne+10":'FS'}
def get_visrad():
    base_dir = '/home/jeff/PRISM/VISRAD_files'
    file = 'myl_C_visrad.dat'
    delim = '# --------------------------------------------------------------------'

    # Read the input file
    # with open(os.path.join(base_dir,file), 'r') as f:
    #     lines = f.readlines()
        
    with open(os.path.join(base_dir,file), 'r') as f:
        chunks = f.read().split(delim)
        # Parse the photon energies and fluxes from the file
        chunks = chunks[1:]
        time_list = []
        time_dict = {}
        flag = True
        for chunk in chunks:
            photon_energies = []
            fluxes = []
            if len(chunk)==1:
                continue # ignore empty lines
            else:
                chunk = chunk[1:-1]# Drop the first and last empty lines
            chunk = chunk.split('\n')
            line = chunk[0]
            if "#  INCIDENT SPECTRUM:                   at time =" in line:
                time = float(line.split("=")[-1].strip())
                time_list.append(time)
            chunk = chunk[6:] #Skip the rest of the header
            # for line in chunk: # Now collect the data in each chunk
            #     tokens = line.split()
            #     photon_energies.append(float(tokens[0]))
            #     fluxes.append(float(tokens[1]))
            for i,line in enumerate(chunk): # Now collect the data in each chunk
                if i%2 == 0 : #skip every other row
                    continue
                tokens = line.split()
                photon_energies.append(float(tokens[0]))
                fluxes.append(float(tokens[1]))
            time_dict[time]=fluxes
            if flag:
                flag=False
                photon_energies_b = photon_energies
    integ_fluxes = []
    xaxis = np.asarray(photon_energies_b)
    for i, time in enumerate(time_list):
        flux = np.asarray(time_dict[time])
        integ_fluxes.append(np.trapz(y=flux,x=xaxis))
        
    return {'fluxes':fluxes,'energy':photon_energies_b,'integ fluxes':integ_fluxes,'time':np.asarray(time_list)}

#%%Eigenvalues average 
# gen = 'gen1c'
# pressure = 7.5
# params = {
#     'name':'ne_cell',
#     'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v13/{gen}/{pressure}torr/',
#     'pressure': pressure,
#     'gen': gen}
# testobj = CretinData(**params)
# testobj.read_data()
# zones = 51

# # print(testobj.time_data.keys())
# sum_arr = np.zeros(len(testobj.time_data['time']))
# for i in np.arange(1,zones+1):
#     eig_str = f"eigenv_r 1,{i},0,-1"
#     for j,val in enumerate(testobj.time_data[eig_str]):
#         sum_arr[j] = sum_arr[j] + val
        
# eigenvalues_space_avg = -sum_arr/zones


save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v13'
for pressure in [7.5,15,30]:
    params = {
        'name':'ne_cell',
        'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v13/{gen}/{pressure}torr/',
        'pressure': pressure,
        'gen': gen}
    testobj = CretinData(**params)
    testobj.read_data()
    
    sum_arr = np.zeros(len(testobj.time_data['time']))
    for i in np.arange(1,zones+1):
        eig_str = f"eigenv_r 1,{i},0,-1"
        for j,val in enumerate(testobj.time_data[eig_str]):
            sum_arr[j] = sum_arr[j] + val
            
    eigenvalues_space_avg = -sum_arr/zones
    slow_time = 1/eigenvalues_space_avg
    # fast_time = 1/np.array(np.abs(testobj.time_data[max_eig_str]))
    plt.plot(testobj.time_data['time']*1e9,slow_time*1e9,label=f"{pressure} torr")
    # plt.plot(testobj.time_data['time'],fast_time,label='fast')
    # plt.yscale('log')
    plt.ylabel('AK timescale [ns]')
    plt.xlabel('exp time [ns]')
    # plt.xlim([60,108])
    plt.grid(linestyle=':', linewidth=0.5)
    plt.legend()
    
plt.title("Cretin Slow atomic kinetics timescale\n 51 zone space avg, gen1c")
plt.tight_layout()
plt.savefig(os.path.join(save_dir,f'cretin_AK_slow_timescale'+ext), transparent=False)

#%%
keyword_list = []
keyword_list.append(
{
    'edit'        : 'NET HEAT vs R',
    'xvar'        : 'r',
    'yvar'        : 'heatt'
    }
    )
gen = 'gen1c'
pressure = 30
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/steady_state_latest/sw36_off/gen1c/30torr/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()
for keyword in keyword_list:
    testobj.make_animation(**keyword)


#%% Analysis for line transfer 
# keyword_list = []
# keyword_list.append(
# {
#     'edit'        : 'NET HEAT vs R',
#     'xvar'        : 'r',
#     'yvar'        : 'heatt'
#     }
#     )
gen = 'gen1c'
pressure = 7.5
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/lines/v9/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()

# print(testobj.data.keys())
# print(testobj.time_data.keys())
# print(testobj.data['TAU vs ENERGY'])
plt.plot(testobj.time_data['time']*1e9,testobj.time_data['jbarcgs 1,1'],label = r'$H-\alpha$')
plt.plot(testobj.time_data['time']*1e9,testobj.time_data['jbarcgs 2,1'],label = r'$He-\alpha$')
plt.plot(testobj.time_data['time']*1e9,testobj.time_data['jbarcgs 3,1'],label = r'$Li-\alpha$')
# plt.xlim([80,110])
plt.legend(loc='lower left')
plt.grid(linestyle=':', linewidth=0.5)  
plt.xlabel('Experiment time [ns]') 
plt.ylabel('Line strength $[erg/cm^2/sec/Hz/ster]$') 
plt.title(f'Cretin line strength in zone 1 \n {gen} {pressure} torr')
plt.tight_layout()

plt.show()
keyword_list = []
keyword_list.append(
{
    'edit'        : 'SPECTRUM vs ENERGY',
    'xvar'        : 'sp_energy',
    'yvar'        : 'jsp 0,51'
    }
    )
keyword_list.append(
{
    'edit'        : 'TAU vs ENERGY',
    'xvar'        : 'energy',
    'yvar'        : 'taukap 0,51'
    }
    )
for keyword in keyword_list:
    testobj.make_animation(**keyword)
#%%
# keyword_list = []
# keyword_list.append(
# {
#     'edit'        : 'NET HEAT vs R',
#     'xvar'        : 'r',
#     'yvar'        : 'heatt'
#     }
#     )
gen = 'gen1c'
pressure = 30
params = {
    'name':'ne_cell',
    # 'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v7/',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/spectrum/v1/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()
# testobj.data.keys()
keyword_list = []
keyword_list.append(
{
    'edit'        : 'NET HEAT vs R',
    'xvar'        : 'r',
    'yvar'        : 'heatt'
    }
    )
keyword_list.append(
{
    'edit'        : 'ERAD vs R',
    'xvar'        : 'r',
    'yvar'        : 'erad'
    }
    )
keyword_list.append(
{
    'edit'        : 'JNU vs E',
    'xvar'        : 'energy',
    'yvar'        : 'jnu 0,1'
    }
    )
for keyword in keyword_list:
    testobj.make_animation(**keyword)

# keyword = {
#     'edit1'         : 'TEV, TIV, TRADV vs R',
#     'edit2'         : 'NE vs R',
#     'xvar'          : 'r',
#     'yvar1'         : 'tev',
#     'yvar2'         : 'ne'
#     }
keyword = {
    'edit1'         : 'TEV, TIV, TRADV vs R',
    'edit2'         : 'NET HEAT vs R',
    'xvar'          : 'r',
    'yvar1'         : 'tev',
    'yvar2'         : 'heatt'
    }
testobj.coplot_movie(**keyword)
#%% destruction timescale
import matplotlib.cm as cm
cmap = cm.get_cmap('plasma') 
gen1c_drive = get_visrad()
time_data_dict = {}
fig, ax = plt.subplots() 
for i in range(0,11):
    fstr = f"gammatot<r> 1,0,{i},1"
    dest_rate = np.array(testobj.time_data[fstr])
    dest_time = np.reciprocal(dest_rate)
    dest_time[np.isinf(dest_time)] = np.nan
    time_data_dict[fstr]= dest_time
    ion_charge = 10 - i
    if i%2==0:
        ax.plot(testobj.time_data['time']*1e9-100,dest_time,label = ne_ion_dict[f'Ne+{ion_charge}']+'-Like',color = cmap(i/11))
    else:
        ax.plot(testobj.time_data['time']*1e9-100,dest_time,label =ne_ion_dict[ f'Ne+{ion_charge}']+'-Like',linestyle='dashed',color = cmap(i/11))
ax.axhline(y=3e-9, color='b', linestyle='dashdot',label='Z FWHM')
ax2 = ax.twinx()
ax2.plot(gen1c_drive['time']*1e9-100,gen1c_drive['integ fluxes'],linestyle = ':',color='black')
ax2.set_ylabel('Flux $[J/cm^2/s]$')

ax.set_yscale('log')
# ax.set_ylim([1e-12,1e-3])
ax.legend(loc='lower left')
ax.grid(linestyle=':', linewidth=0.5)  
ax.set_xlabel('Relative Time [ns]') 
ax.set_ylabel('Destruction Timescale [s]') 
plt.title(r'Cretin Destruction timescale Ne cell 30 torr, SiN windows, sw(36)=1')
plt.tight_layout()
# labelLines(plt.gca().get_lines(), align=False,backgroundcolor='none', fontsize=10)
plt.savefig(os.path.join(params['path'],'cretin_destr_timescale'+ext), transparent=False)    
#%% USES A LOT OF MEMORY
gen = 'gen1c'
pressure = 30
params1 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
data1 = CretinData(**params1)
data1.read_data()
time_ax1 = data1.time_data['time'] * 1e9

params2 = {
    'name': 'ne_cell',
   
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/steady_state_latest/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
data2 = CretinData(**params2)
data2.read_data()
labels = {
    'label1' : 'TD',
    'label2' : 'SS',
    'titlestr': '\n SS vs TD  gen1c 30 Torr'
    }
two_csd_animation(data1,data2,f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/steady_state_latest/{gen}/{pressure}torr/two_csd_test.mp4',**labels)

#%% Get two separate CRETIN results, interpolate the first to match the time axis of the second for co-plotting


gen = 'gen1c'
pressure = 30
params1 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
data1 = CretinData(**params1)
data1.read_data()
time_ax1 = data1.time_data['time'] * 1e9

params2 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/lines/v6/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
data2 = CretinData(**params2)
data2.read_data()
new_x = data2.time_data['time'] * 1e9

csd_dict1 = {}
csd_interp_dict = {}
for itime, time in enumerate(time_ax1):
    csd_list1 = []
    for ion in data1.ne_ions:
        ionfrac = data1.time_data[data1.yiso_labels[ion]][itime]
        csd_list1.append(ionfrac)
    csd_dict1[time] = csd_list1
for ion in data1.ne_ions:
    ionfrac     = data1.time_data[data1.yiso_labels[ion]]
    interp_ion  = np.interp(new_x, time_ax1, ionfrac)
    csd_interp_dict[ion] = interp_ion
    
csd_dict2 = {}
csd_interp_dict2 = {}
for itime, time in enumerate(new_x):
    csd_list2 = []
    for ion in data2.ne_ions:
        ionfrac = data2.time_data[data2.yiso_labels[ion]][itime]
        csd_list2.append(ionfrac)
    csd_dict2[time] = csd_list2
for ion in data2.ne_ions:
    ionfrac     = data2.time_data[data2.yiso_labels[ion]]
    interp_ion  = ionfrac
    csd_interp_dict2[ion] = interp_ion



for index, time in enumerate(new_x):
    if round(time) == 90:
        start_t_index = index


save_dir = os.path.join(params2['path'],'CSD_plots')
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

for itime, time in enumerate(new_x):
    if itime< start_t_index:
        continue
    csd1,csd2 = [],[]
    
    for ion in ne_ions:
        csd1.append(csd_interp_dict[ion][itime])
        csd2.append(csd_interp_dict2[ion][itime])
    time  = round(time)
    fig, ax = plt.subplots()
    plt.plot(np.arange(len(ne_ions)),csd1,label = 'TD')
    plt.plot(np.arange(len(ne_ions)),csd2,label = 'TD w/ lines')
    plt.title(f"Cretin Line transfer comparison CSD,\n P = {pressure} torr, t = {time} ns")
    plt.xticks(np.arange(len(ne_ions)),ne_ion_dict.values())
    # plt.xlim([2,10])
    plt.ylim([0,0.9])
    plt.grid(linestyle=':', linewidth=0.5)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir,f'CSD_{pressure}torr_cretin_lines_{time}ns'+ext), transparent=False)
    plt.show()
    plt.close()

# interpolated_dict = {}

# keys = list(csd_dict1.keys())
# values = list(csd_dict1.values())

# num_ions = len(values[0])  # Assuming all lists have the same length
# for itime,time in enumerate(keys):
#     csd_list = values[itime]
#     for ion in range(num_ions):
#         csd_list = [values[j][i] for j in range(len(keys))]
#         interpolated_y = np.interp(new_x, keys, y_values)
#         interpolated_dict[new_x[i]] = interpolated_y
#%%
gen = 'gen1c'
pressure = 30
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/gen1c/30torr/more_edits/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()
keyword_list = []
# keyword_list.append(
# {
#     'edit'        : 'TEV, TIV, TRADV vs R',
#     'xvar'        : 'r',
#     'yvar'        : 'tev'
#     }
#     )
keyword_list.append(
{
    'edit'        : 'RAD ABS HEAT vs R',
    'xvar'        : 'r',
    'yvar'        : 'heatj 1'
    }
    )
keyword_list.append(
{
    'edit'        : 'RAD EMIS COOL vs R',
    'xvar'        : 'r',
    'yvar'        : 'coolj 1'
    }
    )
keyword_list.append(
{
    'edit'        : 'ECOL COOL vs R',
    'xvar'        : 'r',
    'yvar'        : 'coolc 1'
    }
    )
keyword_list.append(
{
    'edit'        : 'ICOL COOL vs R',
    'xvar'        : 'r',
    'yvar'        : 'cooli 1'
    }
    )
# print(testobj.data['RAD ABS HEAT vs R'])
for keyword in keyword_list:
    testobj.make_animation(**keyword)
#%% keyword list
keyword_list = []
keyword_list.append(
{
    'edit'        : 'TEV, TIV, TRADV vs R',
    'xvar'        : 'r',
    'yvar'        : 'tev'
    }
    )
# keyword_list.append(
# {
#     'edit'        : 'TEV, TIV, TRADV vs R',
#     'xvar'        : 'r',
#     'yvar'        : 'tiv'
#     }
#     )
# keyword_list.append(
# {
#     'edit'        : 'JNU vs E',
#     'xvar'        : 'energy',
#     'yvar'        : 'jnu 0,1'
#     }
#     )
# keyword_list.append(
# {
#     'edit'        : 'JNU vs E',
#     'xvar'        : 'energy',
#     'yvar'        : 'jnu 0,26'
#     }
#     )
# keyword_list.append(
# {
#     'edit'        : 'JNU vs E',
#     'xvar'        : 'energy',
#     'yvar'        : 'jnu 0,51'
#     }
#     )
keyword_list.append(
{
    'edit'        : 'ERAD vs R',
    'xvar'        : 'r',
    'yvar'        : 'erad'
    }
    )
# keyword_list.append(
# {
#     'edit'        : 'TAU vs E',
#     'xvar'        : 'energy',
#     'yvar'        : 'taukap'
#     }
#     )
# keyword_list.append(
# {
#     'edit'        : 'ZBAR vs R',
#     'xvar'        : 'r',
#     'yvar'        : 'zbar'
#     }
#     )

keyword_list.append(
{
    'edit'        : 'NE vs R',
    'xvar'        : 'r',
    'yvar'        : 'ne'
    }
    )
# keyword_list.append(
# {
#     'edit'        : 'TAU vs R',
#     'xvar'        : 'r',
#     'yvar'        : 'taukap'
#     }
#     )
# keyword_list.append(

#     )
#%%
for gen in ['gen1c']:
# for gen in ['gen1c']:
    for pressure in [7.5,15,30]:
        params = {
            'name':'ne_cell',
            'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_hydro_temp_eq/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen}
        testobj = CretinData(**params)
        testobj.read_data()
        for keywords in keyword_list:
            testobj.make_animation(**keywords)
        
        keyword = {
            'edit1'         : 'TEV, TIV, TRADV vs R',
            'edit2'         : 'NE vs R',
            'xvar'          : 'r',
            'yvar1'         : 'tev',
            'yvar2'         : 'ne'
            }
        # keyword = {
        #     'edit1'         : 'TEV, TIV, TRADV vs R',
        #     'edit2'         : 'ERAD vs R',
        #     'xvar'          : 'r',
        #     'yvar1'         : 'tev',
        #     'yvar2'         : 'erad'
        #     }
        testobj.coplot_movie(**keyword)
#%%
gen = 'gen1c'
pressure = 30
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/lines/v12/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()

keyword = {
    'edit'         : 'TEV, TIV, TRADV vs R',
    'xvar'          : 'r',
    'yvar'         : 'tev',
    }
testobj.make_derivative(**keyword)
print(testobj.data['DERIVATIVES'])

keyword = {
    'edit'         : 'DERIVATIVES',
    'xvar'          : 'r',
    'yvar'         : 'dtev_dr',
    }
testobj.make_animation(**keyword)
#%% Now do a keyword list for more derivatives
gen = 'gen1c'
pressure = 30
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/new_bc_v3/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()

keyword_list = []
keyword_list.append(
{
    'edit'         : 'TEV, TIV, TRADV vs R',
    'xvar'          : 'r',
    'yvar'         : 'tev',
    }
    )
keyword_list.append(
{
    'edit'        : 'NE vs R',
    'xvar'        : 'r',
    'yvar'        : 'ne'
    }
    )
for keywords in keyword_list:
    testobj.make_derivative(**keywords)

keyword_list2 = []
keyword_list2.append(
    {
        'edit'         : 'dtev_dr',
        'xvar'          : 'r',
        'yvar'         : 'dtev_dr',
        }
    )
keyword_list2.append(
    {
        'edit'         : 'dne_dr',
        'xvar'          : 'r',
        'yvar'         : 'dne_dr',
        }
    )
for keyword in keyword_list2:
    testobj.make_animation(**keyword)
#%%

keyword = {
    'edit1'         : 'dtev_dr',
    'edit2'         : 'dne_dr',
    'xvar'          : 'r',
    'yvar1'         : 'dtev_dr',
    'yvar2'         : 'dne_dr'
    }
# keyword = {
#     'edit1'         : 'TEV, TIV, TRADV vs R',
#     'edit2'         : 'ERAD vs R',
#     'xvar'          : 'r',
#     'yvar1'         : 'tev',
#     'yvar2'         : 'erad'
#     }
testobj.coplot_movie(**keyword)
#%%
for gen in ['gen1c','gen2c','gen2f']:
# for gen in ['gen1c']:
    for pressure in [30]:
        params = {
            'name':'ne_cell',
            'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/sw36_off/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen}
        testobj = CretinData(**params)
        testobj.read_data()
        for keywords in keyword_list:
            testobj.make_animation(**keywords)
        
        keyword = {
            'edit1'         : 'TEV, TIV, TRADV vs R',
            'edit2'         : 'NE vs R',
            'xvar'          : 'r',
            'yvar1'         : 'tev',
            'yvar2'         : 'ne'
            }
        # keyword = {
        #     'edit1'         : 'TEV, TIV, TRADV vs R',
        #     'edit2'         : 'ERAD vs R',
        #     'xvar'          : 'r',
        #     'yvar1'         : 'tev',
        #     'yvar2'         : 'erad'
        #     }
        testobj.coplot_movie(**keyword)
#%% make animations for keyword list

for gen in ['gen2c','gen2f']:
    for pressure in [7.5,15,30]:
        params = {
            'name':'ne_cell',
            'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/new_bc_v2/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen}
        testobj = CretinData(**params)
        testobj.read_data()
        for keywords in keyword_list:
            testobj.make_animation(**keywords)
            
#%% Do single run case for one cretin model
keyword_list = []
keyword_list.append(
{
    'edit'        : 'YISO vs R 1',
    'xvar'        : 'r',
    'yvar'        : 'yisofrac 1'
    }
    )
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/new_bc_v2/gen1c/30torr/',
    'pressure': 30,
    'gen': 'gen1c'}
testobj = CretinData(**params)
testobj.read_data()
for keywords in keyword_list:
    testobj.make_animation(**keywords)
#%% do t-dep fractional population plots
for gen in ['gen1c','gen2c','gen2f']:
    for pressure in [7.5,15,30]:
        params = {
            'name':'ne_cell',
            'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen}
        testobj = CretinData(**params)
        testobj.read_data()
        fig, ax = plt.subplots()
        for ion in testobj.ne_ions:
            save_dir = f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/{gen}/{pressure}torr/'
            plt.plot(testobj.time_data['time']*1e9,testobj.time_data[testobj.yiso_labels[ion]],label = ion)
        # plt.xlim([60,108])
        # plt.yscale('log')
        plt.xlabel('Time [ns]')
        plt.ylabel('Fraction')
        plt.title(f'Cretin Ne Cell {testobj.gen} {testobj.pressure} torr')
        # plt.legend()
        labelLines(plt.gca().get_lines(), align=False,backgroundcolor='none', fontsize=10)
        plt.grid(linestyle=':', linewidth=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(save_dir,f'{testobj.gen}_{testobj.pressure}_cretin_yiso'+ext),transparent=False)
        plt.show()

#%% T-dep CSD plots with zbar
def zbar(csd_list):
    zbar = 0
    for charge, pop in enumerate(csd_list):
        # the index is the free charge, and the frac pop is the value
        zbar = zbar + charge * pop
    return zbar

for gen in ['gen1c', 'gen2c', 'gen2f']:
    for pressure in [7.5, 15, 30]:
        print(f"starting {gen} {pressure}")
        params = {
            'name': 'ne_cell',
            'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen
        }
        testobj = CretinData(**params)
        testobj.read_data()

        # Prepare the list of figure objects
        fig_list = []
        time_ax = testobj.time_data['time'] * 1e9
        for itime, time in enumerate(time_ax):
            csd_list = []
            fig, ax = plt.subplots()
            for ion in testobj.ne_ions:
                ionfrac = testobj.time_data[testobj.yiso_labels[ion]][itime]
                csd_list.append(ionfrac)
            plt.plot(range(11), csd_list)
            zavg = zbar(csd_list)
            plt.plot([zavg, zavg], [0, 0.1], color='red', linestyle='--')

            plt.title(f'Cretin Ne Cell {testobj.gen} {testobj.pressure} torr CSD at t = {time:.2f} ns, ')
            plt.ylim([0, 1])
            plt.grid(linestyle=':', linewidth=0.5)
            plt.tight_layout()
            fig_list.append(fig)
            plt.close()

        # Create the animation
        anifig, ax = plt.subplots()  # Create a figure for the animation
        line, = ax.plot(range(11), np.zeros(11))  # Initialize an empty line plot
        lines = [line]

        def animate(itime):
            # Update the existing line plot with the new data
            line.set_ydata([testobj.time_data[testobj.yiso_labels[ion]][itime] for ion in testobj.ne_ions])
            zavg = zbar([testobj.time_data[testobj.yiso_labels[ion]][itime] for ion in testobj.ne_ions])
            ax.set_title(f'Cretin Ne Cell {testobj.gen} {testobj.pressure} torr CSD at t = {time_ax[itime]:.2f} ns, <Z> = {zavg:.2f}')
            ax.set_ylim([-0.05,1])
            ax.set_xticks(range(11))
            ax.set_xlabel('$Ne^{+q}$')
            plt.grid(linestyle=':', linewidth=0.5)
        ani = FuncAnimation(anifig, animate, frames=len(time_ax), blit=False)

        # Save the animation
        save_dir = testobj.path
        ani.save(save_dir + testobj.title_str + '_CSD.mp4', writer='ffmpeg')
#%%
gen = 'gen1c'
pressure = 30
params1 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
obj1 = CretinData(**params1)
obj1.read_data()

params2 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2_SS/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
# params2 = params1
obj2 = CretinData(**params2)
obj2.read_data()

labels = {
    'label1' : 'TD No hydro',
    'label2' : 'SS',
    'titlestr': 'SS vs TD'
    }

two_csd_animation(obj1,obj2,'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/two_csd_test.mp4',**labels)

#%% Get comparison 
for gen in ['gen1c', 'gen2c', 'gen2f']:
    for pressure in [7.5, 15, 30]:
        print(f"starting {gen} {pressure}")
        params1 = {
            'name': 'ne_cell',
            'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen
        }
        obj1 = CretinData(**params1)
        obj1.read_data()
        
        params2 = {
            'name': 'ne_cell',
            'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_w_hydro/{gen}/{pressure}torr/',
            'pressure': pressure,
            'gen': gen
        }
        # params2 = params1
        obj2 = CretinData(**params2)
        obj2.read_data()
        
        labels = {
            'label1' : 'TD No hydro',
            'label2' : 'TD W Hydro',
            'titlestr': 'Hydro On/off'
            }
        
        two_csd_animation(obj1,obj2,f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/{gen}_{pressure}_CSD_hydro_comp.mp4',**labels)
        

# print(testobj.data.keys())
# print(testobj.data['NE vs R'].keys())
# print(testobj.time_data.keys())
# print(testobj.time_data['NE'])
#%%

# import matplotlib.patches as patches
plt.rcParams["figure.dpi"] = 250
plt.rcParams["lines.linewidth"] = 1
plt.rcParams.update({"font.size": 10})
plt.rcParams.update({"legend.fontsize": 10})
plt.rcParams["figure.figsize"] = (4, 4)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)

time_ax     = np.asarray(list(testobj.data['NE vs R'].keys()))
space_ax    = list(testobj.data['NE vs R'][time_ax[0]]['r'])
space_ax    = np.asarray([float(r) for r in space_ax])
eden        = np.array([testobj.data['NE vs R'][time]['ne'].values for time in time_ax])
eden        = eden.astype(float)

plt.figure(figsize=(10, 6))
plt.imshow(eden, aspect='auto', extent=[0, len(time_ax), space_ax[-1], space_ax[0]],
            interpolation='none', cmap='viridis', norm=LogNorm())

plt.colorbar(label='ne')
plt.xlabel('Time')
plt.ylabel('r')
plt.title('Heatmap Plot')
# plt.xticks(np.arange(len(time_ax)), time_ax, rotation=45)
plt.gca().invert_yaxis()

plt.show()
#%%
edit        = 'ZBAR vs R'
time_ax     = np.asarray(list(testobj.data[edit].keys()))
space_ax    = list(testobj.data[edit][time_ax[0]]['r'])
space_ax    = np.asarray([float(r) for r in space_ax])
z           = np.array([testobj.data[edit][time]['zbar'].values for time in time_ax])
z           = z.astype(float)

plt.figure(figsize=(10, 6))
plt.imshow(z, aspect='auto', extent=[0,time_ax[-1], space_ax[-1], space_ax[0]],
            interpolation='none', cmap='viridis', norm=LogNorm())

plt.colorbar(label='zbar')
plt.xlabel('Time')
plt.ylabel('r')
plt.title(edit)
# plt.xticks(np.arange(len(time_ax)), time_ax, rotation=45)
plt.gca().invert_yaxis()

plt.show()
#%%
edit        = 'TEV, TIV, TRADV vs R'
time_ax     = np.asarray(list(testobj.data[edit].keys()))
space_ax    = list(testobj.data[edit][time_ax[0]]['r'])
space_ax    = np.asarray([float(r) for r in space_ax])
z           = np.array([testobj.data[edit][time]['tev'].values for time in time_ax])
z           = z.astype(float)
print(z)
time_ax2 = 1e9*time_ax

plt.figure(figsize=(10, 6))
plt.imshow(z, aspect='auto', extent=[0,time_ax2[-1], space_ax[-1], space_ax[0]],
            interpolation='none', cmap='viridis', norm=Normalize(vmin=min(z[0]), vmax=max(z[-1])))
# plt.imshow(z,cmap='viridis', aspect='auto')
plt.colorbar(label='tev [eV]')
plt.xlabel('Time [ns]')
plt.ylabel('R')
plt.title(edit)
# plt.xticks(np.arange(len(time_ax)), time_ax, rotation=45)
plt.gca().invert_yaxis()

plt.show()
#%%
from matplotlib.animation import FuncAnimation
edit        = 'TEV, TIV, TRADV vs R'
time_ax     = np.asarray(list(testobj.data[edit].keys()))
space_ax    = list(testobj.data[edit][time_ax[0]]['r'])
space_ax    = np.asarray([float(r) for r in space_ax])
# z           = np.array([testobj.data[edit][time]['tev'].values for time in time_ax])
# z           = z.astype(float)

# time_ax2 = 1e9*time_ax
# for time in time_ax:
#     plt.figure(figsize=(3, 4))
#     z           = np.array(testobj.data[edit][time]['tev'].values)
#     z           = z.astype(float)
#     plt.plot(space_ax,z,label = 'Te')
#     plt.ylim([0,27])
#     plt.xlabel('R [cm]')
#     plt.ylabel('Te [eV]')
#     plt.title(str(time*1e9) + 'ns')
#     plt.show()

# Prepare the figure and axis
fig, ax = plt.subplots(figsize=(3, 4))
line, = ax.plot([], [], label='Te')
ax.set_ylim(0, 27)
ax.set_xlim([0,1.4])
ax.set_xlabel('R [cm]')
ax.set_ylabel('Te [eV]')

# Initialize the plot elements
def init():
    line.set_data([], [])
    ax.set_title('')
    return line,

# Animation update function
def animate(time):
    z = np.array(testobj.data[edit][time]['tev'].values).astype(float)
    line.set_data(space_ax, z)
    ax.set_title('Cretin $T_e$ at t = '+str(round(time*1e9)) + 'ns')
    return line,

# Create the animation
ani = FuncAnimation(fig, animate, init_func=init, frames=time_ax, blit=True)

save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/'
# Show the animation
ani.save(save_dir+'animation.mp4', writer='ffmpeg')
plt.show()


# # plt.imshow(z,cmap='viridis', aspect='auto')
# plt.colorbar(label='tev [eV]')
# plt.xlabel('Time [ns]')
# plt.ylabel('R')
# plt.title(edit)
# # plt.xticks(np.arange(len(time_ax)), time_ax, rotation=45)
# plt.gca().invert_yaxis()

plt.show()
#%%
keywords = {
    'save_dir'    : '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/',
    'edit'        : 'TEV, TIV, TRADV vs R',
    'xvar'        : 'r',
    'yvar'        : 'tev'
    }
testobj.make_animation(**keywords)
#%%
print(testobj.data.keys())
keywords = {
    'save_dir'    : '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/',
    'edit'        : 'JNU vs E',
    'xvar'        : 'energy',
    'yvar'        : 'jnu 0,1'
    }
testobj.make_animation(**keywords)
#%%

plt.rcParams["figure.dpi"] = 250
plt.rcParams["lines.linewidth"] = 1
plt.rcParams.update({"font.size": 10})
plt.rcParams.update({"legend.fontsize": 10})
plt.rcParams["figure.figsize"] = (7, 4)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)

fig, ax = plt.subplots()
for ion in testobj.ne_ions:
    save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_trans_grid_v2/'
    plt.plot(testobj.time_data['time']*1e9,testobj.time_data[testobj.yiso_labels[ion]],label = ion)
plt.xlim([60,108])
# plt.yscale('log')
plt.xlabel('Time [ns]')
plt.ylabel('Fraction')
plt.title(f'Cretin Ne Cell {testobj.gen} {testobj.pressure} torr')
# plt.legend()
labelLines(plt.gca().get_lines(), align=False,backgroundcolor='none', fontsize=10)
plt.grid(linestyle=':', linewidth=0.5)

plt.tight_layout()
plt.savefig(os.path.join(save_dir,'Cretin_yiso_log'+ext),transparent=False)
plt.show()

#%%
from scipy.interpolate import interp1d
gen = 'gen1f'
pressure = 60
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/gen1f_georges/gen1f/60torr/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()
eden_df_dict = testobj.data['NE vs R']
times = np.array(list(eden_df_dict.keys()))
pdv_positions = [.2348,.5328,.8308,1.1287]
pdv_data_dict = {1:[],2:[],3:[],4:[]}
for df in eden_df_dict.values():
    r = df['r'].to_numpy(dtype=np.float64)
    ne = df['ne'].to_numpy(dtype=np.float64)
    linear_interp = interp1d(r, ne, kind='linear')
    for index, value in enumerate(pdv_positions):
        pdv_data_dict[index+1].append(linear_interp(value))

for i in range(len(pdv_data_dict)):
    plt.plot(times*1e9-100,pdv_data_dict[i+1],label = f'Probe {i+1}')
       
plt.title(f"Cretin $n_e$ vs time 15 torr, mylar windows, far pos.")
plt.xlim([-40,10])
plt.legend()
plt.ylabel('$n_e$ [$cm^{-3}$]')
plt.xlabel('Experiment time [ns]')
plt.grid(linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.savefig(os.path.join(params['path'],'cretin_4probe_ne_vs_time_peak'+ext), transparent=False)   

save_df = pd.DataFrame(pdv_data_dict)
save_df['time'] = times
save_df.to_csv(os.path.join(params['path'],f'cretin_4probe_PDV_{pressure}torr.csv'),index=False)
#%% CSD Plots for 3 objects
# Define parameters for three objects
gen = 'gen1c'
pressure = 7.5
params1 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
params2 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/const_td_temp/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
params3 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/steady_state_latest/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}

params1 = {
    'name': 'ne_cell',
    'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen
}
data1 = CretinData(**params1)
data1.read_data()
time_ax1 = data1.time_data['time'] * 1e9

data2 = CretinData(**params2)
data2.read_data()
new_x = data2.time_data['time'] * 1e9

data3 = CretinData(**params3)
data3.read_data()
time_ax3 = data3.time_data['time'] * 1e9

csd_dict1 = {}
csd_interp_dict = {}
for itime, time in enumerate(time_ax1):
    csd_list1 = []
    for ion in data1.ne_ions:
        ionfrac = data1.time_data[data1.yiso_labels[ion]][itime]
        csd_list1.append(ionfrac)
    csd_dict1[time] = csd_list1
for ion in data1.ne_ions:
    ionfrac     = data1.time_data[data1.yiso_labels[ion]]
    interp_ion  = np.interp(new_x, time_ax1, ionfrac)
    csd_interp_dict[ion] = interp_ion
    
csd_dict2 = {}
csd_interp_dict2 = {}
for itime, time in enumerate(new_x):
    csd_list2 = []
    for ion in data2.ne_ions:
        ionfrac = data2.time_data[data2.yiso_labels[ion]][itime]
        csd_list2.append(ionfrac)
    csd_dict2[time] = csd_list2
for ion in data2.ne_ions:
    ionfrac     = data2.time_data[data2.yiso_labels[ion]]
    interp_ion  = ionfrac
    csd_interp_dict2[ion] = interp_ion
    
csd_dict3 = {}
csd_interp_dict3 = {}
for itime, time in enumerate(time_ax3):
    csd_list3 = []
    for ion in data3.ne_ions:
        ionfrac = data3.time_data[data3.yiso_labels[ion]][itime]
        csd_list3.append(ionfrac)
    csd_dict3[time] = csd_list3
for ion in data3.ne_ions:
    ionfrac     = data3.time_data[data3.yiso_labels[ion]]
    interp_ion  = np.interp(new_x, time_ax3, ionfrac)
    csd_interp_dict3[ion] = interp_ion

td_temp_interp = np.interp(new_x, time_ax1, data1.time_data['TEV'])
ss_temp_interp = np.interp(new_x, time_ax3, data3.time_data['TEV'])

start_t_index = 94
save_dir = os.path.join(params2['path'],'CSD_plots')
if not os.path.exists(save_dir):
    os.mkdir(save_dir)
for itime, time in enumerate(new_x):
    if itime< start_t_index:
        continue
    csd1,csd2,csd3 = [],[],[]
    
    for ion in ne_ions:
        csd1.append(csd_interp_dict[ion][itime])
        csd2.append(csd_interp_dict2[ion][itime])
        csd3.append(csd_interp_dict3[ion][itime])
    td_temp = round(td_temp_interp[itime],1)
    ss_temp = round(ss_temp_interp[itime],1)
    time  = round(time)
    fig, ax = plt.subplots()
    plt.plot(np.arange(len(ne_ions)),csd1,label = f'TD $T_e$ = {td_temp} eV')
    plt.plot(np.arange(len(ne_ions)),csd2,label = f'SS w/ TD $T_e$ = {td_temp} eV')
    plt.plot(np.arange(len(ne_ions)),csd3,label = f'SS self det $T_e$ = {ss_temp} eV')
    plt.title(f"Cretin SS vs TD CSD w/ Rad. Trans,\n P = {pressure} torr, t = {time} ns")
    plt.xticks(np.arange(len(ne_ions)),ne_ion_dict.values())
    # plt.xlim([2,10])
    plt.ylim([0,0.9])
    plt.grid(linestyle=':', linewidth=0.5)
    plt.legend(loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir,f'CSD_{pressure}torr_cretin_SS_TD_{time}ns'+ext), transparent=False)
#%%
min_eig_str = 'eigenv_r 1,1,0,-1'
max_eig_str = 'eigenv_r 1,1,0,1'
gen = 'gen1c'
save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v12'
for pressure in [7.5,15,30]:
    params = {
        'name':'ne_cell',
        'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v12/{gen}/{pressure}torr/',
        'pressure': pressure,
        'gen': gen}
    testobj = CretinData(**params)
    testobj.read_data()
    slow_time = 1/np.array(np.abs(testobj.time_data[min_eig_str]))
    fast_time = 1/np.array(np.abs(testobj.time_data[max_eig_str]))
    plt.plot(testobj.time_data['time']*1e9,slow_time*1e9,label=f"{pressure} torr")
    # plt.plot(testobj.time_data['time'],fast_time,label='fast')
    # plt.yscale('log')
    plt.ylabel('AK timescale [ns]')
    plt.xlabel('exp time [ns]')
    plt.grid(linestyle=':', linewidth=0.5)
    plt.legend()
    plt.tight_layout()
plt.title("Cretin 1 zone Slow atomic kinetics timescale, gen1c")
plt.savefig(os.path.join(save_dir,f'cretin_AK_slow_timescale'+ext), transparent=False)
#%%
min_eig_str = 'eigenv_r 1,1,0,-1'
max_eig_str = 'eigenv_r 1,1,0,1'
gen = 'gen1c'
save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v12'
for pressure in [7.5,15,30]:
    params = {
        'name':'ne_cell',
        'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v12/{gen}/{pressure}torr/',
        'pressure': pressure,
        'gen': gen}
    testobj = CretinData(**params)
    testobj.read_data()
    slow_time = 1/np.array(np.abs(testobj.time_data[min_eig_str]))
    fast_time = 1/np.array(np.abs(testobj.time_data[max_eig_str]))
    # plt.plot(testobj.time_data['time']*1e9,slow_time*1e9,label=f"{pressure} torr")
    plt.plot(testobj.time_data['time']*1e9,fast_time,label=f"{pressure} torr")
    # plt.yscale('log')
    plt.ylabel('AK timescale [s]')
    plt.xlabel('exp time [ns]')
    plt.grid(linestyle=':', linewidth=0.5)
    plt.legend()
    plt.tight_layout()
plt.title("Cretin 1 zone Fast atomic kinetics timescale, gen1c")
plt.savefig(os.path.join(save_dir,f'cretin_AK_fast_timescale'+ext), transparent=False)
#%% Thermal timescale work
from scipy.integrate import trapz
gen = 'gen1c'
pressure = 7.5
params = {
    'name':'ne_cell',
    'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/thermal_timescale/default/{gen}/{pressure}torr/',
    'pressure': pressure,
    'gen': gen}
testobj = CretinData(**params)
testobj.read_data()
print(testobj.data.keys())

cv_r_dict = testobj.data['Specific HEAT vs R']
cv_dict ={}

for time, df in cv_r_dict.items():
    yval = df['cv'].astype(float)
    xval = df['r'].astype(float)
    result = trapz(yval, x=xval)/1.35
    cv_dict[time] = result

te_r_dict = testobj.data['TEV, TIV, TRADV vs R']
te_dict = {}

for time, df in te_r_dict.items():
    yval = df['tev'].astype(float)
    xval = df['r'].astype(float)
    result = trapz(yval, x=xval)/1.35
    te_dict[time] = result

heat_r_dict = testobj.data['NET HEAT vs R']

heat_dict = {}

for time, df in heat_r_dict.items():
    yval = df['heatjt'].astype(float)
    xval = df['r'].astype(float)
    result = trapz(yval, x=xval)/1.35
    heat_dict[time] = result

timescale_dict = {}

for time in te_dict.keys():
    cv      = cv_dict[time]
    te      = te_dict[time]
    heat    = heat_dict[time]
    
    timescale_dict[time] = 6.25e11*cv*te/heat
    
# plt.plot(np.array(list(timescale_dict.keys()))*1e9,np.array(list(timescale_dict.values()))*1e9,label = str(pressure)+'torr')
# plt.xlim([80,110])
# # plt.ylim([-10,300])
# plt.ylabel('timescale [ns]')
# plt.xlabel('exp time [ns]')
# plt.grid(linestyle=':', linewidth=0.5)
# plt.legend()
# plt.tight_layout()
# plt.title("Cretin Thermal timescale, gen1c")
#%%
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import numpy as np

gen = 'gen1c'
pressures = [7.5, 15, 30]

for pressure in pressures:
    params = {
        'name': 'ne_cell',
        'path': f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/thermal_timescale/default/{gen}/{pressure}torr/',
        'pressure': pressure,
        'gen': gen
    }

    testobj = CretinData(**params)
    testobj.read_data()

    dnedt_r_dict = testobj.data['dNEdt vs R']
    dnedt_dict = {}

    for time, df in dnedt_r_dict.items():
        yval = df['dnedt'].astype(float)
        xval = df['r'].astype(float)
        result = trapz(yval, x=xval) / 1.35
        dnedt_dict[time] = result

    ne_r_dict = testobj.data['NE vs R']
    ne_dict = {}

    for time, df in ne_r_dict.items():
        yval = df['ne'].astype(float)
        xval = df['r'].astype(float)
        result = trapz(yval, x=xval) / 1.35
        ne_dict[time] = result

    te_r_dict = testobj.data['TEV, TIV, TRADV vs R']
    te_dict = {}

    for time, df in te_r_dict.items():
        yval = df['tev'].astype(float)
        xval = df['r'].astype(float)
        result = trapz(yval, x=xval) / 1.35
        te_dict[time] = result

    heat_r_dict = testobj.data['NET HEAT vs R']
    heat_dict = {}

    for time, df in heat_r_dict.items():
        yval = df['heatjt'].astype(float)
        xval = df['r'].astype(float)
        result = trapz(yval, x=xval) / 1.35
        heat_dict[time] = result

    timescale_dict = {}

    for time in ne_dict.keys():
        ne = ne_dict[time]
        dnedt = dnedt_dict[time]
        te = te_dict[time]
        heat = heat_dict[time]

        inv_timescale = 2 / 3 * heat / te - 1 / ne * dnedt
        timescale_dict[time] = 1 / inv_timescale

    plt.plot(np.abs(list(timescale_dict.keys())) * 1e9, np.abs(list(timescale_dict.values())) * 1e9,
             label=str(pressure) + 'torr')

plt.xlim([80, 110])
# plt.ylim([-10,800])
plt.ylabel('timescale [ns]')
plt.xlabel('exp time [ns]')
plt.grid(linestyle=':', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.title("Cretin Thermal timescale, gen1c\n absolute value")
save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/thermal_timescale/default'
plt.savefig(os.path.join(save_dir,f'cretin_thermal_timescale_ABS_VAL'+ext), transparent=False)
