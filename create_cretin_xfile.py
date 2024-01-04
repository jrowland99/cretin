#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:15:35 2023

@author: jeff

turn a visrad radiation file into a format readable by cretin

Assume the photon energy array stays unchanged
this will form part of the header in the cretin file
"""

import numpy as np
import os 
import re
import pandas as pd
#%%gen1c

base_dir = '/home/jeff/PRISM/VISRAD_files/'
file = 'GPLe_all_C_1p0000.dat'
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

# separate_dir = '/home/jeff/PRISM/VISRAD_files/GPLe_all_C_1p0000_pandas'
# sep_file = 't{}.txt'

# for i in range(len(time_list)):
#     photon_energies = []
#     fluxes = []
#     with open(os.path.join(separate_dir,sep_file.format(i)),'r') as f:
#         lines = f.readlines()
#         lines = lines[1:]# delete first row3
#         for line in lines:
#             tokens = line.split()
#             photon_energies.append(float(tokens[0]))
#             fluxes.append(float(tokens[1]))
#     time_dict[time_list[i]] = fluxes

# def write_cretin_card(filename, header, data):
#     with open(filename, 'w') as f:
#         # Write header
#         f.write(header + '\n')

#         # Write data rows
#         for i, item in enumerate(data, start=1):
#             f.write(f'{item:11.5E}')
#             if i % 10 == 0:
#                 f.write('\n')
#             else:
#                 f.write(' ')
#         f.write('\n')


#%%
# Write the output file in the format expected by CRETIN
photon_energies = np.asarray(photon_energies)
a =1
flux_to_energy_density = a*1e7/2.42e14/3e10 #(Joules to ergs)/(ev to Hz)*c/4pi
# flux_to_energy_density = 1
with open(os.path.join(base_dir,'GPLe_all_C_xfile.txt'), 'w') as f:
    # Print the energy structure for the following radiation intensities
    f.write('c ..............................................................................\n')
    header = 'e             {}\n  '.format(len(photon_energies_b))
    # print(header)
    f.write(header)
    for i,item in enumerate(photon_energies_b,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
    f.write('\n')
    f.write('c ..............................................................................\n')
    # Now write the photon energies
    for i, time in enumerate(time_list,start=1): # list the radiation intensities for each time step
        f.write('\n')
    # These can be spatially dependent, but we will write it for boundary        node1          node2
        f.write("\nproblem  "+'  ne_cell.gen  '+' {:12.6e} '.format(time)+' {}'.format(i)+' {}'.format(1)+' {}\n\n'.format(1)) 
        fluxes = np.asarray(time_dict[time])*flux_to_energy_density
        f.write('pbcgs            1\n  ')
        for i,item in enumerate(fluxes,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
        
        f.write('\ndone\n')
        f.write('c ..............................................................................\n')
#%%
integ_fluxes = []
xaxis = np.asarray(photon_energies_b)
for i, time in enumerate(time_list):
    flux = np.asarray(time_dict[time])
    integ_fluxes.append(np.trapz(y=flux,x=xaxis))
    print(str(time)+'\t {:.5e}'.format(integ_fluxes[i]))

#%% Extra interpolation of radiation field before write xfile

def interpolate_dataframe(df, new_times):
    interpolated_df = pd.DataFrame()

    columns = df.columns.tolist()

    for i in range(len(columns) - 1):
        col1 = columns[i]
        col2 = columns[i + 1]

        time1 = float(col1)
        time2 = float(col2)

        for new_time in new_times:
            if new_time in interpolated_df.columns or new_time < time1 or new_time > time2:
                continue

            interpolated_value = df[col1] + (new_time - time1) * (df[col2] - df[col1]) / (time2 - time1)
            interpolated_df[new_time] = interpolated_value

    new_df = pd.concat([df, interpolated_df], axis=1)
    new_columns = sorted(new_df.columns.tolist(), key=float)
    new_df = new_df[new_columns]
    # print(list(new_df.columns))
    return new_df


new_times = 1e-9*(np.linspace(0,106,107))
times = []
for nt in new_times:
    nt = round(nt,12)
    if nt not in time_list:
        times.append(nt)

df = pd.DataFrame(time_dict)
interpolated_df = interpolate_dataframe(df,times)

#%% Write the output file in the format expected by CRETIN
photon_energies = np.asarray(photon_energies)
a =1
flux_to_energy_density = a*1e7/2.42e14/3e10 #(Joules to ergs)/(ev to Hz)*c/4pi
# flux_to_energy_density = 1
interp_tdict = interpolated_df.to_dict()

with open(os.path.join(base_dir,'myl_c_xfile_interp.txt'), 'w') as f:
    # Print the energy structure for the following radiation intensities
    f.write('c ..............................................................................\n')
    header = 'e             {}\n  '.format(len(photon_energies_b))
    # print(header)
    f.write(header)
    for i,item in enumerate(photon_energies_b,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
    f.write('\n')
    f.write('c ..............................................................................\n')
    # Now write the photon energies
    for i, time in enumerate(list(interp_tdict.keys()),start=1): # list the radiation intensities for each time step
        f.write('\n')
    # These can be spatially dependent, but we will write it for boundary        node1          node2
        f.write("\nproblem  "+'  ne_cell.gen  '+' {:12.6e} '.format(time)+' {}'.format(i)+' {}'.format(1)+' {}\n\n'.format(1)) 
        fluxes = np.asarray(list(interp_tdict[time].values()))*flux_to_energy_density
        f.write('pbcgs            1\n  ')
        for i,item in enumerate(fluxes,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
        
        f.write('\ndone\n')
        f.write('c ..............................................................................\n')
        
#%% Gen1F
base_dir = '/home/jeff/PRISM/Runs/Helios_runs/Gen1F_windows_only/'
save_dir = '/home/jeff/Research/Export_Control/Cretin/'
file = 'myl_f.ppd'

delim = '#\n'
with open(os.path.join(base_dir,file), 'r') as f:
    chunks = f.read().split(delim)
# # Parse the photon energies and fluxes from the file
chunks = chunks[1:-1]
time_list = []
time_dict = {}
flag = True
for chunk in chunks:
    
    photon_energies = []
    fluxes = []
    chunk = chunk.split('\n')
    line = chunk[0]
    
    if "# Dataset:" in line:
        time = round(float(line.split()[-2].strip()),1)
        time_list.append(time)
    
    chunk = chunk[3:-1] #Skip the rest of the header
    for i in range(0,len(chunk),2): # Iterate through the chunk in steps of 2, there is repeated data
        line = chunk[i]
        tokens = line.split()
        photon_energies.append(float(tokens[0]))
        fluxes.append(float(tokens[1]))
    time_dict[time]=fluxes
    if flag:
        flag=False
        photon_energies_b = photon_energies
# Write the output file in the format expected by CRETIN
photon_energies = np.asarray(photon_energies)
a =1
flux_to_energy_density = a*1e12*1e7/2.42e14/3e10 #(TW to ergs)/(ev to Hz)*c/4pi
# flux_to_energy_density = 1
with open(os.path.join(save_dir,'myl_only_all_F.txt'), 'w') as f:
    # Print the energy structure for the following radiation intensities
    f.write('c ..............................................................................\n')
    header = 'e             {}\n  '.format(len(photon_energies_b))
    # print(header)
    f.write(header)
    for i,item in enumerate(photon_energies_b,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
    f.write('\n')
    f.write('c ..............................................................................\n')
    # Now write the photon energies
    for i, time in enumerate(time_list,start=1): # list the radiation intensities for each time step
        f.write('\n')
        wtime = time*1e-9
    # These can be spatially dependent, but we will write it for boundary        node1          node2
        f.write("\nproblem  "+'  ne_cell.gen  '+' {:12.6e} '.format(wtime)+' {}'.format(i)+' {}'.format(1)+' {}\n\n'.format(1)) 
        fluxes = np.asarray(time_dict[time])*flux_to_energy_density
        f.write('pbcgs            1\n  ')
        for i,item in enumerate(fluxes,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
        
        f.write('\ndone\n')
        f.write('c ..............................................................................\n')

#%% Gen2c read & write xfile

base_dir = '/home/jeff/Research/Code_Projects/Cloudy_Projects/ss_ini/gen2c/'
file = 'SiN_only_all_C.txt'
save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/allgens_51zone/gen2c/'

delim = '#\n'
with open(os.path.join(base_dir,file), 'r') as f:
    chunks = f.read().split(delim)
# # Parse the photon energies and fluxes from the file
chunks = chunks[1:-1]
time_list = []
time_dict = {}
flag = True
for chunk in chunks:
    
    photon_energies = []
    fluxes = []
    chunk = chunk.split('\n')
    line = chunk[0]
    
    if "# Dataset:" in line:
        time = round(float(line.split()[-2].strip()),1)
        time_list.append(time)
    
    chunk = chunk[3:-1] #Skip the rest of the header
    for i in range(0,len(chunk),2): # Iterate through the chunk in steps of 2, there is repeated data
        line = chunk[i]
        tokens = line.split()
        photon_energies.append(float(tokens[0]))
        fluxes.append(float(tokens[1]))
    time_dict[time]=fluxes
    if flag:
        flag=False
        photon_energies_b = photon_energies
# Write the output file in the format expected by CRETIN
photon_energies = np.asarray(photon_energies)
a =1
flux_to_energy_density = a*1e12*1e7/2.42e14/3e10 #(TW to ergs)/(ev to Hz)*c/4pi
# flux_to_energy_density = 1
with open(os.path.join(save_dir,file), 'w') as f:
    # Print the energy structure for the following radiation intensities
    f.write('c ..............................................................................\n')
    header = 'e             {}\n  '.format(len(photon_energies_b))
    # print(header)
    f.write(header)
    for i,item in enumerate(photon_energies_b,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
    f.write('\n')
    f.write('c ..............................................................................\n')
    # Now write the photon energies
    for i, time in enumerate(time_list,start=1): # list the radiation intensities for each time step
        f.write('\n')
        wtime = time*1e-9
    # These can be spatially dependent, but we will write it for boundary        node1          node2
        f.write("\nproblem  "+'  ne_cell.gen  '+' {:12.6e} '.format(wtime)+' {}'.format(i)+' {}'.format(1)+' {}\n\n'.format(1)) 
        fluxes = np.asarray(time_dict[time])*flux_to_energy_density
        f.write('pbcgs            1\n  ')
        for i,item in enumerate(fluxes,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
        
        f.write('\ndone\n')
        f.write('c ..............................................................................\n')
#%% Gen2f read & write xfile

base_dir = '/home/jeff/Research/Code_Projects/Cloudy_Projects/ss_ini/gen2f/'
file = 'SiN_only_F_all.txt'
save_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/allgens_51zone/gen2f/'

delim = '#\n'
with open(os.path.join(base_dir,file), 'r') as f:
    chunks = f.read().split(delim)
# # Parse the photon energies and fluxes from the file
chunks = chunks[1:-1]
time_list = []
time_dict = {}
flag = True
for chunk in chunks:
    
    photon_energies = []
    fluxes = []
    chunk = chunk.split('\n')
    line = chunk[0]
    
    if "# Dataset:" in line:
        time = round(float(line.split()[-2].strip()),1)
        time_list.append(time)
    
    chunk = chunk[3:-1] #Skip the rest of the header
    for i in range(0,len(chunk),2): # Iterate through the chunk in steps of 2, there is repeated data
        line = chunk[i]
        tokens = line.split()
        photon_energies.append(float(tokens[0]))
        fluxes.append(float(tokens[1]))
    time_dict[time]=fluxes
    if flag:
        flag=False
        photon_energies_b = photon_energies
# Write the output file in the format expected by CRETIN
photon_energies = np.asarray(photon_energies)
a =1
flux_to_energy_density = a*1e12*1e7/2.42e14/3e10 #(TW to ergs)/(ev to Hz)*c/4pi
# flux_to_energy_density = 1
with open(os.path.join(save_dir,'SiN_only_all_F.txt'), 'w') as f:
    # Print the energy structure for the following radiation intensities
    f.write('c ..............................................................................\n')
    header = 'e             {}\n  '.format(len(photon_energies_b))
    # print(header)
    f.write(header)
    for i,item in enumerate(photon_energies_b,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
    f.write('\n')
    f.write('c ..............................................................................\n')
    # Now write the photon energies
    for i, time in enumerate(time_list,start=1): # list the radiation intensities for each time step
        f.write('\n')
        wtime = time*1e-9
    # These can be spatially dependent, but we will write it for boundary        node1          node2
        f.write("\nproblem  "+'  ne_cell.gen  '+' {:12.6e} '.format(wtime)+' {}'.format(i)+' {}'.format(1)+' {}\n\n'.format(1)) 
        fluxes = np.asarray(time_dict[time])*flux_to_energy_density
        f.write('pbcgs            1\n  ')
        for i,item in enumerate(fluxes,start=1):
            f.write(f'{item:11.5E}')
            if i % 10 == 0:
                f.write('\n  ')
            else:
                f.write(' ')
        
        f.write('\ndone\n')
        f.write('c ..............................................................................\n')