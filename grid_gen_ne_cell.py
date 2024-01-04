#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 11:00:22 2023

@author: jeff

write cretin model for all cases of neon cell
"""
import numpy as np
import os
import shutil
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt
import sys
import subprocess
import time as Time
from queue import Queue
from threading import Thread
sys.path.append("/home/jeff/Research/Code_Projects/git/main/cretin/")
sys.path.append("/home/jeff/Research/Code_Projects/git/main/cloudy/")
sys.path.append("/home/jeff/Research/Code_Projects/git/main/prism/")
from CretinData import *


def cretin_thread(directory,name): # directory points to folder containing .gen file, name is "run.gen"
    os.chdir(directory)
    print(f"Inside {directory}")
    with open(os.devnull,'w') as f:
        full_command = "./cretin "+name
        t1 = Time.time()
        subprocess.call(full_command,shell=True,stdout=f)

        t2 = Time.time()
        print (f'{directory}\n'+"finish in {} seconds.\n".format((t2-t1)))
        
def run_cretin(dir_list,num_threads):
    command_queue = Queue()
    com_list = []
    # Walk through all directories and run all '.in' files
    for directory in dir_list:
        for root, dirs, files in os.walk(directory):
            # print(dirs)
            for filename in files:
                if filename.endswith('.gen'):
                    name = filename
                    com_list.append([root,name])
        
     
    def worker():
        while True:
            command = command_queue.get()
            print('starting {}'.format(str(command[0][0])))
            cretin_thread(directory=str(command[0][0]),name = str(command[0][1]))
            command_queue.task_done()
            
    for i in range(num_threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()
    	
    count = 0
    for c in com_list:
        command_queue.put([c,count])
        count += 1
    	
    command_queue.join()
    

#%%
base_dir    ='/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v13'
com         = 'c ------------------------------------------------------------\n'
pressures   = [7.5,15,30]# using 2.5e17,5e17,1e18 as densities
densities   = [2.5e17,5e17,1e18]
gens        = ['gen1c','gen2c','gen2f']
g2xfile     = {'gen1c':'myl_c_xfile_interp.txt','gen1f':'myl_only_all_F.txt',  'gen2c':'SiN_only_all_C.txt','gen2f':'SiN_only_all_F.txt'}
term        = '/home/jeff/Research/Export_Control/Cretin/xfiledir/term10.dat'
cretin_ex   = '/home/jeff/Research/Export_Control/Cretin/cretin'
onezone     = False
rad         = False
use_escape  = True

# gens = ['gen1f']
# pressures = [15]
# densities = [5e17]
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)
for gen in gens:
    directory = os.path.join(base_dir,gen)
    xfile = g2xfile[gen]
    if not os.path.isdir(directory):
        os.mkdir(directory)
    for p in range(len(pressures)):
        dens = densities[p]
        runpath = os.path.join(directory,str(pressures[p])+"torr")
        if not os.path.isdir(runpath):
            os.mkdir(runpath)
        os.chdir(runpath)
        fname = 'ne_cell.gen'
        shutil.copy(term, runpath)
        shutil.copy(cretin_ex, runpath)
        f = open(fname,"w")
        # f.write("\n")
        f.write("c        *** Neon Slab No Windows ***\n")
        f.write(f"alias N_neon	{dens}\t\t\t! neon number density\n")
        if onezone:
            f.write("alias Nnode	1\n")
        else:
            f.write("alias Nnode	51\n")
        f.write("\n")
        f.write("alias T0	0.024\n") # 276K in eV
        f.write("alias T1	0.1\n")
        f.write("\n")
        f.write("alias SIZE	1.35\n")
        f.write("alias DTMIN 	1e-12\n")
        f.write("alias DTMAX 	1e-7\n")
        f.write("\n")
        # f.write("alias N0        1\n")
        # f.write("alias N1        Nnode+/2\n")
        # f.write("alias N2        Nnode\n")
        # f.write("alias DN        25\n")
        # f.write("\n")
        # f.write("alias R0        0.\n")
        # f.write("alias R1        R0 + SIZE/2.\n")
        # f.write("alias R2        R0 + SIZE\n")
        # f.write("alias DR        0.01\n")
        # f.write("\n")
        # f.write("alias S31       1\n")
        # f.write("alias S44       10\n")
        f.write("\n")
        f.write(com)
        f.write("c   Materials\n")
        f.write(com)
        f.write("atoms hydrogenic ne\n")
        f.write("  modeltype dca term\n")
        f.write("\n")
        f.write("region  1 Nnode  T0		! Neon Slab\n")
        f.write("  element  1  N_neon\n")
        f.write(com)
        f.write("c   Geometry\n")
        f.write(com)
        f.write("geometry slab\n")
        
        if not onezone:
            f.write("rlin 1 Nnode 0 SIZE 		! Equally spaced zones\n")
        f.write(com)
        f.write("c   Radiation\n")
        f.write(com)
        f.write("angles 3\n")
        f.write("\n")
        f.write("alias emin  1.e-1\n")
        f.write("alias emax  4.e3\n")
        f.write("\n")
        f.write(f"xfile ix=12 {xfile}   ! define xfile 12\n")
        f.write("\n")
        f.write("source jbndry 12 pbins 1 1.          ! define jbndry from xfile 12\n")
        f.write("xfilebc 12 -1 0 1. 1                  ! use jbndry from xfile 12\n")
        f.write(com)
        f.write("c   Controls\n")
        f.write(com)
        f.write("tstart 0.\n")
        f.write("tquit  1.08e-07\n")
        f.write("restart\n")
        f.write(com)
        f.write("c   Switches and Parameters\n")
        f.write(com)

        if rad:
            f.write("switch 100 1\n")
        else:
            f.write("switch 31  1                    ! do temperature calculation TD = 1, SS = -1,rad hydro = 4\n")
            f.write("switch 36 1			! radiation transport on = 1, off = 0\n")
        if use_escape:
            f.write("switch 33 1			! use excape factors for all photoexcitations if >0\n")
        if not onezone:
            f.write("switch 2 0			! hydrodynamics 0=off\n")
        f.write("switch 11  1                    ! make .plt file\n")
        f.write("switch 20  1                  	! calculate NLTE populations using rate matrices\n")
        # f.write("switch 21  1                    ! account for transition energies\n")
        f.write("switch 22  1                    ! recalculate destruction rates on first iteration\n")
        f.write("switch 25  1                    ! time dependent kinetics = 1, SS = 0\n")
        f.write("switch 28  2                    ! 2: initialize in steady-state kinetics w/o radiation transfer\n")
        f.write("switch 29  2                    ! use variable timesteps\n")
        f.write("switch 30 20                    ! dump every n timesteps\n")
        f.write("switch 44 S44                   ! max iterations per timestep\n")
        # f.write("switch 49 1                   ! Electron thermal conduction on = 1\n")
        # f.write("switch 51 0                   ! Electron thermal conduction coeff Spitzer-Harm")
        f.write("\n")
        f.write("param  40 DTMAX 		! time between edits\n")
        f.write("param  41 DTMIN 		! initial timestep\n")
        f.write("param  44 DTMIN 		! minimum timestep\n")
        f.write("param  45 DTMAX 		! maximum timestep\n")
        f.write("param  42  0.1                 	! max change in Zbar\n")
        f.write("param  46  0.05                 ! maximum frac. change in temperature\n")
        f.write("param  61  1.e-8               	! iso-sequence population threshold\n")
        f.write("param  102  T1              	! min Te\n")
        f.write(com)
        f.write("c   Edits\n")
        f.write(com)

        f.write("\n")
        
        f.write("plot \"TEV, TIV, TRADV vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tev\n")
        f.write("  yvar  tiv\n")
        f.write("  yvar  tradv\n")
        f.write("\n")
        
        f.write("plot \"NE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  ne\n")
        f.write("\n")
        
        f.write("plot \"dNEdt vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  dnedt \n")
        f.write("\n")
        
        f.write("plot \"ZBAR vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  zbar\n")
        f.write("\n")
        
        f.write("plot \"ERAD vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  erad\n")
        f.write("\n")
        
        f.write("plot \"TAUEE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauee\n")
        f.write("\n")
        
        f.write("plot \"TAUEI vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauei\n")
        f.write("\n")
        
        f.write("plot \"TAUII vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauii\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("xvar r\n")
        f.write("yvar heatt\n")
        f.write("\n")        
        
        f.write("plot \"GAMMATOT vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  gammatot 1 0 0:5 1\n")
        f.write("\n")

        f.write("plot \"YISO vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 0:5\n")
        f.write("\n")
        
        f.write("plot \"GAMMATOT vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  gammatot 1 0 6:10 1\n")
        f.write("\n")

        f.write("plot \"YISO vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 6:10\n")
        f.write("\n")
        

        f.write("plot \"JNU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jnu   0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"TAU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"TAU vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        # f.write("plot \"EIGENVAL vs TIME\"\n")
        # f.write("  xvar  time\n")
        # f.write("  yvar  eigenval 1 1:-1:DN 0 1:4\n")
        # f.write("\n")
        
        # f.write("plot \"TEV, TIV, TRADV vs TIME\"\n")
        # f.write("  xvar  time\n")
        # f.write("  yvar  tev   0 1:-1:DN\n")
        # f.write("  yvar  tiv   0 1:-1:DN\n")
        # f.write("  yvar  tradv 0 1:-1:DN\n")
        # f.write("\n")
        
        f.write("plot \"NE vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  ne 0 1:-1:DN\n")
        f.write("\n")

        f.write("plot \"ZBAR vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  zbar 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"ERAD vs TIME\"\n")
        f.write(" xvar time\n")
        f.write(" yvar erad 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"GAMMATOT vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  gammatot 1 0 0:10 1\n")
        f.write("\n")
        
        f.write("plot \"EIGENVAL vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  eigenv_r 1 1:51 0 -1:\n")
        f.write("\n")
        
        # f.write("plot \"EIGENVAL vs TIME desc\"\n")
        # f.write("  xvar  time\n")
        # f.write("  yvar  eigenv_r 1 1 0 1:2:1\n")
        # f.write("\n")
        
        f.write("plot \"TEV, TIV vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  tev 0 1\n")
        f.write("  yvar  tiv 0 1 \n")
        f.write("\n")
        
        f.write("plot \"YISO vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar yisofrac 1 1 0:10\n")
        f.write("\n")

        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar ntry\n")
        f.write("\n")
        
        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar dtime\n")
        f.write("\n")
        
        f.write("#ifdef DISPLAY\n")
        f.write("display 1\n")
        f.write("display 2\n")
        f.write("display 3\n")
        f.write("display 4\n")
        f.write("display 5\n")
        f.write("#endif\n")
        f.close()  
# #%%    
# # base_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/all_gen_v2/'    
# run_cretin([base_dir],9)        
#%% Constant TD Temp run
base_dir    ='/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/const_td_temp/'
com         = 'c ------------------------------------------------------------\n'
pressures   = [7.5,15,30]# using 2.5e17,5e17,1e18 as densities
densities   = [2.5e17,5e17,1e18]
gens        = ['gen1c']#['gen1c','gen2c','gen2f']
g2xfile     = {'gen1c':'myl_c_xfile_interp.txt','gen1f':'myl_only_all_F.txt',  'gen2c':'SiN_only_all_C.txt','gen2f':'SiN_only_all_F.txt'}
term        = '/home/jeff/Research/Export_Control/Cretin/xfiledir/term10.dat'
cretin_ex   = '/home/jeff/Research/Export_Control/Cretin/cretin'
onezone     = False
rad         = False
use_escape  = True

temperature_dir = '/home/jeff/Research/Export_Control/Cretin/Neon_Cell/td_latest/'
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)
for gen in gens:
    directory = os.path.join(base_dir,gen)
    xfile = g2xfile[gen]
    if not os.path.isdir(directory):
        os.mkdir(directory)
    for p in range(len(pressures)):
        dens = densities[p]
        runpath = os.path.join(directory,str(pressures[p])+"torr")
        if not os.path.isdir(runpath):
            os.mkdir(runpath)
        os.chdir(runpath)
        # get temperature for history
        td_temp_path = os.path.join(temperature_dir,gen,str(pressures[p])+'torr')
        params = {
            'name':'ne_cell',
            # 'path':f'/home/jeff/Research/Export_Control/Cretin/Neon_Cell/eigenvals/v7/',
            'path':td_temp_path,
            'pressure': pressures[p],
            'gen': gen}
        obj = CretinData(**params)
        obj.read_data()
        td_temps = obj.time_data['TEV']
        td_times = obj.time_data['time']
        
        fname = 'ne_cell.gen'
        shutil.copy(term, runpath)
        shutil.copy(cretin_ex, runpath)
        f = open(fname,"w")
        # f.write("\n")
        f.write("c        *** Neon Slab No Windows ***\n")
        f.write(f"alias N_neon	{dens}\t\t\t! neon number density\n")
        if onezone:
            f.write("alias Nnode	1\n")
        else:
            f.write("alias Nnode	51\n")
        f.write("\n")
        f.write("alias T0	0.024\n") # 276K in eV
        f.write("alias T1	0.1\n")
        f.write("\n")
        f.write("alias SIZE	1.35\n")
        f.write("alias DTMIN 	1e-10\n")
        f.write("alias DTMAX 	1e-9\n")
        f.write("\n")
        f.write("alias N0        1\n")
        f.write("alias N1        Nnode+/2\n")
        f.write("alias N2        Nnode\n")
        f.write("alias DN        25\n")
        f.write("\n")
        f.write("alias R0        0.\n")
        f.write("alias R1        R0 + SIZE/2.\n")
        f.write("alias R2        R0 + SIZE\n")
        f.write("alias DR        0.01\n")
        f.write("\n")
        f.write("alias S31       1\n")
        f.write("alias S44       10\n")
        f.write("\n")
        f.write(com)
        f.write("c   Materials\n")
        f.write(com)
        f.write("atoms hydrogenic ne\n")
        f.write("  modeltype dca term\n")
        f.write("\n")
        f.write("region  1 Nnode  T0		! Neon Slab\n")
        f.write("  element  1  N_neon\n")
        f.write(com)
        f.write("c   Geometry\n")
        f.write(com)
        f.write("geometry slab\n")
        
        if not onezone:
            f.write("rlin 1 Nnode 0 SIZE 		! Equally spaced zones\n")
        f.write(com)
        f.write("c   Radiation\n")
        f.write(com)
        f.write("angles 3\n")
        f.write("\n")
        f.write("alias emin  1.e-1\n")
        f.write("alias emax  4.e3\n")
        f.write("\n")
        f.write(f"xfile ix=12 {xfile}   ! define xfile 12\n")
        f.write("\n")
        f.write("source jbndry 12 pbins 1 1.          ! define jbndry from xfile 12\n")
        f.write("xfilebc 12 -1 0 1. 1                  ! use jbndry from xfile 12\n")
        f.write(com)
        # Constant temperature source command
        f.write("source te value history 15 \n")
        f.write("history 15 \n")
        for i in range(len(td_times)):
            f.write(f"tv {td_times[i]} {td_temps[i]}\n")
        f.write(com)
        f.write("c   Controls\n")
        f.write(com)
        f.write("tstart 0.\n")
        f.write("tquit  1.08e-07\n")
        f.write("restart\n")
        f.write(com)
        f.write("c   Switches and Parameters\n")
        f.write(com)

        if rad:
            f.write("switch 100 1\n")
        else:
            f.write("switch 31  0                    ! do temperature calculation TD = 1, SS = -1, off = 0\n")
            f.write("switch 36 1			! radiation transport on = 1, off = 0\n")
        if use_escape:
            f.write("switch 33 1			! use excape factors for all photoexcitations if >0\n")
        if not onezone:
            f.write("switch 2 0			! hydrodynamics 0=off\n")
        f.write("switch 11  1                    ! make .plt file\n")
        f.write("switch 20  1                  	! calculate NLTE populations using rate matrices\n")
        # f.write("switch 21  1                    ! account for transition energies\n")
        f.write("switch 22  1                    ! recalculate destruction rates on first iteration\n")
        f.write("switch 25  0                    ! time dependent kinetics = 1, SS = 0\n")
        f.write("switch 28  2                    ! 2: initialize in steady-state kinetics w/o radiation transfer\n")
        f.write("switch 29  2                    ! use variable timesteps\n")
        f.write("switch 30 20                    ! dump every n timesteps\n")
        f.write("switch 44 S44                   ! max iterations per timestep\n")
        # f.write("switch 49 1                   ! Electron thermal conduction on = 1\n")
        # f.write("switch 51 0                   ! Electron thermal conduction coeff Spitzer-Harm")
        f.write("\n")
        f.write("param  40 DTMAX 		! time between edits\n")
        f.write("param  41 DTMIN 		! initial timestep\n")
        f.write("param  44 DTMIN 		! minimum timestep\n")
        f.write("param  45 DTMAX 		! maximum timestep\n")
        f.write("param  42  0.1                 	! max change in Zbar\n")
        f.write("param  46  0.05                 ! maximum frac. change in temperature\n")
        f.write("param  61  1.e-8               	! iso-sequence population threshold\n")
        f.write("param  102  T1              	! min Te\n")
        f.write(com)
        f.write("c   Edits\n")
        f.write(com)

        f.write("\n")
        
        f.write("plot \"TEV, TIV, TRADV vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tev\n")
        f.write("  yvar  tiv\n")
        f.write("  yvar  tradv\n")
        f.write("\n")
        
        f.write("plot \"NE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  ne\n")
        f.write("\n")
        
        f.write("plot \"dNEdt vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  dnedt \n")
        f.write("\n")
        
        f.write("plot \"ZBAR vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  zbar\n")
        f.write("\n")
        
        f.write("plot \"ERAD vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  erad\n")
        f.write("\n")
        
        f.write("plot \"TAUEE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauee\n")
        f.write("\n")
        
        f.write("plot \"TAUEI vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauei\n")
        f.write("\n")
        
        f.write("plot \"TAUII vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauii\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("xvar r\n")
        f.write("yvar heatt\n")
        f.write("\n")        
        
        f.write("plot \"GAMMATOT vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  gammatot 1 0 0:5 1\n")
        f.write("\n")

        f.write("plot \"YISO vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 0:5\n")
        f.write("\n")
        
        f.write("plot \"GAMMATOT vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  gammatot 1 0 6:10 1\n")
        f.write("\n")

        f.write("plot \"YISO vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 6:10\n")
        f.write("\n")
        

        f.write("plot \"JNU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jnu   0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"TAU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"TAU vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        # f.write("plot \"EIGENVAL vs TIME\"\n")
        # f.write("  xvar  time\n")
        # f.write("  yvar  eigenval 1 1:-1:DN 0 1:4\n")
        # f.write("\n")
        
        # f.write("plot \"TEV, TIV, TRADV vs TIME\"\n")
        # f.write("  xvar  time\n")
        # f.write("  yvar  tev   0 1:-1:DN\n")
        # f.write("  yvar  tiv   0 1:-1:DN\n")
        # f.write("  yvar  tradv 0 1:-1:DN\n")
        # f.write("\n")
        
        f.write("plot \"NE vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  ne 0 1:-1:DN\n")
        f.write("\n")

        f.write("plot \"ZBAR vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  zbar 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"ERAD vs TIME\"\n")
        f.write(" xvar time\n")
        f.write(" yvar erad 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"GAMMATOT vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  gammatot 1 0 0:10 1\n")
        f.write("\n")
        
        f.write("plot \"EIGENVAL vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  eigenval_r 1 0 0:10 1:3\n")
        f.write("\n")
        
        f.write("plot \"TEV, TIV vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  tev 0 1\n")
        f.write("  yvar  tiv 0 1 \n")
        f.write("\n")
        
        f.write("plot \"YISO vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar yisofrac 1 1 0:10\n")
        f.write("\n")

        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar ntry\n")
        f.write("\n")
        
        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar dtime\n")
        f.write("\n")
        
        f.write("#ifdef DISPLAY\n")
        f.write("display 1\n")
        f.write("display 2\n")
        f.write("display 3\n")
        f.write("display 4\n")
        f.write("display 5\n")
        f.write("#endif\n")
        f.close()        
#%% Spectrum option
base_dir    ='/home/jeff/Research/Export_Control/Cretin/Neon_Cell/thermal_timescale/rad_hydro/'
com         = 'c ------------------------------------------------------------\n'
pressures   = [7.5,15,30]# using 2.5e17,5e17,1e18 as densities
densities   = [2.5e17,5e17,1e18]
gens        = ['gen1c']
g2xfile     = {'gen1c':'myl_c_xfile_interp.txt','gen1f':'myl_only_all_F.txt',  'gen2c':'SiN_only_all_C.txt','gen2f':'SiN_only_all_F.txt'}
term        = '/home/jeff/Research/Export_Control/Cretin/xfiledir/term10.dat'
cretin_ex   = '/home/jeff/Research/Export_Control/Cretin/cretin'
onezone     = False
rad         = False
use_escape  = True

# gens = ['gen1f']
# pressures = [15]
# densities = [5e17]
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)
for gen in gens:
    directory = os.path.join(base_dir,gen)
    xfile = g2xfile[gen]
    if not os.path.isdir(directory):
        os.mkdir(directory)
    for p in range(len(pressures)):
        dens = densities[p]
        runpath = os.path.join(directory,str(pressures[p])+"torr")
        if not os.path.isdir(runpath):
            os.mkdir(runpath)
        os.chdir(runpath)
        fname = 'ne_cell.gen'
        shutil.copy(term, runpath)
        shutil.copy(cretin_ex, runpath)
        f = open(fname,"w")
        # f.write("\n")
        f.write("c        *** Neon Slab No Windows ***\n")
        f.write(f"alias N_neon	{dens}\t\t\t! neon number density\n")
        if onezone:
            f.write("alias Nnode	1\n")
        else:
            f.write("alias Nnode	51\n")
        f.write("\n")
        f.write("alias T0	0.024\n") # 276K in eV
        f.write("alias T1	0.1\n")
        f.write("\n")
        f.write("alias SIZE	1.35\n")
        f.write("alias DTMIN 	1e-12\n")
        f.write("alias DTMAX 	1e-7\n")
        f.write("\n")
        # f.write("alias N0        1\n")
        # f.write("alias N1        Nnode+/2\n")
        # f.write("alias N2        Nnode\n")
        # f.write("alias DN        25\n")
        # f.write("\n")
        # f.write("alias R0        0.\n")
        # f.write("alias R1        R0 + SIZE/2.\n")
        # f.write("alias R2        R0 + SIZE\n")
        # f.write("alias DR        0.01\n")
        f.write("\n")
        f.write("alias S31       1\n")
        f.write("alias S44       10\n")
        f.write("\n")
        f.write(com)
        f.write("c   Materials\n")
        f.write(com)
        f.write("atoms hydrogenic ne\n")
        f.write("  modeltype dca term\n")
        f.write("\n")
        f.write("region  1 Nnode  T0		! Neon Slab\n")
        f.write("  element  1  N_neon\n")
        f.write(com)
        f.write("c   Geometry\n")
        f.write(com)
        f.write("geometry slab\n")
        
        if not onezone:
            f.write("rlin 1 Nnode 0 SIZE 		! Equally spaced zones\n")
        f.write(com)
        f.write("c   Radiation\n")
        f.write(com)
        f.write("angles 3\n")
        f.write("\n")
        f.write("alias emin  1.e-1\n")
        f.write("alias emax  4.e3\n")
        f.write("\n")
        f.write(f"xfile ix=12 {xfile}   ! define xfile 12\n")
        f.write("\n")
        f.write("source jbndry 12 pbins 1 1.          ! define jbndry from xfile 12\n")
        f.write("xfilebc 12 -1 0 1. 1                  ! use jbndry from xfile 12\n")
        f.write("background 0. 1.e15 ! add background electrons to fix sw 31 = 4\n")
        f.write("ebins 1e3 emin emax\n")
        
        f.write(com)
        f.write("c   Controls\n")
        f.write(com)
        f.write("tstart 0.\n")
        f.write("tquit  1.08e-07\n")
        f.write("restart\n")
        # f.write("nltedump spectrum transmission\n")
        f.write(com)
        f.write("c   Switches and Parameters\n")
        f.write(com)

        if rad:
            f.write("switch 100 1\n")
        else:
            f.write("switch 31  4                    ! do temperature calculation TD = 1, SS = -1,rad hydro = 4\n")
            f.write("switch 36 1			! radiation transport on = 1, off = 0\n")
        if use_escape:
            f.write("switch 33 1			! use excape factors for all photoexcitations if >0\n")
        if not onezone:
            f.write("switch 2 1			! hydrodynamics 0 = off 1 = on\n")
        f.write("switch 11  1                    ! make .plt file\n")
        f.write("switch 20  1                  	! calculate NLTE populations using rate matrices\n")
        # f.write("switch 21  1                    ! account for transition energies\n")
        f.write("switch 22  1                    ! recalculate destruction rates on first iteration\n")
        f.write("switch 25  1                    ! time dependent kinetics = 1, SS = 0\n")
        f.write("switch 28  2                    ! 2: initialize in steady-state kinetics w/o radiation transfer\n")
        f.write("switch 29  2                    ! use variable timesteps\n")
        f.write("switch 30 20                    ! dump every n timesteps\n")
        f.write("switch 37 1                     ! do line transfer\n")
        f.write("switch 44 S44                   ! max iterations per timestep\n")
        # f.write("switch 49 1                   ! Electron thermal conduction on = 1\n")
        # f.write("switch 51 0                   ! Electron thermal conduction coeff Spitzer-Harm")
        f.write("\n")
        f.write("param  40 DTMAX 		! time between edits\n")
        f.write("param  41 DTMIN 		! initial timestep\n")
        f.write("param  44 DTMIN 		! minimum timestep\n")
        f.write("param  45 DTMAX 		! maximum timestep\n")
        f.write("param  42  0.1                 	! max change in Zbar\n")
        f.write("param  46  0.05                 ! maximum frac. change in temperature\n")
        f.write("param  61  1.e-8               	! iso-sequence population threshold\n")
        f.write("param  102  T1              	! min Te\n")
        f.write(com)
        f.write("c   Edits\n")
        f.write(com)

        f.write("\n")
        
        f.write("plot \"TEV, TIV, TRADV vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tev\n")
        f.write("  yvar  tiv\n")
        f.write("  yvar  tradv\n")
        f.write("\n")
        
        f.write("plot \"NE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  ne\n")
        f.write("\n")
        
        f.write("plot \"dNEdt vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  dnedt \n")
        f.write("\n")

        
        f.write("plot \"ERAD vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  erad\n")
        f.write("\n")
        
        f.write("plot \"TAUEE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauee\n")
        f.write("\n")
        
        f.write("plot \"TAUEI vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauei\n")
        f.write("\n")
        
        f.write("plot \"TAUII vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauii\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("xvar r\n")
        f.write("yvar heatt\n")
        f.write("\n")        
        

        f.write("plot \"YISO vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 0:5\n")
        f.write("\n")


        f.write("plot \"YISO vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 6:10\n")
        f.write("\n")
        

        f.write("plot \"JNU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jnu 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"TAU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"JSP vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jsp 0 1 0\n")
        f.write("\n")
        
        # f.write("plot \"SPECTRUM vs ENERGY\"\n")
        # f.write("  xvar sp_energy\n")
        # f.write("  yvar jsp 0 1\n")
        # f.write("  yvar jsp 0 51\n")

        f.write("plot \"JSP vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jsp 0 1 0\n")
        f.write("\n")        

        f.write("plot \"TAU vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  heatjt\n")
        f.write("\n")
        
        f.write("plot \"NE vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  ne 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"ERAD vs TIME\"\n")
        f.write(" xvar time\n")
        f.write(" yvar erad 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"EIGENVAL vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  eigenv_r 1 1:51 0 -1:\n")
        f.write("\n")

        f.write("plot \"TEV, TIV vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  tev 0 1\n")
        f.write("  yvar  tiv 0 1 \n")
        f.write("\n")
        
        f.write("plot \"YISO vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar yisofrac 1 1 0:10\n")
        f.write("\n")

        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar ntry\n")
        f.write("\n")
        
        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar dtime\n")
        f.write("\n")
        
        f.write("#ifdef DISPLAY\n")
        f.write("display 1\n")
        f.write("display 2\n")
        f.write("display 3\n")
        f.write("display 4\n")
        f.write("display 5\n")
        f.write("#endif\n")
        f.close()  
#%% Lines w/ rad hydro 
base_dir    ='/home/jeff/Research/Export_Control/Cretin/Neon_Cell/rad_hydro_w_lines/v1/'
com         = 'c ------------------------------------------------------------\n'
pressures   = [7.5,15,30]# using 2.5e17,5e17,1e18 as densities
densities   = [2.5e17,5e17,1e18]
gens        = ['gen1c']
g2xfile     = {'gen1c':'myl_c_xfile_interp.txt','gen1f':'myl_only_all_F.txt',  'gen2c':'SiN_only_all_C.txt','gen2f':'SiN_only_all_F.txt'}
term        = '/home/jeff/Research/Export_Control/Cretin/xfiledir/term10.dat'
cretin_ex   = '/home/jeff/Research/Export_Control/Cretin/cretin'
onezone     = False
rad         = False
use_escape  = True

# gens = ['gen1f']
# pressures = [15]
# densities = [5e17]
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)
for gen in gens:
    directory = os.path.join(base_dir,gen)
    xfile = g2xfile[gen]
    if not os.path.isdir(directory):
        os.mkdir(directory)
    for p in range(len(pressures)):
        dens = densities[p]
        runpath = os.path.join(directory,str(pressures[p])+"torr")
        if not os.path.isdir(runpath):
            os.mkdir(runpath)
        os.chdir(runpath)
        fname = 'ne_cell.gen'
        shutil.copy(term, runpath)
        shutil.copy(cretin_ex, runpath)
        f = open(fname,"w")
        # f.write("\n")
        f.write("c        *** Neon Slab No Windows ***\n")
        f.write(f"alias N_neon	{dens}\t\t\t! neon number density\n")
        if onezone:
            f.write("alias Nnode	1\n")
        else:
            f.write("alias Nnode	51\n")
        f.write("\n")
        f.write("alias T0	0.024\n") # 276K in eV
        f.write("alias T1	0.1\n")
        f.write("\n")
        f.write("alias SIZE	1.35\n")
        f.write("alias DTMIN 	1e-10\n")
        f.write("alias DTMAX 	1e-7\n")
        # f.write("\n")
        # f.write("alias N0        1\n")
        # f.write("alias N1        Nnode+/2\n")
        # f.write("alias N2        Nnode\n")
        # f.write("alias DN        25\n")
        # f.write("\n")
        # f.write("alias R0        0.\n")
        # f.write("alias R1        R0 + SIZE/2.\n")
        # f.write("alias R2        R0 + SIZE\n")
        # f.write("alias DR        0.01\n")
        # f.write("\n")
        # f.write("alias S31       1\n")
        # f.write("alias S44       10\n")
        f.write("\n")
        f.write(com)
        f.write("c   Materials\n")
        f.write(com)
        f.write("atoms hydrogenic ne\n")
        f.write("  modeltype dca term\n")
        f.write("\n")
        f.write("region  1 Nnode  T0		! Neon Slab\n")
        f.write("  element  1  N_neon\n")
        f.write(com)
        f.write("c   Geometry\n")
        f.write(com)
        f.write("geometry slab\n")
        
        if not onezone:
            f.write("rlin 1 Nnode 0 SIZE 		! Equally spaced zones\n")
        f.write(com)
        f.write("c   Radiation\n")
        f.write(com)
        f.write("angles 3\n")
        f.write("\n")
        f.write("alias emin  1.e-1\n")
        f.write("alias emax  4.e3\n")
        
        f.write("ebins 1e3 emin emax\n")
        f.write("spectrum 50000 .01 900 1\n")
        f.write("spectrum 500 900 1300 1\n")
        
        f.write("spectrum 500 1300 4000 1\n")
        f.write("\n")
        f.write(f"xfile ix=12 {xfile}   ! define xfile 12\n")
        f.write("\n")
        f.write("source jbndry 12 pbins 1 1.          ! define jbndry from xfile 12\n")
        f.write("xfilebc 12 -1 0 1. 1                  ! use jbndry from xfile 12\n")
        f.write(com)
        f.write("\n")
        f.write("linedefault crd \n")
        f.write("\n")
        f.write("line 1 1 1 1 1 2     \t! H-like Ly-alpha\n")
        f.write("    lbins 25 5. 1.02\n")
        f.write("line 2 1 2 1 2 2       \t! He-like Ly-alpha\n")
        f.write("    lbins 25 5. 1.02\n")
        f.write("line 3 1 3 1 3 2       \t! Li-like Ly-alpha\n")
        f.write("    lbins 50 10 1.02\n")
        f.write("line 4 1 2 1 2 3       \t! He-like Ly-beta\n")
        f.write("    lbins 25 5. 1.02\n")
        f.write("line 5 1 2 1 2 4       \t! He-like Ly-gamma\n")
        f.write("    lbins 25 5. 1.02\n")
        f.write("line 6 1 3 1 3 3       \t! Li-like Ly-beta\n")
        f.write("    lbins 50 10 1.02\n")
        f.write(com)
        f.write("c   Controls\n")
        f.write(com)
        f.write("tstart 0.\n")
        f.write("tquit  1.08e-07\n")
        f.write("restart\n")
        f.write("nltedump spectrum transmission\n")
        f.write(com)
        f.write("c   Switches and Parameters\n")
        f.write(com)
        f.write("switch 2 0			! hydrodynamics 0=off\n")
        f.write("switch 11  1                    ! make .plt file\n")
        f.write("switch 20  1                  	! calculate NLTE populations using rate matrices\n")
        # f.write("switch 21  1                    ! account for transition energies\n")
        f.write("switch 22  1                    ! recalculate destruction rates on first iteration\n")
        f.write("switch 25  1                    ! time dependent kinetics = 1, SS = 0\n")
        f.write("switch 28  2                    ! 2: initialize in steady-state kinetics w/o radiation transfer\n")
        f.write("switch 29  2                    ! use variable timesteps\n")
        f.write("switch 30 10                    ! dump every n timesteps\n")
        f.write("switch 31  1                    ! do temperature calculation TD = 1, SS = -1,rad hydro = 4\n")
        f.write("switch 33 1			! use excape factors for all photoexcitations if >0\n")
        f.write("switch 36 1			         ! radiation transport on = 1, off = 0\n")
        f.write("switch 37 1                     ! do line transfer = 1\n")
        f.write("switch 44 S44                   ! max iterations per timestep\n")
        # f.write("switch 49 1                   ! Electron thermal conduction on = 1\n")
        # f.write("switch 51 0                   ! Electron thermal conduction coeff Spitzer-Harm")
        f.write("switch 52 -1                    ! do Stark broadening for everything\n")
        f.write("switch 55 1                     ! do continuum lowering\n")
        f.write("\n")
        f.write("param  40 DTMAX 		! time between edits\n")
        f.write("param  41 DTMIN 		! initial timestep\n")
        f.write("param  44 DTMIN 		! minimum timestep\n")
        f.write("param  45 DTMAX 		! maximum timestep\n")
        f.write("param  42  0.1                 	! max change in Zbar\n")
        f.write("param  46  0.05                 ! maximum frac. change in temperature\n")
        f.write("param  61  1.e-8               	! iso-sequence population threshold\n")
        f.write("param  102  T1              	! min Te\n")
        f.write(com)
        f.write("c   Edits\n")
        f.write(com)

        f.write("\n")
        
        f.write("plot \"TEV, TIV, TRADV vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tev\n")
        f.write("  yvar  tiv\n")
        f.write("  yvar  tradv\n")
        f.write("\n")
        
        f.write("plot \"NE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  ne\n")
        f.write("\n")
        
        f.write("plot \"dNEdt vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  dnedt \n")
        f.write("\n")

        
        f.write("plot \"ERAD vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  erad\n")
        f.write("\n")
        
        f.write("plot \"TAUEE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauee\n")
        f.write("\n")
        
        f.write("plot \"TAUEI vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauei\n")
        f.write("\n")
        
        f.write("plot \"TAUII vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauii\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("xvar r\n")
        f.write("yvar heatt\n")
        f.write("\n")        
        

        f.write("plot \"YISO vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 0:5\n")
        f.write("\n")


        f.write("plot \"YISO vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 6:10\n")
        f.write("\n")
        

        f.write("plot \"JNU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jnu 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"TAU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"JBAR vs TIME\"\n")
        # f.write("averaged\n")
        f.write("  xvar time\n")
        f.write("  yvar jbarcgs 1 1\n")
        f.write("  yvar jbarcgs 2 1\n")
        f.write("  yvar jbarcgs 3 1\n")
        f.write("\n")
        
        f.write("plot \"JSAT vs TIME\"\n")
        # f.write("averaged\n")
        f.write("  xvar time\n")
        f.write("  yvar jsat 1 1\n")
        f.write("  yvar jsat 2 1\n")
        f.write("  yvar jsat 3 1\n")
        f.write("\n")
        
        # f.write("plot \"TAULINE vs EVLINE 1\"\n")
        # f.write("  xvar evline 1\n")
        # f.write("  yvar tauline 1 Nnode \n")
        # f.write("\n")
        
        # f.write("plot \"TAULINE vs EVLINE 2\"\n")
        # f.write("  xvar evline 2\n")
        # f.write("  yvar tauline 2 Nnode \n")
        # f.write("\n")
        
        # f.write("plot \"TAULINE vs EVLINE 3\"\n")
        # # f.write("averaged\n")
        # f.write("  xvar evline 3\n")
        # f.write("  yvar tauline 3 Nnode \n")
        # f.write("\n")
        
        # f.write("plot \"JSP vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jsp 0 1 0\n")
        # f.write("\n")
        
        f.write("plot \"SPECTRUM vs ENERGY\"\n")
        f.write("  xvar sp_energy\n")
        f.write("  yvar jsp 0 1\n")
        f.write("  yvar jsp 0 51\n")

        # f.write("plot \"JSP vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jsp 0 1 0\n")
        # f.write("\n")        

        f.write("plot \"TAU vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"TAU vs ENERGY\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap 0 Nnode\n")
        f.write("\n")
        
        f.write("plot \"NE vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  ne 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"ERAD vs TIME\"\n")
        f.write(" xvar time\n")
        f.write(" yvar erad 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"EIGENVAL vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  eigenv_r 1 1:51 0 -1:\n")
        f.write("\n")

        f.write("plot \"TEV, TIV vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  tev 0 1\n")
        f.write("  yvar  tiv 0 1 \n")
        f.write("\n")
        
        f.write("plot \"YISO vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar yisofrac 1 1 0:10\n")
        f.write("\n")

        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar ntry\n")
        f.write("\n")
        
        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar dtime\n")
        f.write("\n")
        
        f.write("#ifdef DISPLAY\n")
        f.write("display 1\n")
        f.write("display 2\n")
        f.write("display 3\n")
        f.write("display 4\n")
        f.write("display 5\n")
        f.write("#endif\n")
        f.close()  
#%% Lines 
base_dir    ='/home/jeff/Research/Export_Control/Cretin/Neon_Cell/lines/v12/'
com         = 'c ------------------------------------------------------------\n'
pressures   = [7.5,15,30]# using 2.5e17,5e17,1e18 as densities
densities   = [2.5e17,5e17,1e18]
gens        = ['gen1c']
g2xfile     = {'gen1c':'myl_c_xfile_interp.txt','gen1f':'myl_only_all_F.txt',  'gen2c':'SiN_only_all_C.txt','gen2f':'SiN_only_all_F.txt'}
term        = '/home/jeff/Research/Export_Control/Cretin/xfiledir/term10.dat'
cretin_ex   = '/home/jeff/Research/Export_Control/Cretin/cretin'
onezone     = False
rad         = False
use_escape  = True

# gens = ['gen1f']
# pressures = [15]
# densities = [5e17]
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)
for gen in gens:
    directory = os.path.join(base_dir,gen)
    xfile = g2xfile[gen]
    if not os.path.isdir(directory):
        os.mkdir(directory)
    for p in range(len(pressures)):
        dens = densities[p]
        runpath = os.path.join(directory,str(pressures[p])+"torr")
        if not os.path.isdir(runpath):
            os.mkdir(runpath)
        os.chdir(runpath)
        fname = 'ne_cell.gen'
        shutil.copy(term, runpath)
        shutil.copy(cretin_ex, runpath)
        f = open(fname,"w")
        # f.write("\n")
        f.write("c        *** Neon Slab No Windows ***\n")
        f.write(f"alias N_neon	{dens}\t\t\t! neon number density\n")
        if onezone:
            f.write("alias Nnode	1\n")
        else:
            f.write("alias Nnode	51\n")
        f.write("\n")
        f.write("alias T0	0.024\n") # 276K in eV
        f.write("alias T1	0.1\n")
        f.write("\n")
        f.write("alias SIZE	1.35\n")
        f.write("alias DTMIN 	1e-15\n")
        f.write("alias DTMAX 	1e-7\n")
        # f.write("\n")
        # f.write("alias N0        1\n")
        # f.write("alias N1        Nnode+/2\n")
        # f.write("alias N2        Nnode\n")
        # f.write("alias DN        25\n")
        # f.write("\n")
        # f.write("alias R0        0.\n")
        # f.write("alias R1        R0 + SIZE/2.\n")
        # f.write("alias R2        R0 + SIZE\n")
        # f.write("alias DR        0.01\n")
        # f.write("\n")
        # f.write("alias S31       1\n")
        # f.write("alias S44       10\n")
        f.write("\n")
        f.write(com)
        f.write("c   Materials\n")
        f.write(com)
        f.write("atoms hydrogenic ne\n")
        f.write("  modeltype dca \n")
        f.write("\n")
        f.write("region  1 Nnode  T0		! Neon Slab\n")
        f.write("  element  1  N_neon\n")
        f.write("\n")
        f.write("background 0. 1.e15 \t\t! add background electrons to fix sw 31 = 4\n")
        f.write(com)
        f.write("c   Geometry\n")
        f.write(com)
        f.write("geometry slab\n")
        
        if not onezone:
            f.write("rlin 1 Nnode 0 SIZE 		! Equally spaced zones\n")
        f.write(com)
        f.write("c   Radiation\n")
        f.write(com)
        f.write("angles 3\n")
        f.write("\n")
        f.write("alias emin  1.e-1\n")
        f.write("alias emax  4.e3\n")
        
        f.write("ebins 1e3 emin emax\n")
        f.write("spectrum 500 .01 900 1\n")
        f.write("spectrum 500 900 1300 1\n")
        
        f.write("spectrum 500 1300 4000 1\n")
        f.write("\n")
        f.write(f"xfile ix=12 {xfile}   ! define xfile 12\n")
        f.write("\n")
        f.write("source jbndry 12 pbins 1 1.          ! define jbndry from xfile 12\n")
        f.write("xfilebc 12 -1 0 1. 1                 ! use jbndry from xfile 12\n")
        f.write(com)
        f.write("\n")
        f.write("linedefault crd \n")
        f.write(com)
        f.write("c H-like Lyman series to n=5\n")
        f.write(com)
        lineid = 1
        for i in range(2,6):
            f.write(f"line {lineid} 1 1 1 1 {i}     \t! H-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c He-like Lyman series to n=10\n")
        f.write(com)
        for i in range(2,11):
            f.write(f"line {lineid} 1 2 1 2  {i}     \t! He-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c Li-like Lyman series to n=10\n")
        f.write(com)
        for i in range(2,11):
            f.write(f"line {lineid} 1 3 1 3 {i}      \t! Li-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c Be-like Lyman series to n=5\n")
        f.write(com)
        for i in range(2,6):
            f.write(f"line {lineid} 1 3 1 3 {i}      \t! Be-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c B-like Lyman series to n=5\n")
        f.write(com)
        for i in range(2,6):
            f.write(f"line {lineid} 1 4 1 4 {i}      \t! B-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c C-like Lyman series to n=5\n")
        f.write(com)
        for i in range(2,6):
            f.write(f"line {lineid} 1 5 1 5 {i}      \t! C-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c N-like Lyman series to n=5\n")
        f.write(com)
        for i in range(2,6):
            f.write(f"line {lineid} 1 6 1 6 {i}      \t! N-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        f.write("c O-like Lyman series to n=5\n")
        f.write(com)
        for i in range(2,6):
            f.write(f"line {lineid} 1 6 1 6 {i}      \t! O-like \n")
            f.write("    lbins 50 5. 1.\n")
            lineid = lineid + 1
        f.write(com)
        # f.write("c   Controls\n")
        f.write(com)
        f.write("tstart 0.\n")
        f.write("tquit  1.08e-07\n")
        f.write("restart\n")
        f.write("nltedump spectrum transmission\n")
        f.write(com)
        f.write("c   Switches and Parameters\n")
        f.write(com)
        f.write("switch 2 0			\t! hydrodynamics 0=off\n")
        f.write("switch 11  1                    \t! make .plt file\n")
        f.write("switch 20  1                  	\t! calculate NLTE populations using rate matrices\n")
        # f.write("switch 21  1                    \t! account for transition energies\n")
        f.write("switch 22  1                    \t! recalculate destruction rates on first iteration\n")
        f.write("switch 25  1                    \t! time dependent kinetics = 1, SS = 0\n")
        f.write("switch 28  2                    \t! 2: initialize in steady-state kinetics w/o radiation transfer\n")
        f.write("switch 29  2                    \t! use variable timesteps\n")
        f.write("switch 30 10                    \t! dump every n timesteps\n")
        f.write("switch 31 1                    \t\t! do temperature calculation TD = 1, SS = -1,rad hydro = 4\n")
        f.write("switch 33 1			\t! use excape factors for all photoexcitations if >0\n")
        f.write("switch 36 1			     \t! radiation transport on = 1, off = 0\n")
        f.write("switch 37 1                     \t! do line transfer = 1\n")
        f.write("switch 44 S44                   \t! max iterations per timestep\n")
        # f.write("switch 49 1                   \t! Electron thermal conduction on = 1\n")
        # f.write("switch 51 0                   \t! Electron thermal conduction coeff Spitzer-Harm")
        f.write("switch 52 -1                    \t! do Stark broadening for everything\n")
        f.write("switch 55 1                     \t! do continuum lowering\n")
        f.write("\n")
        f.write("param  40 DTMAX 		\t! time between edits\n")
        f.write("param  41 DTMIN 		\t! initial timestep\n")
        f.write("param  44 DTMIN 		\t! minimum timestep\n")
        f.write("param  45 DTMAX 		\t! maximum timestep\n")
        f.write("param  42  0.1                 	\t! max change in Zbar\n")
        f.write("param  46  0.05                 \t! maximum frac. change in temperature\n")
        f.write("param  61  1.e-8               	\t! iso-sequence population threshold\n")
        f.write("param  102  T1              	\t! min Te\n")
        f.write(com)
        f.write("c   Edits\n")
        f.write(com)

        f.write("\n")
        
        f.write("plot \"TEV, TIV, TRADV vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tev\n")
        f.write("  yvar  tiv\n")
        f.write("  yvar  tradv\n")
        f.write("\n")
        
        f.write("plot \"NE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  ne\n")
        f.write("\n")
        
        f.write("plot \"dNEdt vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  dnedt \n")
        f.write("\n")

        
        f.write("plot \"ERAD vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  erad\n")
        f.write("\n")
        
        f.write("plot \"TAUEE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauee\n")
        f.write("\n")
        
        f.write("plot \"TAUEI vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauei\n")
        f.write("\n")
        
        f.write("plot \"TAUII vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tauii\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("xvar r\n")
        f.write("yvar heatt\n")
        f.write("\n")        
        

        f.write("plot \"YISO vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 0:5\n")
        f.write("\n")


        f.write("plot \"YISO vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 6:10\n")
        f.write("\n")
        

        f.write("plot \"JNU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  jnu 0 1:-1:DN\n")
        f.write("\n")
        
        f.write("plot \"TAU vs E\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"JBAR vs TIME\"\n")
        # f.write("averaged\n")
        f.write("  xvar time\n")
        f.write("  yvar jbarcgs 1 1\n")
        f.write("  yvar jbarcgs 2 1\n")
        f.write("  yvar jbarcgs 3 1\n")
        f.write("\n")
        
        f.write("plot \"JSAT vs TIME\"\n")
        # f.write("averaged\n")
        f.write("  xvar time\n")
        f.write("  yvar jsat 1 1\n")
        f.write("  yvar jsat 2 1\n")
        f.write("  yvar jsat 3 1\n")
        f.write("\n")
        
        # f.write("plot \"TAULINE vs EVLINE 1\"\n")
        # f.write("  xvar evline 1\n")
        # f.write("  yvar tauline 1 Nnode \n")
        # f.write("\n")
        
        # f.write("plot \"TAULINE vs EVLINE 2\"\n")
        # f.write("  xvar evline 2\n")
        # f.write("  yvar tauline 2 Nnode \n")
        # f.write("\n")
        
        # f.write("plot \"TAULINE vs EVLINE 3\"\n")
        # # f.write("averaged\n")
        # f.write("  xvar evline 3\n")
        # f.write("  yvar tauline 3 Nnode \n")
        # f.write("\n")
        
        # f.write("plot \"JSP vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jsp 0 1 0\n")
        # f.write("\n")
        
        f.write("plot \"SPECTRUM vs ENERGY\"\n")
        f.write("  xvar sp_energy\n")
        f.write("  yvar jsp 0 1\n")
        f.write("  yvar jsp 0 51\n")

        # f.write("plot \"JSP vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jsp 0 1 0\n")
        # f.write("\n")        

        f.write("plot \"TAU vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"TAU vs ENERGY\"\n")
        f.write("  xvar  energy\n")
        f.write("  yvar  taukap 0 Nnode\n")
        f.write("\n")
        
        f.write("plot \"NE vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  ne 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"ERAD vs TIME\"\n")
        f.write(" xvar time\n")
        f.write(" yvar erad 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"EIGENVAL vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  eigenv_r 1 1:51 0 -1:\n")
        f.write("\n")

        f.write("plot \"TEV, TIV vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  tev 0 1\n")
        f.write("  yvar  tiv 0 1 \n")
        f.write("\n")
        
        f.write("plot \"YISO vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar yisofrac 1 1 0:10\n")
        f.write("\n")

        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar ntry\n")
        f.write("\n")
        
        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar dtime\n")
        f.write("\n")
        
        f.write("#ifdef DISPLAY\n")
        f.write("display 1\n")
        f.write("display 2\n")
        f.write("display 3\n")
        f.write("display 4\n")
        f.write("display 5\n")
        f.write("#endif\n")
        f.close() 
#%% Thermal Conduction test
base_dir    ='/home/jeff/Research/Export_Control/Cretin/Neon_Cell/thermal_timescale/default/v2/'
com         = 'c ------------------------------------------------------------\n'
pressures   = [7.5,15,30]# using 2.5e17,5e17,1e18 as densities
densities   = [2.5e17,5e17,1e18]
gens        = ['gen1c']
g2xfile     = {'gen1c':'myl_c_xfile_interp.txt','gen1f':'myl_only_all_F.txt',  'gen2c':'SiN_only_all_C.txt','gen2f':'SiN_only_all_F.txt'}
term        = '/home/jeff/Research/Export_Control/Cretin/xfiledir/term10.dat'
cretin_ex   = '/home/jeff/Research/Export_Control/Cretin/cretin'
onezone     = False
rad         = False
use_escape  = True

# gens = ['gen1f']
# pressures = [15]
# densities = [5e17]
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)
for gen in gens:
    directory = os.path.join(base_dir,gen)
    xfile = g2xfile[gen]
    if not os.path.isdir(directory):
        os.mkdir(directory)
    for p in range(len(pressures)):
        dens = densities[p]
        runpath = os.path.join(directory,str(pressures[p])+"torr")
        if not os.path.isdir(runpath):
            os.mkdir(runpath)
        os.chdir(runpath)
        fname = 'ne_cell.gen'
        shutil.copy(term, runpath)
        shutil.copy(cretin_ex, runpath)
        f = open(fname,"w")
        # f.write("\n")
        f.write("c        *** Neon Slab No Windows ***\n")
        f.write(f"alias N_neon	{dens}\t\t\t! neon number density\n")
        if onezone:
            f.write("alias Nnode	1\n")
        else:
            f.write("alias Nnode	51\n")
        f.write("\n")
        f.write("alias T0	0.024\n") # 276K in eV
        f.write("alias T1	0.1\n")
        f.write("\n")
        f.write("alias SIZE	1.35\n")
        f.write("alias DTMIN 	1e-12\n")
        f.write("alias DTMAX 	1e-7\n")
        f.write("\n")
        # f.write("alias N0        1\n")
        # f.write("alias N1        Nnode+/2\n")
        # f.write("alias N2        Nnode\n")
        # f.write("alias DN        25\n")
        # f.write("\n")
        # f.write("alias R0        0.\n")
        # f.write("alias R1        R0 + SIZE/2.\n")
        # f.write("alias R2        R0 + SIZE\n")
        # f.write("alias DR        0.01\n")
        f.write("\n")
        f.write("alias S31       1\n")
        f.write("alias S44       10\n")
        f.write("\n")
        f.write(com)
        f.write("c   Materials\n")
        f.write(com)
        f.write("atoms hydrogenic ne\n")
        f.write("  modeltype dca term\n")
        f.write("\n")
        f.write("region  1 Nnode  T0		! Neon Slab\n")
        f.write("  element  1  N_neon\n")
        f.write(com)
        f.write("c   Geometry\n")
        f.write(com)
        f.write("geometry slab\n")
        
        if not onezone:
            f.write("rlin 1 Nnode 0 SIZE 		! Equally spaced zones\n")
        f.write(com)
        f.write("c   Radiation\n")
        f.write(com)
        f.write("angles 3\n")
        f.write("\n")
        f.write("alias emin  1.e-1\n")
        f.write("alias emax  4.e3\n")
        f.write("\n")
        f.write(f"xfile ix=12 {xfile}   ! define xfile 12\n")
        f.write("\n")
        f.write("source jbndry 12 pbins 1 1.          ! define jbndry from xfile 12\n")
        f.write("xfilebc 12 -1 0 1. 1                  ! use jbndry from xfile 12\n")
        # f.write("background 0. 1.e15 ! add background electrons to fix sw 31 = 4\n")
        f.write("ebins 1e3 emin emax\n")
        
        f.write(com)
        f.write("c   Controls\n")
        f.write(com)
        f.write("tstart 0.\n")
        f.write("tquit  1.08e-07\n")
        f.write("restart\n")
        # f.write("nltedump spectrum transmission\n")
        f.write(com)
        f.write("c   Switches and Parameters\n")
        f.write(com)

        if rad:
            f.write("switch 100 1\n")
        else:
            f.write("switch 31  1                    ! do temperature calculation TD = 1, SS = -1,rad hydro = 4\n")
            f.write("switch 36 1			! radiation transport on = 1, off = 0\n")
        if use_escape:
            f.write("switch 33 1			! use excape factors for all photoexcitations if >0\n")
        if not onezone:
            f.write("switch 2 0			! hydrodynamics 0 = off 1 = on\n")
        f.write("switch 11  1                    ! make .plt file\n")
        f.write("switch 20  1                  	! calculate NLTE populations using rate matrices\n")
        # f.write("switch 21  1                    ! account for transition energies\n")
        f.write("switch 22  1                    ! recalculate destruction rates on first iteration\n")
        f.write("switch 25  1                    ! time dependent kinetics = 1, SS = 0\n")
        f.write("switch 28  2                    ! 2: initialize in steady-state kinetics w/o radiation transfer\n")
        f.write("switch 29  2                    ! use variable timesteps\n")
        f.write("switch 30 20                    ! dump every n timesteps\n")
        f.write("switch 37 1                     ! do line transfer\n")
        f.write("switch 44 S44                   ! max iterations per timestep\n")
        # f.write("switch 49 1                   ! Electron thermal conduction on = 1\n")
        # f.write("switch 51 0                   ! Electron thermal conduction coeff Spitzer-Harm")
        f.write("\n")
        f.write("param  40 DTMAX 		! time between edits\n")
        f.write("param  41 DTMIN 		! initial timestep\n")
        f.write("param  44 DTMIN 		! minimum timestep\n")
        f.write("param  45 DTMAX 		! maximum timestep\n")
        f.write("param  42  0.1                 	! max change in Zbar\n")
        f.write("param  46  0.05                 ! maximum frac. change in temperature\n")
        f.write("param  61  1.e-8               	! iso-sequence population threshold\n")
        f.write("param  102  T1              	! min Te\n")
        f.write(com)
        f.write("c   Edits\n")
        f.write(com)

        f.write("\n")
        
        f.write("plot \"TEV, TIV, TRADV vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  tev\n")
        f.write("  yvar  tiv\n")
        f.write("  yvar  tradv\n")
        f.write("\n")
        
        f.write("plot \"NE vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  ne\n")
        f.write("\n")
        
        f.write("plot \"dNEdt vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  dnedt \n")
        f.write("\n")

        
        f.write("plot \"ERAD vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  erad\n")
        f.write("\n")
        
        # f.write("plot \"TAUEE vs R\"\n")
        # f.write("  xvar  r\n")
        # f.write("  yvar  tauee\n")
        # f.write("\n")
        
        # f.write("plot \"TAUEI vs R\"\n")
        # f.write("  xvar  r\n")
        # f.write("  yvar  tauei\n")
        # f.write("\n")
        
        # f.write("plot \"TAUII vs R\"\n")
        # f.write("  xvar  r\n")
        # f.write("  yvar  tauii\n")
        # f.write("\n")
        
        f.write("plot \"TOTAL HEAT vs R\"\n")
        f.write("xvar r\n")
        f.write("yvar heatt\n")
        f.write("\n")        
        

        f.write("plot \"YISO vs R 1\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 0:5\n")
        f.write("\n")


        f.write("plot \"YISO vs R 2\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar yisofrac 1 0 6:10\n")
        f.write("\n")
        

        # f.write("plot \"JNU vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jnu 0 1:-1:DN\n")
        # f.write("\n")
        
        # f.write("plot \"TAU vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  taukap\n")
        # f.write("\n")
        
        # f.write("plot \"JSP vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jsp 0 1 0\n")
        # f.write("\n")
        
        # f.write("plot \"SPECTRUM vs ENERGY\"\n")
        # f.write("  xvar sp_energy\n")
        # f.write("  yvar jsp 0 1\n")
        # f.write("  yvar jsp 0 51\n")

        # f.write("plot \"JSP vs E\"\n")
        # f.write("  xvar  energy\n")
        # f.write("  yvar  jsp 0 1 0\n")
        # f.write("\n")        

        f.write("plot \"TAU vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  taukap\n")
        f.write("\n")
        
        f.write("plot \"NET HEAT vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  heatjt\n")
        f.write("\n")
        
        f.write("plot \"Specific HEAT vs R\"\n")
        f.write("  xvar  r\n")
        f.write("  yvar  cv\n")
        f.write("\n")
        
        f.write("plot \"NE vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  ne 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"ERAD vs TIME\"\n")
        f.write(" xvar time\n")
        f.write(" yvar erad 0 1:-1:DN\n")
        f.write("\n")


        
        f.write("plot \"EIGENVAL vs TIME\"\n")
        f.write("  xvar  time\n")
        f.write("  yvar  eigenv_r 1 1:51 0 -1:\n")
        f.write("\n")

        f.write("plot \"TEV, TIV vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar  tev 0 1\n")
        f.write("  yvar  tiv 0 1 \n")
        f.write("\n")
        
        f.write("plot \"YISO vs TIME\"\n")
        f.write("averaged\n")
        f.write("  xvar  time\n")
        f.write("  yvar yisofrac 1 1 0:10\n")
        f.write("\n")

        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar ntry\n")
        f.write("\n")
        
        f.write("plot\n")
        f.write(" xvar time\n")
        f.write(" yvar dtime\n")
        f.write("\n")
        
        f.write("#ifdef DISPLAY\n")
        f.write("display 1\n")
        f.write("display 2\n")
        f.write("display 3\n")
        f.write("display 4\n")
        f.write("display 5\n")
        f.write("#endif\n")
        f.close()