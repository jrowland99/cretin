#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 12:45:03 2023

@author: jeff
"""
import pandas as pd
import os
from labellines import labelLines
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
plt.rcParams["figure.dpi"] = 350
plt.rcParams["lines.linewidth"] = 1
plt.rcParams.update({"font.size": 10})
plt.rcParams.update({"legend.fontsize": 10})
plt.rcParams["figure.figsize"] = (9, 4)
ext = ".png"
matplotlib.rcParams.update(
    {"text.usetex": False, "font.family": "stixgeneral", "mathtext.fontset": "stix",}
)

class CretinData:
    def __init__(self, **kwargs):
        self.name           = kwargs.get('name', None)
        self.path           = kwargs.get('path', None) # Contains the root run directory
        parent_folder       = os.path.basename(os.path.dirname(self.path))
        grandparent_folder  = os.path.basename(os.path.dirname(os.path.dirname(self.path)))
        self.title_str      = f"{grandparent_folder}_{parent_folder}"
        self.gen            = kwargs.get('gen', None)
        self.pressure       = kwargs.get('pressure', None)
        self.ext            = '.png'
        # time data holds spatially integrated, 1-D time-dependent data
        self.time_data      = {} # header times are rounded to 3 decimals
        self.r              = [] # common space axis
        # Data holds space & time or energy & time data, 2-D data in general
        self.data           = {} # Store data that is only in dict
        self.ne_ions        = ["Ne", "Ne+", "Ne+2", "Ne+3", "Ne+4", "Ne+5", "Ne+6", "Ne+7", "Ne+8", "Ne+9", "Ne+10"]
        self.yiso_labels    = { # number of bound electrons is the last number in the id string
            "Ne":'yisofrac<r> 1,1,10',
            "Ne+":'yisofrac<r> 1,1,9',
            "Ne+2":'yisofrac<r> 1,1,8',
            "Ne+3":'yisofrac<r> 1,1,7',
            "Ne+4":'yisofrac<r> 1,1,6',
            "Ne+5":'yisofrac<r> 1,1,5',
            "Ne+6":'yisofrac<r> 1,1,4',
            "Ne+7":'yisofrac<r> 1,1,3',
            "Ne+8":'yisofrac<r> 1,1,2',
            "Ne+9":'yisofrac<r> 1,1,1',
            "Ne+10":'yisofrac<r> 1,1',
            }

    def read_data(self):
        with open(os.path.join(self.path,self.name+'.plt'), 'r') as f:

            chunks = f.read().split('#')
            chunks = chunks[1:]
            for i,chunk in enumerate(chunks):
                # print(f"starting chunk {i}")
                lines = chunk.split('\n')
                info1 = lines[1]  # Lines[0] is # header line, these are column label names
                info2 = lines[2]  # this is the index info line
                lines = [line for line in lines if not line.startswith('$')]
                header = lines[0]
                lines = lines[1:]
                if ':' in header.strip(): # This is a section with a timestamp in the header
                    tstep = float(header.split(':')[0].split('=')[-1])
                    header = header.split(':')[-1].strip()
                    cols = info1.split('$')[-1].split()
                    if info2.split() != '$': #if this line is not empty

                        cols2 = info2.split('$')[-1].split()
                        # print(f'running line 48 {cols2}')
                    if cols2:#if we did the previous check
                        cols3 = []
                        for i,col in enumerate(cols2):
                            cols3.append(cols[i+1]+' '+col)# the i+1 will skip the xvar column and append indices to the yvars
                            # print(f'running line 53 {cols3}')
                        if cols3:
                            dummylist = []
                            dummylist.append(cols[0])
                            dummylist.extend(cols3)
                            cols = dummylist
                    temp_dict = {key:[] for key in cols}
                    # print("just initialized ",temp_dict)
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            for i,col in enumerate(cols):
                                temp_dict[col].append(items[i]) # Populate temporary dictionary with the values in each row
                    # print(temp_dict)
                    dummy_df = pd.DataFrame(temp_dict)
                    if header not in self.data: # Insert a new timestep df into self.data
                        self.data[header] = {}
                        self.data[header][tstep] = dummy_df
                    elif header in self.data:
                        self.data[header][tstep] = dummy_df
                    else:
                        print("This shouldnt happen but check your data entry")
                # Check for different data types
                elif header.strip() == 'TEV, TIV vs TIME':
                    # print('exec TEV, TIV vs TIME')
                    self.time_data['TEV'], self.time_data['TIV'],self.time_data['TIME'] = [], [],[]
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            time = float(items[0]) * 1e9
                            tev = float(items[1])
                            tiv = float(items[2])
                            self.time_data['TIME'].append(time)
                            self.time_data['TEV'].append(tev)
                            self.time_data['TIV'].append(tiv)

                elif header.strip() == 'NE vs TIME':
                    # print('exec NE vs TIME')
                    self.time_data['NE'] = []
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            time = float(items[0]) * 1e9
                            ne = float(items[1])
                            self.time_data['NE'].append(ne)
                elif header.strip() == 'YISO vs TIME':
                    # print('exec YISO vs TIME')
                    cols = info1.split('$')[-1].split()
                    if info2.split() != '$': #if this line is not empty
                        cols2 = info2.split('$')[-1].split()
                        # print(f'running line 103 {cols2}')
                    if cols2:#if we did the previous check
                        cols3 = []
                        for i,col in enumerate(cols2):
                            cols3.append(cols[i+1]+' '+col)# the i+1 will skip the xvar column and append indices to the yvars
                            # print(f'running line 108 {cols3}')
                        if cols3:
                            dummylist = []
                            dummylist.append(cols[0])
                            dummylist.extend(cols3)
                            cols = dummylist
                    temp_dict = {key:[] for key in cols}
                    # print("just initialized ",temp_dict)
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            for i,col in enumerate(cols):
                                temp_dict[col].append(items[i]) # Populate temporary dictionary with the values in each row
                    # print("just populated yiso",temp_dict)
                    if 'TIME' not in self.time_data:
                        self.time_data['TIME'] = []
                    for key, val in temp_dict.items():
                        val = [float(v) for v in val]
                        self.time_data[key]=np.array(val)
                # self.time_data = {key: float(value) for key, value in self.time_data.items()}# make sure floats are floats
                elif header.strip() == 'GAMMATOT vs TIME':
                    print('exec GAMMATOT vs TIME')
                    cols = info1.split('$')[-1].split()
                    if info2.split() != '$': #if this line is not empty
                        cols2 = info2.split('$')[-1].split()
                        # print(f'running line 103 {cols2}')
                    if cols2:#if we did the previous check
                        cols3 = []
                        for i,col in enumerate(cols2):
                            cols3.append(cols[i+1]+' '+col)# the i+1 will skip the xvar column and append indices to the yvars
                            # print(f'running line 108 {cols3}')
                        if cols3:
                            dummylist = []
                            dummylist.append(cols[0])
                            dummylist.extend(cols3)
                            cols = dummylist
                    temp_dict = {key:[] for key in cols}
                    # print("just initialized ",temp_dict)
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            for i,col in enumerate(cols):
                                temp_dict[col].append(items[i]) # Populate temporary dictionary with the values in each row
                    # print("just populated yiso",temp_dict)
                    if 'TIME' not in self.time_data:
                        self.time_data['TIME'] = []
                    for key, val in temp_dict.items():
                        val = [float(v) for v in val]
                        self.time_data[key]=np.array(val)
                
                
                elif header.strip() == 'EIGENVAL vs TIME':
                    print('exec EIGENVAL vs TIME')
                    cols = info1.split('$')[-1].split()
                    if info2.split() != '$': #if this line is not empty
                        cols2 = info2.split('$')[-1].split()
                        # print(f'running line 103 {cols2}')
                    if cols2:#if we did the previous check
                        cols3 = []
                        for i,col in enumerate(cols2):
                            cols3.append(cols[i+1]+' '+col)# the i+1 will skip the xvar column and append indices to the yvars
                            # print(f'running line 108 {cols3}')
                        if cols3:
                            dummylist = []
                            dummylist.append(cols[0])
                            dummylist.extend(cols3)
                            cols = dummylist
                    temp_dict = {key:[] for key in cols}
                    # print("just initialized ",temp_dict)
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            for i,col in enumerate(cols):
                                temp_dict[col].append(items[i]) # Populate temporary dictionary with the values in each row
                    # print("just populated yiso",temp_dict)
                    if 'TIME' not in self.time_data:
                        self.time_data['TIME'] = []
                    for key, val in temp_dict.items():
                        val = [float(v) for v in val]
                        self.time_data[key]=np.array(val)
                        
                elif header.strip() == 'EIGENVAL vs TIME desc':
                    print('exec EIGENVAL vs TIME desc')
                    cols = info1.split('$')[-1].split()
                    if info2.split() != '$': #if this line is not empty
                        cols2 = info2.split('$')[-1].split()
                        # print(f'running line 103 {cols2}')
                    if cols2:#if we did the previous check
                        cols3 = []
                        for i,col in enumerate(cols2):
                            cols3.append(cols[i+1]+' '+col)# the i+1 will skip the xvar column and append indices to the yvars
                            # print(f'running line 108 {cols3}')
                        if cols3:
                            dummylist = []
                            dummylist.append(cols[0])
                            dummylist.extend(cols3)
                            cols = dummylist
                    temp_dict = {key:[] for key in cols}
                    # print("just initialized ",temp_dict)
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            for i,col in enumerate(cols):
                                temp_dict[col].append(items[i]) # Populate temporary dictionary with the values in each row
                    # print("just populated yiso",temp_dict)
                    if 'TIME' not in self.time_data:
                        self.time_data['TIME'] = []
                    for key, val in temp_dict.items():
                        val = [float(v) for v in val]
                        self.time_data[key]=np.array(val)
                        
                elif header.strip() == 'JBAR vs TIME':
                    print('exec JBAR vs TIME')
                    cols = info1.split('$')[-1].split()
                    if info2.split() != '$': #if this line is not empty
                        cols2 = info2.split('$')[-1].split()
                        # print(f'running line 103 {cols2}')
                    if cols2:#if we did the previous check
                        cols3 = []
                        for i,col in enumerate(cols2):
                            cols3.append(cols[i+1]+' '+col)# the i+1 will skip the xvar column and append indices to the yvars
                            # print(f'running line 108 {cols3}')
                        if cols3:
                            dummylist = []
                            dummylist.append(cols[0])
                            dummylist.extend(cols3)
                            cols = dummylist
                    temp_dict = {key:[] for key in cols}
                    # print("just initialized ",temp_dict)
                    for line in lines:
                        if len(line) > 0:
                            items = line.split()
                            for i,col in enumerate(cols):
                                temp_dict[col].append(items[i]) # Populate temporary dictionary with the values in each row
                    # print("just populated yiso",temp_dict)
                    if 'TIME' not in self.time_data:
                        self.time_data['TIME'] = []
                    for key, val in temp_dict.items():
                        val = [float(v) for v in val]
                        self.time_data[key]=np.array(val)
                
    def make_animation(self, **kwargs):
        from matplotlib.animation import FuncAnimation
        save_dir = kwargs.get('save_dir', self.path) # if no save dir defined do it in run dir
        edit = kwargs.get('edit', None)
        xvar = kwargs.get('xvar', None) # string variable for x-axis, like 'r' or 'E'
        yvar = kwargs.get('yvar', None) # the parameter to plot like Te or <Z> etc
        if edit == None or xvar == None or yvar == None:
            print('Must define all variables for the plot')
            return
        print(f"Animating {xvar} {yvar}")
        data_obj    = self.data[edit]
        frames      = np.asarray(list(data_obj.keys())) # These are points in time
        x_axis      = list(data_obj[frames[0]][xvar])
        x_axis      = np.asarray([float(x) for x in x_axis])
        y_min       = min(np.array(list(data_obj[frames[0]][yvar])).astype(float))
        y_max = 0
        for time in frames:
            new_max = max(np.array(list(data_obj[time][yvar])).astype(float))
            new_min = min(np.array(list(data_obj[time][yvar])).astype(float))
            if new_max > y_max:
                y_max = new_max
            if new_min < y_min:
                y_min = new_min

        print(y_min)
        print(y_max)

        fig, ax = plt.subplots(figsize=(5, 5))
        lines = []  # List to store line objects
        ax.set_ylim([y_min, y_max+0.05*y_max])
        # ax.set_ylim([y_min, 1e7])
        ax.set_xlim([min(x_axis),max(x_axis)])
        if xvar == 'energy':
            ax.set_xlim([0,2000])
        ax.set_xlabel(xvar)
        ax.set_ylabel(yvar)
        def init():
            for line in lines:
                line.set_data([], [])
            ax.set_title('')
            return lines
        def animate(time):
            y = np.array(data_obj[time][yvar].values).astype(float)

            # Add new line to the list and plot it
            line, = ax.plot(x_axis, y, label=str(time*1e9) + 'ns',color='red')
            lines.append(line)

            # Remove old lines if there are more than 10
            if len(lines) > 3:
                old_line = lines.pop(0)
                old_line.remove()

            # Decrease opacity for previous lines
            for i, l in enumerate(lines[:-1]):
                l.set_alpha(0.1 + (0.9 / len(lines)) * (i + 1))
                l.set_linestyle('--')
                l.set_color('slategrey')

            ax.set_title(f"Cretin {yvar} vs {xvar} at t = {time*1e9:.2f}"+'ns'+f' {self.gen} {self.pressure} torr',loc='left')
            return lines


        time_ax     = np.asarray(list(self.data[edit].keys()))
        ani = FuncAnimation(fig, animate, init_func=init, frames=time_ax, blit=True)
        plt.grid(linestyle=':', linewidth=0.5)
        # plt.tight_layout()
        if edit == 'DERIVATIVES':
            ani.save(save_dir+f"{yvar}_{xvar}_time"+'.mp4', writer='ffmpeg')
        else:
            ani.save(save_dir+self.title_str+f"_{yvar}_{xvar}_time"+'.mp4', writer='ffmpeg')
        return
    
    def coplot_movie(self, **kwargs):
        from matplotlib.animation import FuncAnimation
        save_dir = kwargs.get('save_dir', self.path)
        edit1 = kwargs.get('edit1', None)  # Edit variable for the first y-variable
        edit2 = kwargs.get('edit2', None)  # Edit variable for the second y-variable
        xvar = kwargs.get('xvar', None)
        yvar1 = kwargs.get('yvar1', None)
        yvar2 = kwargs.get('yvar2', None)
        
        if edit1 is None or edit2 is None or xvar is None or yvar1 is None or yvar2 is None:
            print('Must define all variables for the plot')
            return
        
        data_obj1 = self.data[edit1]
        data_obj2 = self.data[edit2]
        
        frames = np.asarray(list(data_obj1.keys()))
        x_axis = list(data_obj1[frames[0]][xvar])
        x_axis = np.asarray([float(x) for x in x_axis])
        
        y1_min = min(np.array(list(data_obj1[frames[0]][yvar1])).astype(float))
        y1_max = 0
        y2_min = min(np.array(list(data_obj2[frames[0]][yvar2])).astype(float))
        y2_max = 0
        
        for time in frames:
            new_y1_max = max(np.array(list(data_obj1[time][yvar1])).astype(float))
            new_y2_max = max(np.array(list(data_obj2[time][yvar2])).astype(float))
            
            if new_y1_max > y1_max:
                y1_max = new_y1_max
            if new_y2_max > y2_max:
                y2_max = new_y2_max
        
        print(y1_min)
        print(y1_max)
        print(y2_min)
        print(y2_max)
        
        fig, ax1 = plt.subplots(figsize=(5, 5))
        ax2 = ax1.twinx()
        
        lines1 = []
        lines2 = []
        
        ax1.set_ylim([y1_min, y1_max + 0.05 * y1_max])
        ax2.set_ylim([y2_min, y2_max + 0.05 * y2_max])
        
        ax1.set_xlim([min(x_axis), max(x_axis)])
        if xvar == 'energy':
            ax1.set_xlim([0, 2000])
        
        ax1.set_xlabel(xvar)
        ax1.set_ylabel(yvar1, color='tab:blue')
        ax2.set_ylabel(yvar2, color='tab:red')
        
        def init():
            for line in lines1 + lines2:
                line.set_data([], [])
            ax1.set_title('')
            return lines1 + lines2
        
        def animate(time):
            y1 = np.array(data_obj1[time][yvar1].values).astype(float)
            y2 = np.array(data_obj2[time][yvar2].values).astype(float)
            
            for line in lines1:
                line.remove()
            for line in lines2:
                line.remove()
            
            lines1.clear()
            lines2.clear()
            
            line1, = ax1.plot(x_axis, y1, label=str(time * 1e9) + 'ns', color='blue')
            line2, = ax2.plot(x_axis, y2, label=str(time * 1e9) + 'ns', color='red')
            
            lines1.append(line1)
            lines2.append(line2)
            
            ax1.set_title(f"Cretin {yvar1} vs {xvar} at t = {time * 1e9:.2f}" + 'ns' + f' {self.gen} {self.pressure} torr', loc='left')
            return lines1 + lines2
        
        time_ax = np.asarray(list(data_obj1.keys()))
        ani = FuncAnimation(fig, animate, init_func=init, frames=time_ax, blit=True)
        plt.grid(linestyle=':', linewidth=0.5)
        
        ani.save(save_dir + self.title_str + f"_{yvar1}_{yvar2}_{xvar}_time" + '.mp4', writer='ffmpeg')
        return
    def make_derivative(self,**kwargs):
        edit = kwargs.get('edit', None)
        xvar = kwargs.get('xvar', None) # string variable for x-axis, like 'r' or 'E'
        yvar = kwargs.get('yvar', None) # the parameter to plot like Te or <Z> etc
        yvar = np.asarray(yvar)
        
        if edit == None or xvar == None or yvar == None:
            print('Must define all variables for the plot')
            return
        print(f"calculating finite difference derivative of {yvar} wrt {xvar}")
        data_obj    = self.data[edit]
        frames      = np.asarray(list(data_obj.keys())) # These are points in time
        x_axis      = list(data_obj[frames[0]][xvar])
        x_axis      = np.asarray([float(x) for x in x_axis])
        midpoints = (x_axis[:-1] + x_axis[1:]) / 2
        # print(f"xaxis {x_axis}")
        # self.data['DERIVATIVES'] = {}
        self.data[f"d{yvar}_d{xvar}"] = {}
        for time in frames:
            y_axis = np.array(data_obj[time][yvar].values).astype(float)
            # print(f"yaxis {y_axis}")
            # Calculate finite differences
            dx = np.diff(x_axis)  # Differences between x values
            dy = np.diff(y_axis)  # Differences between y values
            
            # Calculate derivative dy/dx using finite differences
            numerical_derivative = dy / dx
            derivative_at_midpoints = np.interp(x_axis, midpoints, numerical_derivative)
            dummy_dict = {xvar:x_axis,f"d{yvar}_d{xvar}":derivative_at_midpoints}
            dummy_df = pd.DataFrame(dummy_dict)
            # self.data['DERIVATIVES'][time] = dummy_df
            self.data[f"d{yvar}_d{xvar}"][time] = dummy_df
        return
    # def coplot_movie(self, **kwargs):
        #this version keeps the shaded traces, messy looking plots
    #     from matplotlib.animation import FuncAnimation
    #     save_dir = kwargs.get('save_dir', self.path) 
    #     edit = kwargs.get('edit', None)
    #     xvar = kwargs.get('xvar', None)
    #     yvar1 = kwargs.get('yvar1', None)  # New parameter for the first y-variable
    #     yvar2 = kwargs.get('yvar2', None)  # New parameter for the second y-variable
    #     if edit is None or xvar is None or yvar1 is None or yvar2 is None:
    #         print('Must define all variables for the plot')
    #         return
        
    #     data_obj = self.data[edit]
    #     frames = np.asarray(list(data_obj.keys()))
    #     x_axis = list(data_obj[frames[0]][xvar])
    #     x_axis = np.asarray([float(x) for x in x_axis])
    #     y1_min = min(np.array(list(data_obj[frames[0]][yvar1])).astype(float))
    #     y1_max = 0
    #     y2_min = min(np.array(list(data_obj[frames[0]][yvar2])).astype(float))
    #     y2_max = 0
    #     for time in frames:
    #         new_y1_max = max(np.array(list(data_obj[time][yvar1])).astype(float))
    #         new_y2_max = max(np.array(list(data_obj[time][yvar2])).astype(float))
    #         if new_y1_max > y1_max:
    #             y1_max = new_y1_max
    #         if new_y2_max > y2_max:
    #             y2_max = new_y2_max
        
    #     print(y1_min)
    #     print(y1_max)
    #     print(y2_min)
    #     print(y2_max)
        
    #     fig, ax1 = plt.subplots(figsize=(5, 5))
    #     ax2 = ax1.twinx()  # Create a secondary y-axis
        
    #     lines1 = []  # List to store line objects for the first y-variable
    #     lines2 = []  # List to store line objects for the second y-variable
        
    #     ax1.set_ylim([y1_min, y1_max + 0.05 * y1_max])
    #     ax2.set_ylim([y2_min, y2_max + 0.05 * y2_max])
        
    #     ax1.set_xlim([min(x_axis), max(x_axis)])
    #     if xvar == 'energy':
    #         ax1.set_xlim([0, 2000])
        
    #     ax1.set_xlabel(xvar)
    #     ax1.set_ylabel(yvar1, color='tab:blue')
    #     ax2.set_ylabel(yvar2, color='tab:red')
        
    #     def init():
    #         for line in lines1 + lines2:
    #             line.set_data([], [])
    #         ax1.set_title('')
    #         return lines1 + lines2
    #     # Uses the fading lines 
    #     # def animate(time):
    #     #     y1 = np.array(data_obj[time][yvar1].values).astype(float)
    #     #     y2 = np.array(data_obj[time][yvar2].values).astype(float)
            
    #     #     line1, = ax1.plot(x_axis, y1, label=str(time * 1e9) + 'ns', color='blue')
    #     #     line2, = ax2.plot(x_axis, y2, label=str(time * 1e9) + 'ns', color='red')
            
    #     #     lines1.append(line1)
    #     #     lines2.append(line2)
            
    #     #     if len(lines1) > 10:
    #     #         old_line1 = lines1.pop(0)
    #     #         old_line1.remove()
    #     #     if len(lines2) > 10:
    #     #         old_line2 = lines2.pop(0)
    #     #         old_line2.remove()
            
    #     #     for i, l in enumerate(lines1[:-1]):
    #     #         l.set_alpha(0.1 + (0.9 / len(lines1)) * (i + 1))
    #     #         l.set_linestyle('--')
    #     #         l.set_color('slategrey')
    #     #     for i, l in enumerate(lines2[:-1]):
    #     #         l.set_alpha(0.1 + (0.9 / len(lines2)) * (i + 1))
    #     #         l.set_linestyle('--')
    #     #         l.set_color('lightcoral')
            
    #     #     ax1.set_title(f"Cretin {yvar1}, {yvar2} vs {xvar} at t = {time * 1e9:.2f}" + 'ns' + f' {self.gen} {self.pressure} torr', loc='left')
    #     #     return lines1 + lines2
    #     def animate(time):
    #         y1 = np.array(data_obj[time][yvar1].values).astype(float)
    #         y2 = np.array(data_obj[time][yvar2].values).astype(float)
            
    #         for line in lines1:
    #             line.remove()  # Remove all lines from the previous frame
    #         for line in lines2:
    #             line.remove()
            
    #         lines1.clear()  # Clear the list of old lines
    #         lines2.clear()
            
    #         line1, = ax1.plot(x_axis, y1, label=str(time * 1e9) + 'ns', color='blue')
    #         line2, = ax2.plot(x_axis, y2, label=str(time * 1e9) + 'ns', color='red')
            
    #         lines1.append(line1)  # Add the current lines to the list
    #         lines2.append(line2)
            
    #         ax1.set_title(f"Cretin {yvar1} vs {xvar} at t = {time * 1e9:.2f}" + 'ns' + f' {self.gen} {self.pressure} torr', loc='left')
    #         return lines1 + lines2
        
    #     time_ax = np.asarray(list(self.data[edit].keys()))
    #     ani = FuncAnimation(fig, animate, init_func=init, frames=time_ax, blit=True)
    #     plt.grid(linestyle=':', linewidth=0.5)
        
    #     ani.save(save_dir + self.title_str + f"_{yvar1}_{yvar2}_{xvar}_time" + '.mp4', writer='ffmpeg')
    #     return
    
    
def two_csd_animation(data1, data2, save_filename, **kwargs):
    # Read in the data labels for the legend, default to 1 and 2
    label1 = kwargs.get('label1', '1')
    label2 = kwargs.get('label2', '2')
    titlestr = kwargs.get('titlestr', '') # Optional additional info for the title
    print(titlestr)
    # This will make an animation of two CSD's and force them to match a common time axis
    # Define the zbar function
    def zbar(csd_list):
        zbar = 0
        for charge, pop in enumerate(csd_list):
            zbar = zbar + charge * pop
        return zbar

    # Determine the longer time axis
    time_ax1 = data1.time_data['time'] * 1e9
    csd_dict1 = {}
    for itime, time in enumerate(time_ax1):
        csd_list1 = []
        for ion in data1.ne_ions:
            ionfrac = data1.time_data[data1.yiso_labels[ion]][itime]
            csd_list1.append(ionfrac)
        csd_dict1[time] = csd_list1
            
    time_ax2 = data2.time_data['time'] * 1e9
    csd_dict2 = {}
    for itime, time in enumerate(time_ax2):
        csd_list2 = []
        for ion in data2.ne_ions:
            ionfrac = data2.time_data[data2.yiso_labels[ion]][itime]
            csd_list2.append(ionfrac)
        csd_dict2[time] = csd_list2

    print(len(time_ax1), len(time_ax2))
    if len(time_ax1) > len(time_ax2):
        print('changing object 2')
        longer_time_ax = time_ax1
        shorter_time_ax = time_ax2
        data1_csd = csd_dict1
        # print("before", csd_dict2[0])
        data2_csd = interpolate_dictionary(data2, longer_time_ax)
        print("len of data1 = " ,len(data1_csd))
        print("len of data2 = " ,len(data2_csd))

    elif len(time_ax2) > len(time_ax1):
        print('changing object 1')
        longer_time_ax = time_ax2
        shorter_time_ax = time_ax1
        data1_csd = interpolate_dictionary(data1, longer_time_ax)
        data2_csd = csd_dict2
    else: #they are the same length no interp needed
        print('no changes, same length')
        longer_time_ax = time_ax1
        data1_csd = csd_dict1
        data2_csd = csd_dict2

    # Prepare the list of figure objects
    fig_list = []
    zavg_data1,zavg_data2 = [],[]
    for itime, time in enumerate(longer_time_ax):
        print(f"itime = {itime}")
        csd_list_data1 = list(data1_csd.values())[itime]
        # print(csd_list_data1)
        csd_list_data2 = list(data2_csd.values())[itime]
        # print(csd_list_data2)
        
        fig, ax = plt.subplots()

        
        plt.plot(range(11), csd_list_data1, label=label1)
        plt.plot(range(11), csd_list_data2, label=label2)
        
        zavg_data1.append(zbar(csd_list_data1))
        zavg_data2.append(zbar(csd_list_data2))
        
        # plt.axvline(x=zavg_data1, color='red', linestyle='--', label='Z-bar Data 1')
        # plt.axvline(x=zavg_data2, color='blue', linestyle='--', label='Z-bar Data 2')
        
        plt.xlabel('Charge state')
        plt.ylabel('Fraction')
        plt.title(f'Cretin Ne Cell Comparison')
        plt.ylim([0, 1])
        plt.grid(linestyle=':', linewidth=0.5)
        plt.legend()
        plt.tight_layout()
        
        fig_list.append(fig)
        plt.close()

    # Create the animation
    anifig, ax = plt.subplots()
    line1, = ax.plot(range(11), np.zeros(11), label=label1)
    line2, = ax.plot(range(11), np.zeros(11), label=label2)
    ax.set_ylim([-0.05, 1])
    ax.set_xticks(range(11))
    ax.set_xlabel('$Ne^{+q}$')
    plt.grid(linestyle=':', linewidth=0.5)
    plt.legend()
    lines = [line1, line2]

    # def animate(itime):
    #     line1.set_ydata([data1.time_data[data1.yiso_labels[ion]][itime] for ion in data1.ne_ions])
    #     line2.set_ydata([data2.time_data[data2.yiso_labels[ion]][itime] for ion in data2.ne_ions])
    #     ax.set_title(f'Cretin Ne Cell Comparison at t = {longer_time_ax[itime]:.2f} ns')
    def animate(itime):
        csd_list_data1 = list(data1_csd.values())[itime]
        csd_list_data2 = list(data2_csd.values())[itime]
    
        line1.set_ydata(csd_list_data1)
        line2.set_ydata(csd_list_data2)
        ax.set_title(f'Cretin Ne Slab at t = {longer_time_ax[itime]:.2f} ns, $<Z_1>$ = {zavg_data1[itime]:.2f},$<Z_2>$ = {zavg_data2[itime]:.2f} '+titlestr)
        plt.legend()
    ani = FuncAnimation(anifig, animate, frames=len(longer_time_ax), blit=False)
    ani.save(save_filename, writer='ffmpeg')
    

    
def interpolate_dictionary(data1, new_x):
    time_ax1 = data1.time_data['time'] * 1e9
    csd_interp_dict = {}
    for ion in data1.ne_ions:
        ionfrac     = data1.time_data[data1.yiso_labels[ion]]
        interp_ion  = np.interp(new_x, time_ax1, ionfrac)
        csd_interp_dict[ion] = interp_ion
    csd_dict1 = {}    
    for itime, time in enumerate(new_x):
        csd_list1 = []
        for ion in data1.ne_ions:
            ionfrac = csd_interp_dict[ion][itime]
            csd_list1.append(ionfrac)
        csd_dict1[time] = csd_list1

    return csd_dict1


