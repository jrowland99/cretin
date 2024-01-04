#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:02:25 2023

@author: jeff

run cretin 
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



def cretin_thread(directory,name): # directory points to folder containing .gen file, name is "run.gen"
    with open(os.devnull,'w') as f:
        full_command = "cd "+directory+"; cretin "+name
        t1 = Time.time()
        subprocess.call(full_command,shell=True,stdout=f)
        t2 = Time.time()
        print (f'{directory}\n'+"finish in {} seconds.\n".format((t2-t1)))
        
def run_cretin(directory,num_threads):
    command_queue = Queue()
    com_list = []
    # Walk through all directories and run all '.in' files
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
if __name__ == "__main__":
    # Check if an argument is provided
    if len(sys.argv) < 2:
        print("Error: Please provide number of threads.")
    else:
        # Get the integer argument from the command line
        arg = sys.argv[1]
    cwd = os.getcwd()
    try:
        number = int(arg)
    except ValueError:
        print("Error: Invalid integer argument.")
    else:
        # Call the function to process the integer
        run_cretin(cwd,number)