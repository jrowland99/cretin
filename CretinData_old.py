#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:00:15 2023

@author: jeff

functions to read a creatin .plt file for a single model
"""
import numpy as np
class CretinData:
    """
    path points to the .plt file
    times are stored in ns
    temperatures in ev
    erad
    """
    def __init__(self, path):
        self.path = path
        self.data = {}
        self.units = {}

    def _get_chunk_lines(self, header):
        with open(self.path, 'r') as f:
            chunks = f.read().split('#')
            chunks = chunks[1:]
            for chunk in chunks:
                lines           = chunk.split('\n')
                lines           = [line for line in lines if not line.startswith('$')]
                chunk_header    = lines[0]
                lines           = lines[1:]
                if chunk_header.strip() == header:
                    return lines
        return None

    def _get_time_and_values(self, lines):
        time            = []
        values          = []
        for line in lines:
            if len(line) == 0:
                continue
            t, val = map(float, line.split()[:2])
            time.append(t * 1e9)
            values.append(val)
        return np.asarray(time), np.asarray(values)

    def get_cretin_temp(self):
        lines = self._get_chunk_lines('TEV, TIV vs TIME')
        if lines:
            time, tev           = self._get_time_and_values(lines)
            tiv                 = [float(line.split()[2]) for line in lines if not len(line) == 0]
            self.data['te']     = np.asarray(tev)
            self.data['ti']     = np.asarray(tiv)
            self.data['time']   = np.asarray(time)
            self.units['time']  = 'ns'
            self.units['tev']   = 'eV'
            self.units['tiv']   = 'eV'
        else:
            print("Warning: 'TEV, TIV vs TIME' data not found in the file.")

    def get_cretin_eden(self):
        lines = self._get_chunk_lines('NE vs TIME')
        if lines:
            time, ne = self._get_time_and_values(lines)
            self.data['ne'] = np.asarray(ne)
            self.units['ne'] = '$cm^{-3}$'
        else:
            print("Warning: 'NE vs TIME' data not found in the file.")

    def get_cretin_erad(self):
        lines = self._get_chunk_lines('ERAD vs TIME')
        if lines:
            time, erad = self._get_time_and_values(lines)
            self.data['erad'] = np.asarray(erad)
            self.units['erad'] = 'ergs $cm^{-3}$'
        else:
            print("Warning: 'ERAD vs TIME' data not found in the file.")

    def load_data(self):
        self.get_cretin_temp()
        self.get_cretin_eden()
        self.get_cretin_erad()

    def get_data(self):
        return self.data
    def get_units(self):
        return self.units
