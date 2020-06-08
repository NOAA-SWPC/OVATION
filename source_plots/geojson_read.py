#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 07:02:05 2020

@author: rodney.viereck
"""

import geojson
# import json
import numpy as np
import matplotlib as mp

open_file = '/home/rodney.viereck/python/ovation_2020/Output/NOWCAST/Text/Gridded/Aurora_latest.gjson'

with open(open_file) as f:
    data = geojson.load(f)

# data = pd.read_json(open_file)
	
for key in data:
	print (key)

	
obs_time = data["Observation Time"]
print (obs_time)

values = data["coordinates"]

# npa = np.asarray(values, dtype=np.float32)
npa = np.asarray(values)

aur = np.zeros((181,360))

for i1 in range(len(npa)):
	aur[90-npa[i1][0],npa[i1][1]] = npa[i1][2]

mp.pyplot.imshow(aur)

# print(values[0,0,0])




