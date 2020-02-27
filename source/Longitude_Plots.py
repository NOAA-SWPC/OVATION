# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 08:07:24 2019

@author: rodney.viereck
"""

#Plot Absorption vs Longitude
#  Must run Ovation_HF_plot first to define the array geo_aur

import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt')        #For seperate plot window
#get_ipython().run_line_magic('matplotlib', 'inline')    #For pllots in Spyder Consoul


plat = 65
plon = 180

#pdata = Kp_9[:,plon,plat]

#t1 = Kp_9[0,:,:]

#t1 = np.roll(t1,180, axis = 0)

fig, p1 = plt.subplots()

for ip in range (6):
	ilat = 80 - ip*5
	
	print (ip,ilat)
	plt_arr = Kp_5[:,plon,ilat]
	roll_amnt = int(plon/15.)
	plt_arr = np.roll(plt_arr,roll_amnt)
	plt.plot(plt_arr, label = ilat)

plt.axis([0,24,10,0])
fig.suptitle('HF Signal\n Kp = 5    Geographic Longitude = '+str(plon))
plt.xlabel('Local Time')
plt.ylabel('HF Impact Severity (Absorption)')
plt.legend(loc='lower right',title = 'Latitude')

plt.xticks(np.arange(0,25,4))

labels = [item.get_text() for item in p1.get_xticklabels()]
labels = ['12','16','20','24','4','8','12']
p1.set_xticklabels(labels)

plt.show()