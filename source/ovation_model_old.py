# ovation_realtime.py:  Python program to make a precipitation map appropriate to a
# specified date and time from season_epoch.py, with seasonal variations incorporated
# This combines all 4 aurora types, saves to text file and plots on map for northern
# hemisphere and southern hemisphere
#
# Author: Diana Morosan
# morosand@tcd.ie
#
# Based on Ovation Prime Model IDL code, 2013 by Patrick Newell

import os
os.environ["PROJ_LIB"] = "C:/Users/rodney.viereck/AppData/Local/Continuum/miniconda2/pkgs/proj4-5.2.0-ha925a31_1/Library/share"

#file:///C:/Users/rodney.viereck/AppData/Local/Continuum/miniconda2/pkgs/proj4-5.2.0-ha925a31_1/Library/share/epsg

import datetime
import numpy as np
import matplotlib 
import datetime as dt

#matplotlib.use('Agg') # for running crontabs
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as colors

#~ import urllib

#from mpl_toolkits.basemap import Basemap
#from matplotlib.mlab import griddata
from scipy.ndimage.filters import uniform_filter, gaussian_filter, maximum_filter, generic_filter
from matplotlib import ticker

from swpc_get_solar_wind_realtime_data_json import swpc_get_solar_wind_realtime_data_json
from solar_wind_coupling import sol_coup
from utilities_functions_season_epoch import prob_estimate, season_weights, mlt_bin, mlat_bin
from season_epoch_new_2 import season_epoch
from custom_cmap import make_colormap
from Hemispheric_Power import Hemispheric_Power
from write_ascii_file import write_ascii_file
from get_solar_wind_historic_data_new import get_solar_wind_historic_data
from get_solar_wind_omni import get_solar_wind_omni
from rem_out import rem_out
from ovation_plot_geomag_new import ovation_plot_geomag


urlpath = 'http://services.swpc.noaa.gov/products/solar-wind/'  # paths needed to get solar wind realtime data

Prime_path = 'C:/Docs/Python/Ovation_New_Realtime_2019/'

SW_Data_path = Prime_path + 'SW_Data/'

dmsp_path = SW_Data_path + 'OP_DMSP_data/'
guvi_path = SW_Data_path + 'OP_GUVI_data/'


Output_Path_text = Prime_path + 'Output/Text/'
Output_Path_imag = Prime_path + 'Output/Images/'

	

# forecast mode: 'FORECAST' or 'NOWCAST'
#  'NOWCAST" uses latest solar wind values to provide a short range forecast 
#  'FORECAST' Requires 3 days of 3 hourly Kp forecast for input (24 Kp Values Total)

mode = 'NOWCAST'

#**********   RT = 1  Then Real-Time  Else Historic  ********************
RT = 1			#  1 = real time

start_date = datetime.datetime(2017,9,7,0,0)   	#YYYY, MM, DD, HH, MM

end_date = datetime.datetime(2017,10,1,0,0)	

Omni = 0

running_date = start_date

cadence = 5.  						#Cadance in minutes

tsec = (end_date - start_date).total_seconds()

nloops = int((tsec/60.)/cadence)

#nloops = 1

print (nloops)


if RT != 1:
#**********   Setup for reading historic data  ************************
	print ("Running in Historic Mode from - to ",start_date, end_date)
	
	Ifile = SW_Data_path + 'sw_data_2017_Sept.dat'
	input_file = open(Ifile, 'r')
	
	ifile_rows = input_file.readlines()
	#~ fhead = input_file.readline()		# Header	
	input_file.close()
	
	
	iline = 1 			# Initiate line counter for input file	
	
	set_date = start_date
	delta_time = dt.timedelta(minutes = cadence)
	last_date = start_date
	
		#~ Creat Arrays
	
	time_a = []
	Bx = []
	By = []
	Bz= []
	Btot = []
	density = []
	speed = []
		
	first_read = 1
	
	#~ nloops = 10

if RT ==1: nloops=1


for iloop in range(nloops):	
	
	if RT == 1:

	# *************   Get Real_time Data  ***********************
		print ('Running in Real-time Mode')
		time = datetime.datetime.utcnow() # can be any time in the past 24hrs

		#~ ***************  Get Realtime Input Data ********************
		sw_avg = swpc_get_solar_wind_realtime_data_json(urlpath, mode, time)
		# ***************************************************************

		# time of latest solar wind
		#~ time_sw = sw_avg['time_latest_solar_wind']

	# running for both hemispheres
	
	else:

	# *************   Get Historic Data  ***********************	


		if Omni == 1:
			#~ ***************  Get Historic Input Data (Omni) ********************
			sw_avg,last_date, Bx, By,Bz, Btot, density, speed = get_solar_wind_omni(input_file, set_date,first_read, last_date, Bx, By,Bz, Btot, density, speed)
			# ***************************************************************
		else:
			#~ ***************  Get Historic Input Data (flat file) ********************
#			sw_avg,last_date, Bx, By,Bz, Btot, density, speed, iline = get_solar_wind_historic_data(input_file, iline, set_date,first_read, last_date, Bx, By,Bz, Btot, density, speed)
			sw_avg,last_date, Bx, By,Bz, Btot, density, speed, iline = get_solar_wind_historic_data(ifile_rows, iline, set_date, last_date, Bx, By,Bz, Btot, density, speed)
			
			# ***************************************************************
		
		first_read = 0
		set_date = set_date + delta_time

		
	#time of latest solar wind and forecast lead time
	time_sw = sw_avg['time_latest_solar_wind']
	#~ lead_time = 1.5e6/speed[-1]
	#~ print ('Lead time (minutes)  ',lead_time/60.)
	#~ forecast_time = time_sw + dt.timedelta(seconds = lead_time)


	for NS in range(1):

		day_of_year = sw_avg['current_time'].timetuple().tm_yday
		
		
	# *********************   Call Main Routine **************************
#		print ('Calling Season Epoch  ')
		#~ Bx, By, Bz, v, ni, je = season_epoch( dmsp_path, guvi_path, sw_avg, mode, day_of_year )
		je = season_epoch( dmsp_path, guvi_path, sw_avg, mode, day_of_year )		
		# ***************************************************************
		
#		print (je)
#		stop()
	#***************   Remove Outliers  *******************************
		for ii1 in range(4):			
			
			#~ plt.figure(figsize=[8,8])
			#~ ax1 = plt.subplot(2,1,1)
			#~ ax1 = plt.imshow(je[ii1,:,:])
			
			je[ii1,:,:] = rem_out(je[ii1,:,:])
			
			#~ ax1 = plt.subplot(2,1,2)
			#~ plt.imshow(je[ii1,:,:])
			#~ plt.show()
			#~ quit()
		#***************************************************************
			
		# combining all 4 auroral types
		je_combined = je[0,:,:] + je[1,:,:] + je[2,:,:]  + je[3,:,:]  
		je_sum =  je[0,:,:] + je[1,:,:] + je[2,:,:] #   Electrons Only


	# ************   Call Routine to Calculate Hemispheric Power  ******
		
		#~ for i1 in range (4): je[i1,48,80:159] = 25.
		#~ je_combined[48,80:159] = 25.
		#~ je_sum[48,80:159] = 25.
	
#		print ('jes  ',je, je_combined,je_sum)
		
		mlt_array, mlat_array, je_array, je_d, je_m, je_w, je_i, je_s, power_hemi = Hemispheric_Power ( NS, mlt_bin, mlat_bin, je, je_combined,je_sum)
	# ***************************************************************


	#*************   Tests  *******************



		#~ print ('2 ',je[0,10:20,100:110])
		
		#~ plt.imshow(je_sum)
		#~ plt.show()

		print ('hemispheric power:  %6.1f ' %power_hemi)
#		print ('Bx = ',sw_avg['Bx'])
#		print ('By = ',sw_avg['By'])
#		print ('Bz = ',sw_avg['Bz'])
#		print ('v  = ',sw_avg['v'])
#		print ('ni = ',sw_avg['ni'])
		
	#  ******************  Write Data to File  *****************************
		
		print ("Time of Day  ",time_sw)
		wf = write_ascii_file(NS,Output_Path_text,time_sw, sw_avg,mlt_array,mlat_array,je_s,je_d,je_m,je_w,je_i,power_hemi)
		
		ovation_plot_geomag(Output_Path_text,wf, Output_Path_imag+wf+'.jpg')

																											
print('Done with program')






