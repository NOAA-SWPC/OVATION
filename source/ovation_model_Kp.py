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
import math as math

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

	

# forecast mode: 'FORECAST' or 'NOWCAST' or 'HISTORIC'
#  'HISTORIC' loops through specific dates
#  'NOWCAST" uses latest solar wind values to provide a short range forecast 
#  'FORECAST' Requires 3 days of 3 hourly Kp forecast for input (24 Kp Values Total)

mode = 'FORECAST'

#**********   RT = 1  Then Real-Time  Else Historic  ********************
nloops = 24
Omni = 0
Kp_1 = 0

if mode == 'HISTORIC':

	start_date = datetime.datetime(2017,9,9,0,0)   	#YYYY, MM, DD, HH, MM
	
	end_date = datetime.datetime(2017,10,1,6,0)	
	
	cadence = 60.  						#Cadance in minutes
	
	tsec = (end_date - start_date).total_seconds()
	
	nloops = int((tsec/60.)/cadence)
	
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



#**************   Reading Forecasted Kp  *************

if mode == 'FORECAST':
	# *************   Get Forecasted Kp Data  ***********************
	print ('Running in 3-Day Forecast Mode')
	
#	time_Kp, Kp = swpc_get_Kp_forecast_data_json()
#	
#	sw_avg = { 'current_time' : time_Kp[0], 'time_latest_solar_wind' : time_Kp[0],  'forecast_time': time_Kp[-1], 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }

#	To force a run status, uncomment the following lines and comment the previous
	
	cur_time = dt.datetime.now()
	lsw_time = dt.datetime(2019,3,21, 12, 00, 00)
	for_time = lsw_time + dt.timedelta(seconds = 900)
	time_Kp = np.array([for_time])
	
	Kp = np.array([9])

	sw_avg = { 'current_time' : cur_time, 'time_latest_solar_wind' : lsw_time,  'forecast_time': for_time, 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }


	nloops = len(Kp)


print ('Number of Loops  ',nloops)



for iloop in range(nloops):	
	
	if mode == 'NOWCAST':

	# *************   Get Real_time Data  ***********************
		print ('Running in Real-time Mode')
		time = datetime.datetime.utcnow() # can be any time in the past 24hrs

		#~ ***************  Get Realtime Input Data ********************
		sw_avg = swpc_get_solar_wind_realtime_data_json(urlpath, mode, time)
		# ***************************************************************
		
	if mode == 'HISTORIC':

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

		

		
#	lead_time = 1.5e6/sw_avg['v'][-1]
#	print ('Lead time (minutes)  ',lead_time/60.)
#	forecast_time = time_sw + dt.timedelta(seconds = lead_time)

	if mode == 'FORECAST': 
		Kp_1 = int(Kp[iloop])
		time_sw = time_Kp[0]
		time_for = time_Kp[iloop]
		day_of_year = time_for.timetuple().tm_yday
		
	else:
					
		day_of_year = sw_avg['current_time'].timetuple().tm_yday
		

		
	for NS in range(1):

		
	# *********************   Call Main Routine **************************
#		print ('Calling Season Epoch  ')
		#~ Bx, By, Bz, v, ni, je = season_epoch( dmsp_path, guvi_path, sw_avg, mode, day_of_year )
		je = season_epoch_Kp( dmsp_path, guvi_path, sw_avg, mode, day_of_year,Kp_1)		
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
	
		if mode != 'FORECAST':

			#time of latest solar wind
			time_sw = sw_avg['time_latest_solar_wind']
			
				#  Calculate Kp  from Hp		
			y0 = 6.9498
			A1 = 2.04984
			t1 = 1.8400
		
			lnval = ((power_hemi/A1)-y0)
			if lnval <= 0.0: lnval = 1e-5
			Kp_1= int(t1*math.log(lnval))
			if Kp_1 < 0.0: Kp_1 = 0.
		
			time_lab=time_sw
			time_for = sw_avg['forecast_time']
		else:
			time_for = time_Kp[iloop]
			time_lab = time_for
			
		
		print ("Time of Day  ",time_sw)
		print ("Forecast Time", time_for)
		
		wf = write_ascii_file(NS,Output_Path_text,time_sw, time_for, time_lab, mlt_array,mlat_array,je_s,je_d,je_m,je_w,je_i,power_hemi,Kp_1)
		
		ovation_plot_geomag(Output_Path_text,wf, Output_Path_imag+wf+'.jpg')

																											
print('Done with program')






