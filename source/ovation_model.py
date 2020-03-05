# ovation_realtime.py:  Python program to make a precipitation map appropriate to a
# specified date and time from season_epoch.py, with seasonal variations incorporated
# This combines all 4 aurora types, saves to text file and plots on map for northern
# hemisphere and southern hemisphere
#
# Author: Diana Morosan
# morosand@tcd.ie

#  Modified by Rodney Viereck
#  rodney.viereck@noaa.gov
#
# Based on Ovation Prime 2013 Model IDL code, by Patrick Newell

#import os
#os.environ["PROJ_LIB"] = "C:/Users/rodney.viereck/AppData/Local/Continuum/miniconda2/pkgs/proj4-5.2.0-ha925a31_1/Library/share"

#file:///C:/Users/rodney.viereck/AppData/Local/Continuum/miniconda2/pkgs/proj4-5.2.0-ha925a31_1/Library/share/epsg

import os
import datetime
# import numpy as np
#import matplotlib
import datetime as dt
import math as math
# import configparser

from swpc_get_solar_wind_realtime_data_json import swpc_get_solar_wind_realtime_data_json
from swpc_get_Kp_forecast_data_json import swpc_get_Kp_forecast_data_json
# from solar_wind_coupling import sol_coup
from utilities_functions_season_epoch import prob_estimate, season_weights, mlt_bin, mlat_bin
from season_epoch_Kp import season_epoch_Kp
#from custom_cmap import make_colormap
from Hemispheric_Power import Hemispheric_Power
from write_ascii_file import write_ascii_file
from get_solar_wind_historic_data_new import get_solar_wind_historic_data
from get_solar_wind_omni import get_solar_wind_omni
from rem_out import rem_out
from ovation_plot_geomag import ovation_plot_geomag

#***************************   Set Run Mode  *******************************

# forecast mode: 'FORECAST' or 'NOWCAST' or 'HISTORIC'
#  'HISTORIC' loops through specific dates in the past... must have access to a file with
#	historic solar wind data.  Could be an OMNI file or it could be an ASCII file with
#   just the required parameters of Velocity, Dentisy, Bz

#  'NOWCAST" uses latest solar wind values to provide a short range forecast.
#	Reads data from real-time JSON files at SWPC

#  'FORECAST' Requires 3 days of 3 hourly Kp forecast for input (24 Kp Values Total)
#	Reads these data from the SWPC online database


mode = os.environ.get('mode', 'NOWCAST')

#************************    Set Variables  *************************
#import configparser1

Home_path = os.environ.get('Home_path','.')
SW_Data_path = os.environ.get('SW_Data_path', '../SW_Data/')
dmsp_path = os.environ.get('dmsp_path', '../SW_Data/OP_DMSP_data/')
guvi_path = os.environ.get('guvi_path', '../SW_Data/OP_GUVI_data/')

Output_path = os.environ.get('Output_path', '../Output/')

input_file_historic = os.environ.get('input_file_historic', 'sw_data_2017_Sept.dat')

urlpath = os.environ.get('urlpath', 'http://services.swpc.noaa.gov/products/solar-wind/')

start_date = os.environ.get('start_date', '09-27-2017 00:00')
start_date = datetime.datetime.strptime(start_date, '%m-%d-%Y %H:%M')
end_date = os.environ.get('end_date', '09-30-2017 00:00')
end_date = datetime.datetime.strptime(end_date, '%m-%d-%Y %H:%M')
cadence = os.environ.get('cadence', 60)

#Flags

NorthSouth = os.environ.get('NorthSouth', True)  #If True then plot both hemisphere
Omni = os.environ.get('Omni', True)  #If, True, then use OMNI data... otherwise use a flat file
image_output = os.environ.get('image_output', False)  #If True, then output images



Output_Path_text = Output_path + mode + '/Text/'
# Create dirs if missing
os.makedirs(Output_Path_text + 'North/', exist_ok=True)
print("making directory path {}".format(Output_Path_text +'North/'))
os.makedirs(Output_Path_text + 'South/', exist_ok=True)
print("making directory path {}".format(Output_Path_text +'South/'))
Output_path_images = Output_path + mode + '/Images/'
os.makedirs(Output_path_images, exist_ok=True)
print("making directory path {}".format(Output_path_images))

#Clean out data from previous forecast

if mode == 'FORECAST': 	os.system('rm ' + Output_Path_text + '/North/*.txt')


nloops = 1

NS_loop = 1
if NorthSouth == True:
	NS_loop = 2  #if NS = 1 then just do norther hemisphere.   Make NS = 2 for both hemispheres

if mode == 'FORECAST': NS_loop = 1       #Don't calculate southern hemisphere

time = datetime.datetime.utcnow() # Set Current Time


Kp_1 = 0

if mode == 'HISTORIC':

# 	sd = start_date
# 	ed = end_date
# 	start_date = datetime.datetime(int(sd[0:4]),int(sd[5:7]),int(sd[8:10]),int(sd[11:13]),int(sd[14:16]))
# 	end_date = datetime.datetime(int(ed[0:4]),int(ed[5:7]),int(ed[8:10]),int(ed[11:13]),int(ed[14:16]))
							#Cadence in minutes)
	Ifile = SW_Data_path + input_file_historic

	print ("Running in Historic Mode from - to ",start_date, end_date)
	print ("Running at ", cadence, " minute cadence")

	tsec = (end_date - start_date).total_seconds()
	nloops = int((tsec/60.)/cadence)

	input_file = open(Ifile, 'r')

	if Omni != 1:
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
# 	print ('Running in 3-Day Forecast Mode')
	nloops = 24   #Loop through 3 days

	time_Kp, Kp = swpc_get_Kp_forecast_data_json()#
	sw_avg = { 'current_time' : time_Kp[0], 'time_latest_solar_wind' : time_Kp[0],  'forecast_time': time_Kp[-1], 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }

#	To force a run status, uncomment the following lines and comment the previous two

#	cur_time = dt.datetime.now()
#	lsw_time = dt.datetime(2019,3,21, 12, 00, 00)
#	for_time = lsw_time + dt.timedelta(seconds = 900)
#	time_Kp = np.array([for_time])
#	Kp = np.array([9])
#	sw_avg = { 'current_time' : cur_time, 'time_latest_solar_wind' : lsw_time,  'forecast_time': for_time, 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }


# 	nloops = len(Kp)   #Loop more


# print ('Number of Loops  ',nloops)



for iloop in range(nloops):

	if mode == 'NOWCAST':

	# *************   Get Real_time Data  ***********************
# 		print ('Running in Real-time Mode')
# 		time = datetime.datetime.utcnow() # can be any time in the past 24hrs

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



	for NS in range(NS_loop):


		if NS == 1:				#For southern hemisphere
			day_of_year = day_of_year +182
			if day_of_year > 365: day_of_year = day_of_year - 365

	# *********************   Call Main Routine **************************

		je = season_epoch_Kp( dmsp_path, guvi_path, sw_avg, mode, day_of_year,Kp_1)
		# ***************************************************************

	#***************   Remove Outliers  *******************************
		for ii1 in range(4):

			je[ii1,:,:] = rem_out(je[ii1,:,:])

		#***************************************************************

# 		# combining all 4 auroral types
# 		je_combined = je[0,:,:] + je[1,:,:] + je[2,:,:]  + je[3,:,:]
# 		je_sum =  je[0,:,:] + je[1,:,:] + je[2,:,:] #   Electrons Only


	# ************   Call Routine to Calculate Hemispheric Power  ******


# 		mlt_array, mlat_array, je_array, je_d, je_m, je_w, je_i, je_s, power_hemi = Hemispheric_Power ( NS, mlt_bin, mlat_bin, je, je_combined,je_sum)
		mlt_array, mlat_array, je_array, je_d, je_m, je_w, je_i, power_hemi = Hemispheric_Power ( NS, mlt_bin, mlat_bin, je)

	# ***************************************************************




		if NS != 1:
			print ('North hemispheric power:  %6.1f ' %power_hemi)
		else:
			print ('South hemispheric power:  %6.1f ' %power_hemi)

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

			time_lab=time
			time_for = sw_avg['forecast_time']
		else:
			time_for = time_Kp[iloop]
			time_lab = time_for

		print ("Current Time ",time)
		print ("Time of Last Solar Wind ",time_sw)
		print ("Forecast Time", time_for)
		print ('NS = ',NS)



		wf = write_ascii_file(mode,NS,Output_Path_text,time_sw, time_for, time_lab, mlt_array,mlat_array,je_d,je_m,je_w,je_i,power_hemi,Kp_1)

		#if image_output == False: ovation_plot_geomag(Output_Path_text,wf, Output_path_images+wf+'.jpg')
		print("image_output, Output_Path_text: {}, 2nd arg: {}".format(Output_Path_text, Output_path_images+wf+'.png'))
		if image_output == False: ovation_plot_geomag(Output_Path_text,wf, Output_path_images+wf+'.png')


# print('Done with program')
