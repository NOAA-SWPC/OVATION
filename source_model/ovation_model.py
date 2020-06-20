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
from write_HP_file_model import write_HP_file

#***************************   Set Run Mode  *******************************

# OVATION mode: 'FORECAST' or 'NOWCAST' or 'HISTORIC'
#  'HISTORIC' loops through specific dates in the past... must have access to a file with
#	historic solar wind data.  Could be an OMNI file or it could be an ASCII file with
#   just the required parameters of Velocity, Dentisy, Bz

#  'NOWCAST" uses latest solar wind values to provide a short range forecast.
#	Reads data from real-time JSON files at SWPC

#  'FORECAST' Requires upt to 3 days of 3 hourly Kp forecast for input (24 Kp Values Total)
#	Reads these data from the SWPC online database


mode = os.environ.get('mode', 'NOWCAST')

print("mode is: {}".format(mode))

#************************    Set Variables  *************************

# Set Paths

Home_path = os.environ.get('Home_path','.')
SW_Data_path = os.environ.get('SW_Data_path', '../SW_Data/')
dmsp_path = os.environ.get('dmsp_path', '../SW_Data/OP_DMSP_data/')
guvi_path = os.environ.get('guvi_path', '../SW_Data/OP_GUVI_data/')
header_path = os.environ.get('header_path', '../SW_Data/Header_Text/')

Output_path = os.environ.get('Output_path', '../output/')

input_file_historic = os.environ.get('input_file_historic', 'Historic_SW_Data/sw_data_2017_Sept.dat')

urlpath = os.environ.get('urlpath', 'http://services.swpc.noaa.gov/products/solar-wind/')


start_date = os.environ.get('start_date', '09-03-2017 00:00')
start_date = datetime.datetime.strptime(start_date, '%m-%d-%Y %H:%M')
end_date = os.environ.get('end_date', '09-10-2017 01:00')
end_date = datetime.datetime.strptime(end_date, '%m-%d-%Y %H:%M')
cadence = os.environ.get('cadence', 30)   #Cadence in Minutes

#Set Run Option Flags

Omni = os.environ.get('Omni', False)  #If, True, then use OMNI data... otherwise use a flat file
NorthSouth = os.environ.get('NorthSouth', True)  #If True then plot both hemisphereal
HPI_output = os.environ.get('HPI_output', True)  # If True then output Hemispheric Power Index to a file
aurora_output = os.environ.get('aurora_output', True)  # If True, the output the aurora ASCII file

if mode == 'FORECAST': HPI_output = False
num_forecast_days = os.environ.get('num_forecast_days', 1)  #Set the number of days to forecast (must be less than 3)

Output_Path_text = Output_path + mode + '/model_output/'
HP_Output_path = Output_path + mode + '/ovation_products/hpi_text/'

# Create directories if missing

os.makedirs(Output_Path_text + 'north/', exist_ok=True)
# print("making directory path (if necessary) {}".format(Output_Path_text +'north/'))
os.makedirs(Output_Path_text + 'south/', exist_ok=True)
# print("making directory path (if necessary) {}".format(Output_Path_text +'south/'))

if HPI_output:
 	os.makedirs(Output_path + mode + '/ovation_products/hpi_text/', exist_ok=True)
 	# print("making directory path (if necessary) {}".format(Output_path + mode + '/ovation_products/hpi_text/'))


#Clean out FORECAST data from previous forecast

if mode == 'FORECAST':
	os.system('rm ' + Output_Path_text + 'north/*.txt')
	os.system('rm ' + Output_Path_text + 'south/*.txt')

if mode == 'HISTORIC':
	print('Do you want to delete existing files from previous runs?')

	answer = None
	while answer not in ("yes", "no"):
		answer = input("Enter yes or no: ")
		if answer == "yes":
			os.system('rm ' + Output_Path_text + 'north/*.txt')
			os.system('rm ' + Output_Path_text + 'south/*.txt')
			os.system('rm ' + HP_Output_path + '*.txt')
		elif answer == "no":
			print("May overwrite or append to existing files")
		else:
			print("Please enter yes or no.")
			
# os.system('rm ' + Output_Path_text + 'north/*.txt')
# os.system('rm ' + Output_Path_text + 'south/*.txt')
# os.system('rm ' + HP_Output_path + '*.txt')

#  nloops is the number of times the model will be run (times 2 for North and South)
#	  HISTORIC mode it will depend on length of period and cadence
# 	  NOWCAST mode nloops will be 1
# 	  FORECAST mode nloops will be 8 times per day times the number of days in the forecast
	
nloops = 1

NS_loop = 1
if NorthSouth == True:
	NS_loop = 2  #if NS = 1 then just do norther hemisphere.   Make NS = 2 for both hemispheres

time_now = datetime.datetime.utcnow() # Set Current Time
time_now = time_now.replace(second=0, microsecond=0)

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
# 	print ('Running in Multi-Day Forecast Mode')

	time_Kp, Kp , time_issue = swpc_get_Kp_forecast_data_json()
	
	sw_avg = { 'current_time' : time_now, 'time_latest_solar_wind' : time_issue,  'forecast_time': time_Kp[-1], 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }

	nloops = num_forecast_days * 8
	if nloops > len(Kp):
		nloops = len(Kp)
		
# 	nloops = 5 #override for testing

#	To force a run status, uncomment the following lines and comment the previous two

#	cur_time = dt.datetime.now()
#	lsw_time = dt.datetime(2019,3,21, 12, 00, 00)
#	for_time = lsw_time + dt.timedelta(seconds = 900)
#	time_Kp = np.array([for_time])
#	Kp = np.array([9])
#	sw_avg = { 'current_time' : cur_time, 'time_latest_solar_wind' : lsw_time,  'forecast_time': for_time, 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }


#   *************  Done with setup....   start running the model  ***************************888888


for iloop in range(nloops):

	if mode == 'NOWCAST':

	# *************   Get Real_time Data  ***********************
# 		print ('Running in Real-time Mode')
# 		time = datetime.datetime.utcnow() # can be any time in the past 24hrs

		#~ ***************  Get Realtime Input Data ********************
		sw_avg = swpc_get_solar_wind_realtime_data_json(urlpath, mode, time_now)
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
		time_sw = time_issue
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



	# ************   Call Routine to Calculate Hemispheric Power  ******


# 		mlt_array, mlat_array, je_array, je_d, je_m, je_w, je_i, je_s, power_hemi = Hemispheric_Power ( NS, mlt_bin, mlat_bin, je, je_combined,je_sum)
		mlt_array, mlat_array, je_array, je_d, je_m, je_w, je_i, power_hemi = Hemispheric_Power ( NS, mlt_bin, mlat_bin, je)

	# ***************************************************************




		if NS != 1:
			print ('North hemispheric power:  %6.1f ' %power_hemi)
		else:
			print ('South hemispheric power:  %6.1f ' %power_hemi)
			

	#  ******************  Write Data to Files  *****************************

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

# 			time_lab=time_now
			time_lab = time_sw
			time_for = sw_avg['forecast_time']
			
		else:
			time_for = time_Kp[iloop]
			time_lab = time_for
			
		if mode == 'NOWCAST': time_lab = time_now

		if HPI_output: write_HP_file(power_hemi, HP_Output_path, header_path,time_lab,time_for, NS, NorthSouth, mode)

		lbl_1 = "Time of Last Solar Wind "
		if mode == 'FORECAST': lbl_1 = "Time Kp Forecast was Issued "
# 		print ("Current Time ",time_now)
# 		print (lbl_1 ,time_sw)
		print ("Forecast Time", time_for)
		print ('NS = ',NS)

		if aurora_output: 
			opath,ofile = write_ascii_file(mode,NS,Output_Path_text,time_sw, time_for, time_now, mlt_array,mlat_array,je_d,je_m,je_w,je_i,power_hemi,Kp_1,sw_avg)


