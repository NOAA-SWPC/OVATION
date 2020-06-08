#  Getting Data from the swdsst database
#Collects Data from the Mag table and the Plasma table to create an output of
# Date     Time      Bx      By     Bz     Bt     Speed     Dens      SatPos_x     SatPos_y

#~ import configparser
#~ import pyodbc 
#~ import datetime as dt
#~ from dateutil.parser import parse
#~ import numpy as np
#~ from datetime import datetime, date, time

def get_solar_wind_historic_data(ifile_rows, iline, set_date, last_date, Bx, By,Bz, Btot, density, speed):
	
	import configparser
	#~ import pyodbc 
	import datetime as dt
	from dateutil.parser import parse
	import numpy as np
	from datetime import datetime, date, time

	#~ start_date = datetime(2015,3,17,5,30)   	#YYYY, MM, DD
	start_date = set_date

	
	start_date_time = (start_date - dt.timedelta(hours = 4))
	
#	print ('Reading next row  ',start_date,start_date_time)

	#~ ipath = 'C:/Docs/Python/Ovation_Real_Time/SW_Data/'


	date_array = []
	
	#~ print ipath
	#~ print input_file

	#~ input_file = open(ipath+input_file, 'r')
#	fhead = input_file.readline()		# Header

#	fdatei = input_file.readline() .strip().split()
	fdatei = ifile_rows[iline] .split()
	del_time =  fdatei[0]+ ' ' + fdatei[1]
	da = datetime.strptime(del_time, '%Y-%m-%d %H:%M:%S')

	print('Reading Next Block of Data ',iline, start_date_time, da)
	
	while da >= start_date_time:  							# Skip Forward until starting time is found
		iline = iline - 60		#Go back one hour at a time
		if iline < 0: ilin = 0

#		fdatei = input_file.readline() .strip().split()
		fdatei = ifile_rows[iline] .split()
		del_time =  fdatei[0]+ ' ' + fdatei[1]
		da =datetime.strptime(del_time, '%Y-%m-%d %H:%M:%S')
#		print(iline,da,start_date_time)
		
	while da < start_date_time:  							# Skip Forward until starting time is found
		iline = iline + 1
#		fdatei = input_file.readline() .strip().split()
		fdatei = ifile_rows[iline].split()
		del_time =  fdatei[0]+ ' ' + fdatei[1]
		da =datetime.strptime(del_time, '%Y-%m-%d %H:%M:%S')
#		print(iline, set_date, start_date_time, da)

	last_date = start_date
	
	
	time_sw = []
	Bx = []
	By = []
	Bz = []
	Btot = []
	density = []
	speed = []
	gse_x = []
	gse_y = []

	while da <= start_date:  									# Read four hours of data into arrays
		iline = iline + 1
#		fdatei = input_file.readline() .strip().split()
		fdatei = ifile_rows[iline] .split()
		del_time =  fdatei[0]+ ' ' + fdatei[1]
		da = datetime.strptime(del_time, '%Y-%m-%d %H:%M:%S')
		time_sw.append(da)

		Bx.append(float(fdatei[2]))
		By.append(float(fdatei[3]))
		Bz.append(float(fdatei[4]))
		Btot.append(float(fdatei[5]))
		speed.append(float(fdatei[6]))
		density.append(float(fdatei[7]))
		#~ gse_x.append(float(fdatei[8]))
		#~ gse_y.append(float(fdatei[9]))

	Bx = np.asarray(Bx)
	By = np.asarray(By)
	Bz = np.asarray(Bz)
	Btot = np.asarray(Btot)
	density = np.asarray(density)
	speed = np.asarray(speed)
	#~ gse_x = np.asarray(gse_x)
	#~ gse_y = np.asarray(gse_y)

	#~ print (len(gse_x), len(gse_y),len(Bx))

	if len(Bx) == 0 or len(density) == 0 or len(speed) == 0:
		print('No realtime solar wind data available... Aborting...')
		#~ return
		# maybe use Kp index here...


	#######################################################################
	# 4 Hourly Weighted Averages
	#######################################################################

	num_avg = len(Bz)
	num_hour = int((num_avg)/4)
	
		
	Bx_avg = np.zeros(4)
	By_avg = np.zeros(4)
	Bz_avg = np.zeros(4)
	Btot_avg = np.zeros(4)
	den_avg = np.zeros(4)
	speed_avg = np.zeros(4)
	

	for i1 in range(4):

		Bx_avg[i1] = np.nanmean(Bx[i1*num_hour:(i1+1)*num_hour])
		By_avg[i1] = np.nanmean(By[i1*num_hour:(i1+1)*num_hour])
		Bz_avg[i1] = np.nanmean(Bz[i1*num_hour:(i1+1)*num_hour])
		Btot_avg[i1] = np.nanmean(Btot[i1*num_hour:(i1+1)*num_hour])
		den_avg[i1] = np.nanmean(density[i1*num_hour:(i1+1)*num_hour])
		speed_avg[i1] = np.nanmean(speed[i1*num_hour:(i1+1)*num_hour])
		
	
	#~ weights = 1. + ((np.arange(0,num_avg))/100.)**1.5
	#~ wm = np.sum(weights)
	#~ weights = weights/(wm)

	#~ Bx_avg = np.sum(Bx * weights)
	#~ By_avg = np.sum(By * weights)
	#~ Bz_avg = np.sum(Bx * weights)
	#~ Btot_avg = np.sum(Btot * weights)
	#~ den_avg = np.sum(density * weights)
	#~ speed_avg = np.sum(speed * weights)
	#~ gse_x_avg = np.sum(gse_x * weights)
	#~ gse_y_avg = np.sum(gse_y* weights)

	current_time = start_date
	time_latest_solar_wind = start_date

	# final solar wind data structure
	Bx_avg = Bx_avg.tolist()
	By_avg = By_avg.tolist()
	Bz_avg  = Bz_avg.tolist()
	Btot_avg = Btot_avg.tolist()
	speed_avg = speed_avg.tolist()
	den_avg = den_avg.tolist()
	
# **********  Calculate Forecast Lead Time ***************
	
	current_time = dt.datetime.now()
	
	time_sw = last_date
	lead_time = int(1.5e6/(60*speed[-1]))
	print ('Lead time (minutes)  ',lead_time)

	
	forecast_time = time_sw + dt.timedelta(minutes = lead_time)
	

	sw_avg = {  'current_time' : current_time, 'time_latest_solar_wind' : last_date,  'forecast_time': forecast_time,  'Bx' : Bx_avg, 'By' : By_avg, 'Bz' : Bz_avg, 'B_average' : Btot_avg, 'v' : speed_avg, 'ni' :den_avg }

	return sw_avg, last_date, Bx, By,Bz, Btot, density, speed, iline


