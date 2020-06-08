#  Getting Data from the swdsst database
#Collects Data from the Mag table and the Plasma table to create an output of
# Date     Time      Bx      By     Bz     Bt     Speed     Dens      SatPos_x     SatPos_y

#~ import configparser
#~ import pyodbc 
#~ import datetime as dt
#~ from dateutil.parser import parse
#~ import numpy as np
#~ from datetime import datetime, date, time

def get_solar_wind_historic_data(Input_file, set_date,end_date):
	
	import configparser
	import pyodbc 
	import datetime as dt
	from dateutil.parser import parse
	import numpy as np
	from datetime import datetime, date, time

	#~ start_date = datetime(2015,3,17,5,30)   	#YYYY, MM, DD
	start_date = set_time

	#~ print start_date

	start_date_time = (start_date - dt.timedelta(hours = 4))

	ipath = 'C:/Docs/Python/Ovation_Real_Time/SW_Data/'
	ifile = Input_file


	date_array = []

	input_file = open(ipath+ifile, 'r')
	fhead = input_file.readline()		# Header

	fdatei = input_file.readline() .strip().split()
	dt =  fdatei[0]+ ' ' + fdatei[1]
	da = datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')

	while da <= start_date_time:  							# Skip Forward until starting time is found
		fdatei = input_file.readline() .strip().split()
		dt =  fdatei[0]+ ' ' + fdatei[1]
		da =datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')

	Bx = []
	By = []
	Bz = []
	Btot = []
	density = []
	speed = []

	while da <= start_date:  									# Read four hours of data into arrays
		fdatei = input_file.readline() .strip().split()
		dt =  fdatei[0]+ ' ' + fdatei[1]
		da = datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')

		Bx.append(float(fdatei[2]))
		By.append(float(fdatei[3]))
		Bz.append(float(fdatei[4]))
		Btot.append(float(fdatei[5]))
		speed.append(float(fdatei[6]))
		density.append(float(fdatei[7]))
		#~ gse_x=(float(fdatei[8]))
		#~ gse_y=(float(fdatei[9]))

	Bx = np.asarray(Bx)
	By = np.asarray(By)
	Bz = np.asarray(Bz)
	Btot = np.asarray(Btot)
	density = np.asarray(density)
	speed = np.asarray(speed)

	if len(Bx) == 0 or len(density) == 0 or len(speed) == 0:
		print('No realtime solar wind data available... Aborting...')
		#~ return
		# maybe use Kp index here...



	#######################################################################
	# Weighted Averages
	#######################################################################

	num_avg = len(Bz)
	num_hour = int((num_avg)/4)
	
	current_time=[]
	Bx_avg = []
	By_avg = []
	Bz_avg = []
	Btot_avg = []
	den_avg = []
	speed_avg = []
	
	#~ Bx_a= np.zeros(4)
	#~ By_a = np.zeros(4)
	#~ Bz_a = np.zeros(4)
	#~ Btot_a = np.zeros(4)
	#~ den_a = np.zeros(4)
	#~ speed_a = np.zeros(4)
	

	#~ for i1 in range(4):

		#~ Bx_a[i1] = np.nanmean(Bx[i1*num_hour:(i1+1)*num_hour])
		#~ By_a[i1] = np.nanmean(By[i1*num_hour:(i1+1)*num_hour])
		#~ Bz_a[i1] = np.nanmean(Bz[i1*num_hour:(i1+1)*num_hour])
		#~ Btot_a[i1] = np.nanmean(Btot[i1*num_hour:(i1+1)*num_hour])
		#~ den_a[i1] = np.nanmean(density[i1*num_hour:(i1+1)*num_hour])
		#~ speed_a[i1] = np.nanmean(speed[i1*num_hour:(i1+1)*num_hour])
		
	
	weights = 1. + ((np.arange(0,num_avg))/100.)**1.5
	wm = np.sum(weights)
	weights = weights/(wm)

	time_sw.append( start_date)
	Bx_avg.append(np.sum(Bx * weights))
	By_avg.append(np.sum(By * weights))
	Bz_avg.append(np.sum(Bx * weights))
	Btot_avg.append(np.sum(Btot * weights))
	den_avg.append(np.sum(density * weights))
	speed_avg.append(np.sum(speed * weights))
	#~ gse_x_avg.append(np.sum(gse_x * weights))
	#~ gse_y_avg.append(np.sum(gse_y* weights))

	#~ time_latest_solar_wind = start_date
	
	while da <= end_date:
		current_end_date = current_end_date + dt.timedelta(minutes = cadence)
		while db <= current_end_date:
			
			fdatei = input_file.readline() .strip().split()
			dt =  fdatei[0]+ ' ' + fdatei[1]
			db = datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')
			
			Bx = shift(Bx,1)
			By = shift(By,1)
			Bz = shift(Bz,1)
			Btot = shift(Btot,1)
			speed = shift(speed,1)
			density = shift(density,1)

			Bx[-1] = (float(fdatei[2]))
			By[-1] = (float(fdatei[3]))
			Bz[-1] = (float(fdatei[4]))
			Btot[-1] = (float(fdatei[5]))
			speed[-1] = (float(fdatei[6]))
			density[-1] = (float(fdatei[7]))
			
		da = db

		time_sw.append(current_end_date)
		Bx_avg.append(np.sum(Bx * weights))
		By_avg.append(np.sum(By * weights))
		Bz_avg.append(np.sum(Bx * weights))
		Btot_avg.append(np.sum(Btot * weights))
		den_avg.append(np.sum(density * weights))
		speed_avg.append(np.sum(speed * weights))
		#~ gse_x_avg.append(np.sum(gse_x * weights))
		#~ gse_y_avg.append(np.sum(gse_y* weights))
				
				

	#~ # final solar wind data structure
	#~ Bx_avg = Bx_avg.tolist()
	#~ By_avg = By_avg.tolist()
	#~ Bz_avg  = Bz_avg.tolist()
	#~ Btot_avg = Btot_avg.tolist()
	#~ speed_avg = speed_avg.tolist()
	#~ den_avg = den_avg.tolist()
	

	sw_avg = { 'current_time' : current_time, 'time_latest_solar_wind' : time_latest_solar_wind, 'Bx' : Bx_avg, 'By' : By_avg, 'Bz' : Bz_avg, 'B_average' : Btot_avg, 'v' : speed_avg, 'ni' :den_avg }

	return sw_data


