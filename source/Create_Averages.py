
import configparser
#~ import pyodbc 
import datetime as dt
from dateutil.parser import parse
import numpy as np
from datetime import datetime, date, time
import matplotlib.pyplot as plt

def create_averages(input_file, start_date,first_read, last_date, Bx, By,Bz, Btot, density, speed):


	initial_date = (start_date - dt.timedelta(hours = 4))

	#~ date_array = []
	
	if first_read ==1:

		fdatei = input_file.readline() .strip().split()
		read_date =  fdatei[0]+ ' ' + fdatei[1]
		da = datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

		while da < initial_date:  							# Skip Forward until starting time is found
			fdatei = input_file.readline() .strip().split()
			read_date =  fdatei[0]+ ' ' + fdatei[1]
			da =datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

		#~ test = []
		#~ test_date = []
		
		while da < start_date:  									# Read four hours of data into arrays making 1 hour averages
			fdatei = input_file.readline() .strip().split()
			read_date =  fdatei[0]+ ' ' + fdatei[1]
			da = datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

			#~ time_a.append(da)
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

		#~ current_end_date = start_date +  dt.timedelta(minutes = cadence)

		#~ test_date.append(current_end_date)
		#~ test.append(float(fdatei[2]))


		if len(Bx) == 0 or len(density) == 0 or len(speed) == 0:
			print('No realtime solar wind data available... Aborting...')
			#~ return
			# maybe use Kp index here...

	else:

		da = last_date
		
		while da < start_date:						#  **   Just roll off first 5 minutes and add 5 more minutes.  
			fdatei = input_file.readline() .strip().split()
			read_date =  fdatei[0]+ ' ' + fdatei[1]
			da = datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')
					
			Bx = np.roll(Bx,-1)
			By = np.roll(By,-1)
			Bz = np.roll(Bz,-1)
			Btot = np.roll(Btot,-1)
			speed = np.roll(speed,-1)
			density = np.roll(density,-1)

			Bx[-1] = (float(fdatei[2]))
			By[-1] = (float(fdatei[3]))
			Bz[-1] = (float(fdatei[4]))
			Btot[-1] = (float(fdatei[5]))
			speed[-1] = (float(fdatei[6]))
			density[-1] = (float(fdatei[7]))
			
			#~ test_date.append(db)
			#~ test.append(float(fdatei[2]))

	#######################################################################
	# Averages
	#######################################################################

	num_avg = len(Bz)
	num_hour = int((num_avg)/4)
	ihour = [0,num_avg/4, num_avg/2, 3*(num_avg/4), num_avg]
	
	Bx_avg = []
	By_avg = []
	Bz_avg = []
	Btot_avg = []
	den_avg = []
	speed_avg = []
	
	for ihour in range(4):
		Bx_avg.append(np.sum(Bx[ihour:ihour+1]))
		By_avg.append(np.sum(By [ihour:ihour+1]))
		Bz_avg.append(np.sum(Bx[ihour:ihour+1]))
		Btot_avg.append(np.sum(Btot[ihour:ihour+1]))
		den_avg.append(np.sum(density[ihour:ihour+1]))
		speed_avg.append(np.sum(speed[ihour:ihour+1]))

	last_date = da 
	
	
	#time of latest solar wind and forecast lead time
	current_time = datetime.now()
	
	time_sw = last_date
	lead_time = int(1.5e6/(60*speed[-1]))
	print ('Lead time (minutes)  ',lead_time/60.)
	
	forecast_time = time_sw + dt.timedelta(minutes = lead_time)
	
	sw_avg = { 'current_time' : current_time, 'time_latest_solar_wind' : last_date,  'forecast_time': forecast_time, 'Bx' : Bx_avg, 'By' : By_avg, 'Bz' : Bz_avg, 'B_average' : Btot_avg, 'v' : speed_avg, 'ni' : den_avg }
	
	
	return sw_avg,	last_date, Bx, By,Bz, Btot, density, speed


	

				
				
