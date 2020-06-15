# Extracts the real-time solar wind data from SWPC FTP servers
# Computes time lag from L1 satellites to Earth bow shock
#
# Author: Diana Morosan (2016)
# morosand@tcd.ie
#
# Based on NOAA Ovation Model Operational

import datetime
import numpy as np
#import datetime as dt

#~ import urllib2
#~ import urllib

import requests

# import json

# proxy settings for use of real time system behind a firewall
#~ proxy  = urllib2.ProxyHandler({'http':'http://webproxy.metoffice.gov.uk:8080'})
#~ opener = urllib2.build_opener(proxy)
#~ urllib2.install_opener(opener)

# realtime solar wind data url - DISCOVR combined with ACE data
urlpath = 'http://services.swpc.noaa.gov/products/solar-wind/'
mode = 'NOWCAST'

# the time from L1 solar wind sattelite to bow shock assuming the 'half-way-in-between model'
# http://omniweb.gsfc.nasa.gov/html/omni2_doc.html#shift ? URL does not work
#
# input parameters:
#   Vsw: in, required, type='numpy.ndarray'
#   jtime: in, required, type='numpy.ndarray', datetime.datetime format
#   gse_x: in, required, type='numpy.ndarray'
#   gse_y: in, required, type='numpy.ndarray'
#
# returns:
#   Array of numbers in seconds the same length as the input parameters, which represents propagation
#   from L1 Solar Wind measurement point to Earth bow shock.
#   type='numpy.ndarray'

def compute_lag(v_sw, jtime, gse_x, gse_y):
    v_earth = 30. # Earth orbital speed (km/s)
    v_parker = 428. # (km/s)
    Re = 6378. # Earth radius (km)

    n = len(v_sw)
    delta_t = []
    
    v_sw = np.asarray(v_sw)
    gse_x = np.asarray(gse_x)
    gse_y = np.asarray(gse_y) # need to be in array format to check for nans
    
    for i in range(n):
        # if i-th Vsw is NaN, then compute an average Vsw from nearest non NaN neighbors
        # assuming not at the start nor the end
        t_sw = v_sw[i] # temporary v_sw variable
        if np.isnan(t_sw) and t_sw!=0 and i!=0 and i!=n-1: #check for nans
            # check the indices of non-nan values -> ~means not a nan
            idx_left = np.argwhere( ~np.isnan(v_sw[0:i-1]) )
            n_left = len(idx_left)
            idx_right = np.argwhere( ~np.isnan(v_sw[i+1:n-1]) )
            n_right = len(idx_right)
            
            if n_left>0 and n_right>0:
                i_left = idx_left[n_left-1]
                i_right = i+1+idx_right[0]
                
                # time needs turning to datetime object:
                # differences are in datetime.timedelta in seconds
                #print 'Difference', (jtime[i]- jtime[i_left]).seconds,(jtime[i_right]-jtime[i_left]).seconds #in seconds
                t_pos = float((jtime[i]- jtime[i_left]).seconds)/float((jtime[i_right]-jtime[i_left]).seconds) #in seconds
                #print 't_pos', t_pos
                t_sw = v_sw[i_left]*t_pos + v_sw[i_right]*( 1 - t_pos ) #has to be in km/s
                #print 't_sw', t_sw

        # from OMNI website, 'half-way-in-between' model
        W = np.tan(0.5*np.arctan( t_sw/v_parker )) # both in km/s
        X = gse_x[i]*Re
        Y = gse_y[i]*Re

        delta_t.append( (X / t_sw)*( ( 1 + Y*W/X )/( 1 - v_earth*W/t_sw ) ) )

    return np.asarray(delta_t)



# get ACE data from SWPC
# solar wind L1 data from SWPC FTP servers
#
# input parameters:
#   ftp_server: swpc data server, type = 'string'
#   ftp_dir: swpc ace realtime data directory, required, type = 'string'
#   ftp_dir2: swpc ace2 location directory, required, type = 'string'
#
# returns: sw_avg
#

def swpc_get_solar_wind_realtime_data_json(urlpath, mode, time):

	#current date of realtime observations
	current_time = time
# 	current_date = current_time.strftime('%Y%m%d') # text file format
	#how long before the observation
	minutes_before = 270 # 4.5 hours;
	start_time = current_time - datetime.timedelta(minutes=minutes_before)
# 	end_time = current_time

##############################################################
# get realtime solar wind data between start_time and end_time
##############################################################

# 1 day data

	file_mag = 'mag-1-day.json'
	file_plasma = 'plasma-1-day.json'

	url_mag = urlpath + file_mag
	url_plasma = urlpath + file_plasma

	response = requests.get(url_mag)
	data_mag = response.json()

	#~ with urllib.request.urlopen(url_mag) as response:
		#~ data_mag = response.read()
		#~ data_mag = json.loads(data_mag)

	time_sw = []
	Bx = []
	By = []
	Bz = []
	n = len(data_mag)

	for i in range(1,n-1): # because first line is header
		temp = str(data_mag[i][0])
		time = (datetime.datetime(int(temp[0:4]), int(temp[5:7]), int(temp[8:10]), int(temp[11:13]), int(temp[14:16])))
		if time > start_time and time < current_time:
			# time array for solar wind data to be extracted
			time_sw.append(time)
			# magnetic field
			Bx.append(float(data_mag[i][1]))
			By.append(float(data_mag[i][2]))
			Bz.append(float(data_mag[i][3]))

	Bx = np.asarray(Bx)
	By = np.asarray(By)
	Bz = np.asarray(Bz)

    # checking for no data values and setting them to nans
    #ind_nodata = np.where(dsflag_mag!=0) #flag 9 means no data
    #Bx[ind_nodata] = np.nan
    #By[ind_nodata] = np.nan
    #Bz[ind_nodata] = np.nan
    
	response = requests.get(url_plasma)
	data_plasma = response.json()
	
	#~ response = request.urlopen(url_plasma)
	#~ data_plasma = response.read()
	#~ data_plasma = json.loads(data_plasma)

	density = []
	speed = []
	n = len(data_plasma)
    
	for i in range(1,n): # because first line is header
		temp = str(data_plasma[i][0])
		time = (datetime.datetime(int(temp[0:4]), int(temp[5:7]), int(temp[8:10]), int(temp[11:13]), int(temp[14:16])))
		if time > start_time and time < current_time:
			# density and speed
			if data_plasma[i][2] is not None and data_plasma[i][1] is not None:
#				print(i,data_plasma[i][1],data_plasma[i][2])
				density.append(float(data_plasma[i][1]))
				speed.append(float(data_plasma[i][2]))
			# break when array sizes are the same.. writing time to files differs by 1 minute
			if len(density) == len(Bx):
				break

	density = np.asarray(density)
	speed = np.asarray(speed)
	
	if len(Bx) > len(speed):
		ls = len(speed)
		Bx = Bx[0:ls]
		By = By[0:ls]
		Bz = Bz[0:ls]
# 	print( 'Array lengths', len(Bx), len(speed), len(density))

    # when no realtime solar wind data available, stop the code
#	print(len(Bx),len(density))
	if len(Bx) == 0 or len(density) == 0 or len(speed) == 0:
		print('No realtime solar wind data available... Aborting...')
		return
		# maybe use Kp index here...


    ##################################################
    # get L1 solar wind satellitle locations - DSCOVER
    ##################################################

	file_loc = 'ephemerides.json'
	url_loc = urlpath + file_loc
	
	response = requests.get(url_loc)
	data_loc = response.json()	

	#~ response = urllib.request.urlopen(url_loc)
	#~ data_loc = response.read()
	#~ data_loc = json.loads(data_loc)

	Re = 6378. # Earth radius (km) for location units in Earth radii
	current_hour = datetime.datetime(current_time.year, current_time.month, current_time.day, current_time.hour)

	time_arr = []
	n = len(data_loc)
	for i in range(1,n-1): # because first line is header
		temp = str(data_loc[i][0])
		time_arr.append( datetime.datetime(int(temp[0:4]), int(temp[5:7]), int(temp[8:10]), int(temp[11:13]), int(temp[14:16])))

    #print 'Current hour/Time array', current_hour, time_a
	if current_hour in time_arr:
		ind = time_arr.index( current_hour )
# 		time_loc = time_arr[ind]
		gse_x = float(data_loc[ind][1])/Re
		gse_y = float(data_loc[ind][2])/Re
			
		# get previous row of solar wind locations; one hour earlier for interpolation purpose
# 		before_time_loc = time_arr[ind-1]
		before_gse_x = float(data_loc[ind-1][1])/Re
		before_gse_y = float(data_loc[ind-1][2])/Re
			
    # if location at current time not available:
    # get latest available location; locations in Earth radii
	else:
		print('No real time location available, getting latest available spacecraft location...')
# 		time_loc = time_arr[-1]
		gse_x = float(data_loc[-1][1])/Re
		gse_y = float(data_loc[-1][2])/Re
			
		# get previous row of solar wind locations; one hour earlier for interpolation purpose
# 		before_time_loc = time_arr[-2]
		before_gse_x = float(data_loc[-2][1])/Re
		before_gse_y = float(data_loc[-2][2])/Re



    #############################
    # Lag to Bow Shock (OMNI way)
    #############################
    
    # interpolate L1 solar wind hourly location data onto 1 minute data
    # location array will be same size as speed array
    #
    # new function here: compute_lag

    # if before and after x coordinates are the same, no interp needed
	if gse_x == before_gse_x:
		gse_x_interp = np.zeros(len(speed)) + gse_x
	else:
		gse_x_interp = np.linspace( gse_x, before_gse_x, len(speed) )
    # if before and after y coordinates are the same, no interp needed
	if gse_y == before_gse_y:
		gse_y_interp = np.zeros(len(speed)) + gse_y
	else:
		gse_y_interp = np.linspace( gse_y, before_gse_y, len(speed) )

	delta_t = compute_lag( speed, time_sw, gse_x_interp, gse_y_interp )
    #print 'lag time', delta_t

    #######################################################################
    # 4 Hourly Averages
    # starting with:
    #   Forecast Mode: most futuristic impact time, i.e. slowest solar wind
    #   Nowcast Mode: latest data
    #######################################################################


    # seconds to impact Earth from Now
	time_diff = []

	for i in time_sw:
		time_diff.append((current_time - i).seconds) #current time changes every second
	time_diff = -1 * np.array(time_diff)

    #  set the averaging epoch time: Forecast or Nowcast
	if mode == 'NOWCAST':
		seconds_epoch = 0  #epoch is now for NOWCAST
		seconds_to_impact = time_diff
	else:
		seconds_to_impact = delta_t + time_diff[0:len(delta_t)]
		seconds_epoch = np.nanmax(seconds_to_impact) # excludes nan values

    # time of latest solar wind observation
 
	time_latest_solar_wind = time_sw[-1]

	print ('Current time ', time, '  Latest solar wind ',time_latest_solar_wind)

	
	n_hours = 4 # data averaged over 4 hourly bins
	Bx_average = []
	By_average = []
	Bz_average = []
	Bmag_avg = []
	v_avg = []
	density_avg = []
	sec_avg = []
	

	
	for i in range(n_hours):
		idx_hour = np.where( np.logical_and( seconds_to_impact >= (seconds_epoch-(i+1)*3600), seconds_to_impact<(seconds_epoch-i*3600) ) )
		if int(int(idx_hour[-1][-1])) >= len(Bx): del(idx_hour[-1][-1])  #Occasionally there will be too many elements in the array
		Bx_average.append( np.nanmean(Bx[idx_hour]) ) # excluding nans
		By_average.append( np.nanmean(By[idx_hour]) )
		Bz_average.append( np.nanmean(Bz[idx_hour]) )
		Bmag_avg.append( np.sqrt( Bx_average[i]**2 + By_average[i]**2 + Bz_average[i]**2 ))
		v_avg.append( np.nanmean(speed[idx_hour]) )
		density_avg.append( np.nanmean(density[idx_hour]) )
		sec_avg.append( np.nanmean(seconds_to_impact[idx_hour]) )

		
# 	**************  Fill in Missing Data with adjacent or with averages  ***************************
#   *       First try to fill missing data with more recent data
#   *       If themost recent hour of solar wind data  are not avaialbe, set the whole array to zeros    
#   ***********************************************************************************************	
	
# 	Bmag_avg[0] = float("nan")   #Uncomment this line to test the Kp backup option
					 
	for i in range(n_hours):
		if np.isnan(Bmag_avg[i]):
			print('Missing Data: Attempting to fill gaps for element ',i)
			if i > 0 and np.isnan(Bmag_avg[i-1]) == False:
				Bx_average[i] = Bx_average[i-1]
				By_average[i] = By_average[i-1]
				Bz_average[i] = Bz_average[i-1]
				Bmag_avg[i] = ( np.sqrt( Bx_average[i]**2 + By_average[i]**2 + Bz_average[i]**2 ))
				v_avg[i] = v_avg[i-1]
				density_avg[i] = density_avg[i-1]
				
	if np.isnan(Bmag_avg[0]) or Bmag_avg[0] == [0]: 
		Bx_average = np.zeros(n_hours)
		By_average = np.zeros(n_hours)
		Bz_average = np.zeros(n_hours)
		Bmag_avg = np.zeros(n_hours)
		v_avg = np.zeros(n_hours)
		density_avg = np.zeros(n_hours)
		forecast_time = time_latest_solar_wind
	else:
		lead_time = int(1.5e6/v_avg[-1])
		print ('Forecast Lead time (minutes)  ',int(lead_time/60.))
		forecast_time = time_sw[-1] + datetime.timedelta(seconds = lead_time)
	
	sw_avg = { 'current_time' : current_time, 'time_latest_solar_wind' : time_latest_solar_wind,  'forecast_time': forecast_time, 'Bx' : Bx_average, 'By' : By_average, 'Bz' : Bz_average, 'B_average' : Bmag_avg, 'v' : v_avg, 'ni' : density_avg }

	return sw_avg

