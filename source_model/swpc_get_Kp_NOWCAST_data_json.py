# Extracts the real-time KP Forecast data from SWPC FTP servers

#
# Author: Rodney Viereck (2019)
#
# Based on NOAA Ovation Model Operational

import datetime as dt
import numpy as np

import requests

import json


def swpc_get_Kp_NOWCAST_data_json():
	
	urlpath = 'http://services.swpc.noaa.gov/products/'
	filename = 'noaa-estimated-planetary-k-index-1-minute.json'
	
	
	##############################################################
	# get Kp forecast data 
	##############################################################
	
	
	url_Kp = urlpath + filename
	
	response = requests.get(url_Kp)
	
	data_Kp = json.loads(response.text)
	
	
	numavg= 30			#Set number of 1-minute Kp values to average
	
	n = len(data_Kp)
	Kp_el = np.zeros(n)
	
	if n <= numavg: return(0)
	
	# for i in range(1,n-1): # because first line is header
	temp = str(data_Kp[n-1][0])
	f_time = (dt.datetime(int(temp[0:4]), int(temp[5:7]), int(temp[8:10]), int(temp[11:13]), int(temp[14:16])))
	
	for i1 in range (1,n-1):Kp_el[i1] = data_Kp[i1][1]
	
	Kp = np.mean(Kp_el[n-1-numavg:n-1])
	
	return(f_time,Kp,f_time)




