# Extracts the real-time KP Forecast data from SWPC FTP servers

#
# Author: Rodney Viereck (2019)
#
# Based on NOAA Ovation Model Operational

import datetime as dt
import numpy as np

import requests

import json


def swpc_get_Kp_forecast_data_json():
    
    urlpath = 'http://services.swpc.noaa.gov/products/'
    filename = 'noaa-planetary-k-index-forecast.json'
    
    
    ##############################################################
    # get Kp forecast data 
    ##############################################################
    
    
    url_Kp = urlpath + filename
    
    response = requests.get(url_Kp)
    
    data_Kp = json.loads(response.text)

    
    time_Kp = []
    Kp = []
    
    n = len(data_Kp)
    Mode = []
    
    for i in range(1,n-1): # because first line is header
        temp = str(data_Kp[i][0])
        f_time = (dt.datetime(int(temp[0:4]), int(temp[5:7]), int(temp[8:10]), int(temp[11:13]), int(temp[14:16])))
        Kp_el = str(data_Kp[i][1])
        Mode_old = Mode
        Mode = str(data_Kp[i][2])
        
        if Mode_old == "estimated" and Mode == "predicted":     #Estimate time that Kp Forecast was made.
            issue_time = f_time - dt.timedelta(hours = 11.5)
    
        if Mode != 'observed':
            time_Kp = np.append(time_Kp,f_time)
            Kp = np.append(Kp,Kp_el)
            
    
#   for ii in range(len(Kp)): print(ii,time_Kp[ii], Kp[ii])
    
    
    return(time_Kp,Kp,issue_time)


