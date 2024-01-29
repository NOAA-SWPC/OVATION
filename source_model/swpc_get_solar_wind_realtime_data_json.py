# Extracts the real-time solar wind data from SWPC FTP servers
# Computes time lag from L1 satellites to Earth bow shock
#
# Author: Rodney Viereck (2020)
# rodney.viereck@noaa.gov
#
# Based on NOAA Ovation Model Operational

import datetime
import numpy as np
import pandas as pd

import requests

import json


# realtime solar wind data url - DISCOVR combined with ACE data
urlpath = 'http://services.swpc.noaa.gov/products/solar-wind/'

def swpc_get_solar_wind_realtime_data_json(urlpath, mode, time):

    #current date of realtime observations
    current_time = time


##############################################################
# get realtime solar wind data between start_time and end_time
##############################################################

# 1 day data

    file_mag = 'mag-1-day.json'
    file_plasma = 'plasma-1-day.json'

    url_mag = urlpath + file_mag
    url_plasma = urlpath + file_plasma
    
    
    magfile = pd.read_json(url_mag)
    new_header = magfile.iloc[0]
    magfile = magfile[1:]
    magfile.columns = new_header
    magfile['Datetime'] = pd.to_datetime(magfile['time_tag'])
    magfile = magfile.set_index('Datetime') 
    
    del magfile['time_tag']                     
    
    plasmafile = pd.read_json(url_plasma)
    new_header = plasmafile.iloc[0]
    plasmafile = plasmafile[1:]
    plasmafile.columns = new_header
    plasmafile['Datetime'] = pd.to_datetime(plasmafile['time_tag'])
    plasmafile = plasmafile.set_index('Datetime')   
    
    del plasmafile['time_tag']  
    
    combined = pd.concat([magfile, plasmafile], axis = 1, sort = True)
    combined = combined.astype(np.float16)
    combined = combined.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    combined["temperature"] = combined["temperature"].astype(int)

    time_latest_solar_wind = combined.index[-1]
    time_epoch = current_time
    

    
    
    # Set up four 1-hour bins
    td = datetime.timedelta(hours = 1)  
    n_hours = 4 

    Bx_average = []
    By_average = []
    Bz_average = []
    Bmag_avg = []
    v_avg = []
    density_avg = []
    sec_avg = []
    
    for i1 in range(n_hours):
    
        et = (n_hours - i1)*3600.
        st = (n_hours - i1+1)*3600.
        mask = (combined.index >= time_epoch-(i1+1)*td) & (combined.index < time_epoch - i1*td) 
                             
        c0 = combined.loc[mask]
        c1 = c0.mean()
        
        Bx_average.append(c1["bx_gsm"])
        By_average.append(c1["by_gsm"])
        Bz_average.append(c1["bz_gsm"])
        Bmag_avg.append(c1["bt"])
        v_avg.append(c1["speed"])
        density_avg.append(c1["density"])

        
#   print(Bx_average)
#   print(By_average)
#   print(Bz_average)
#   print(Bmag_avg)
#   print(v_avg)
#   print(density_avg)

    
#   **************  Fill in Missing Data with adjacent or with averages  ***************************
#   *       First try to fill missing data with more recent data
#   *       If themost recent hour of solar wind data  are not avaialbe, set the whole array to zeros    
#   *********************************************************************************************** 
    
#   Bmag_avg[0] = float("nan")   #Uncomment this line to test the Kp backup option
                     
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
                
    lead_time = int(1.5e6/v_avg[-1])
    print ('Forecast Lead time (minutes)  ',int(lead_time/60.))
    forecast_time = current_time + datetime.timedelta(seconds = lead_time)
    
    sw_avg = { 'current_time' : current_time, 'time_latest_solar_wind' : time_latest_solar_wind,  'forecast_time': forecast_time, 'Bx' : Bx_average, 'By' : By_average, 'Bz' : Bz_average, 'B_average' : Bmag_avg, 'v' : v_avg, 'ni' : density_avg }

    return sw_avg

