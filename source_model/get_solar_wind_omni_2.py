#   This program reads an AACII flat file created by the NASA OMNI High Resolution data web page at...  

# Obtain the historic solar wind data:  
#   The easiest way to do this is at the NASA OMNI web site.  
#   You can select the data you need depending on the cadence of your required output.   
# https://omniweb.gsfc.nasa.gov/form/omni_min.html for the  1 or 5 minute high resolution data

# Once you are on this web site select your preferred parameters…
# Output as a file 
# Resolution (e.g. 5 minute)
# Start and Stop datetimes
# Select the following variables
#   IMF Magnitude 
#   Bx GSE/GSM
#   By GSE
#   Bz GSE
#   Flow Speed
#   Proton Density
# You should choose a start date that is half a day prior to the time you want to start because the model 
# requires the last 4 hours of data.  Then submit your request and you should get a link to a page with all your data.

#import configparser
#~ import pyodbc 
import datetime as dt
#from dateutil.parser import parse
import numpy as np
from datetime import datetime, date, time
#import matplotlib.pyplot as plt

def get_solar_wind_omni(input_file, start_date,first_read, last_date, Bx, By,Bz, Btot, density, speed):


    initial_date = (start_date - dt.timedelta(hours = 4))

    #~ date_array = []
    
    if first_read ==1:
        
        row = input_file.readline()    #Read and discard first line in case its a header line
        
        row = input_file.readline().split()
        
        year = int(row[0])
        doy =int( row[1])
        hour = int(row[2])
        minute = int(row[3])
        
        ldate = dt.datetime(year, 1, 1, 0, 0)
        ldate = ldate + dt.timedelta(days = doy-1)
        da = ldate + dt.timedelta(hours = hour) + dt.timedelta(minutes = minute)
    


        while da < initial_date:        # Skip Forward until starting time is found

            
            row = input_file.readline().split()
            
            year = int(row[0])
            doy =int( row[1])
            hour = int(row[2])
            minute = int(row[3])
            
            ldate = dt.datetime(year, 1, 1, 0, 0)
            ldate = ldate + dt.timedelta(days = doy-1)
            da = ldate + dt.timedelta(hours = hour) + dt.timedelta(minutes = minute)
            
        
        
        while da < start_date:                                      # Read four hours of data into arrays making 1 hour averages

            row = input_file.readline().split()
            
            year = int(row[0])
            doy =int( row[1])
            hour = int(row[2])
            minute = int(row[3])
            
            ldate = dt.datetime(year, 1, 1, 0, 0)
            ldate = ldate + dt.timedelta(days = doy-1)
            da = ldate + dt.timedelta(hours = hour) + dt.timedelta(minutes = minute)


            if float(row[8]) < 10000.:
                Bx.append(float(row[5]))
                By.append(float(row[6]))
                Bz.append(float(row[7]))
                Btot.append(float(row[4]))
                speed.append(float(row[8]))
                density.append(float(row[9]))


        if len(Bx) == 0 or len(density) == 0 or len(speed) == 0:
            print('No realtime solar wind data available... Aborting...')
            return()

    else:

        da = last_date
        
        while da < start_date:                      #  **   Just roll off first 5 minutes and add 5 more minutes.  

            row = input_file.readline().split()
            
            year = int(row[0])
            doy =int( row[1])
            hour = int(row[2])
            minute = int(row[3])
            
            ldate = dt.datetime(year, 1, 1, 0, 0)
            ldate = ldate + dt.timedelta(days = doy-1)
            da = ldate + dt.timedelta(hours = hour) + dt.timedelta(minutes = minute)
            

                
            if float(row[8]) < 10000.:
                Bx.append(float(row[5]))
                By.append(float(row[6]))
                Bz.append(float(row[7]))
                Btot.append(float(row[4]))
                speed.append(float(row[8]))
                density.append(float(row[9]))



    #######################################################################
    # Averages
    #######################################################################

    num_avg = len(Bz)
#   num_hour = int((num_avg)/4)
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

    if len(speed) == 0 or len(Btot) == 0:
        print('No Data to process:\nCheck to be sure that start dates and end dates are both within your dataset')
        return()
    
    time_sw = last_date
    lead_time = int(1.5e6/(60*speed[-1]))
    print ('Lead time (minutes)  ',lead_time)
    
    forecast_time = time_sw + dt.timedelta(minutes = lead_time)
    
    sw_avg = { 'current_time' : current_time, 'time_latest_solar_wind' : last_date,  'forecast_time': forecast_time, 'Bx' : Bx_avg, 'By' : By_avg, 'Bz' : Bz_avg, 'B_average' : Btot_avg, 'v' : speed_avg, 'ni' : den_avg }
    

    
    return sw_avg,  last_date, Bx, By,Bz, Btot, density, speed


    

                
                
