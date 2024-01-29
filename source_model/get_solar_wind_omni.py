
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
    
        
        #~ print (year, doy, hour, minute)
        #~ print (ldatetime)
        #~ quit()

        #~ fdatei = input_file.readline() .strip().split()
        #~ read_date =  fdatei[0]+ ' ' + fdatei[1]
        #~ da = datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

        while da < initial_date:        # Skip Forward until starting time is found
            #~ fdatei = input_file.readline() .strip().split()
            #~ read_date =  fdatei[0]+ ' ' + fdatei[1]
            #~ da =datetime.strptime(read_date, '%Y-%m-%d %H:%M:%S')

            #~ test = []
            #~ test_date = []
            
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
            
            #~ print (row[0:12])
            #~ print (row[13:24])
            #~ print (row[25:36])
            #~ quit()
            
            #~ time_a.append(da)
            print('row   ',row)
            if float(row[21]) < 10000.:
                Bx.append(float(row[14]))
                By.append(float(row[15]))
                Bz.append(float(row[16]))
                Btot.append(float(row[13]))
                speed.append(float(row[21]))
                density.append(float(row[25]))


        #~ Bx = np.asarray(Bx)
        #~ By = np.asarray(By)
        #~ Bz = np.asarray(Bz)
        #~ Btot = np.asarray(Btot)
        #~ density = np.asarray(density)
        #~ speed = np.asarray(speed)

        #~ current_end_date = start_date +  dt.timedelta(minutes = cadence)

        #~ test_date.append(current_end_date)
        #~ test.append(float(fdatei[2]))


        if len(Bx) == 0 or len(density) == 0 or len(speed) == 0:
            print('No realtime solar wind data available... Aborting...')
            #~ return
            # maybe use Kp index here...

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
            
            #~ print (row[0:12])
            #~ print (row[13:24])
            #~ print (row[25:36])
            #~ quit()
            
            #~ time_a.append(da)
#           print (type(Bx))
#           quit()
                
            if float(row[21]) < 10000.:
                Bx.append(float(row[14]))
                By.append(float(row[15]))
                Bz.append(float(row[16]))
                Btot.append(float(row[13]))
                speed.append(float(row[21]))
                density.append(float(row[25]))

            
            #~ test_date.append(db)
            #~ test.append(float(fdatei[2]))

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
    print ('Speed  ', speed)
    lead_time = int(1.5e6/(60*speed[-1]))
    print ('Lead time (minutes)  ',lead_time/60.)
    
    forecast_time = time_sw + dt.timedelta(minutes = lead_time)
    
    sw_avg = { 'current_time' : current_time, 'time_latest_solar_wind' : last_date,  'forecast_time': forecast_time, 'Bx' : Bx_avg, 'By' : By_avg, 'Bz' : Bz_avg, 'B_average' : Btot_avg, 'v' : speed_avg, 'ni' : den_avg }
    
    
#   print (sw_avg)
#   quit()
    
    return sw_avg,  last_date, Bx, By,Bz, Btot, density, speed


    

                
                
