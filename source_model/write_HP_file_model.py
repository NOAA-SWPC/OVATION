#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

write_HP_file 
appends the current HP file with the latest HP values
Creates a new HP file at the beginning of each day.


Created on Wed Mar 18 08:45:52 2020

@author: rodney.viereck
"""



from os import path
import datetime as dt

def write_HP_file(HP, HP_Output_path, Header_path, time_lab, time_for, NS, NorthSouth, mode):
	
	syear = time_lab.year
	smonth = time_lab.month
	sday = time_lab.day
	shour = time_lab.hour
	sminute = time_lab.minute
	
	fyear = time_for.year
	fmonth = time_for.month
	fday = time_for.day
	fhour = time_for.hour
	fminute = time_for.minute
	   
	Date_string = "%4i%2.2i%2.2i"%(syear,smonth,sday)
	if mode == "HISTORIC": Date_string = "%4i%2.2i"%(syear,smonth)   #Output Monthly files instead of Daily
	Date_string2 = "%4i-%2.2i-%2.2i"%(syear,smonth,sday)
	Date_Time_string = "%4i-%2.2i-%2.2i_%2.2i:%2.2i"%(syear,smonth,sday,shour,sminute)
	FDate_Time_string = "%4i-%2.2i-%2.2i_%2.2i:%2.2i"%(fyear,fmonth,fday,fhour,fminute)
	
	HP_str = "    %4i"%(HP)
	HP_file = HP_Output_path + "SWPC_aurora_power_" + Date_string + ".txt"	

	
	if path.exists(HP_file) == False:
		with open(Header_path + 'HP_Header_1.txt') as h1:
			hd1 = h1.read()
		with open(Header_path + 'HP_Header_2.txt') as h2:
			hd2 = h2.read()
		hd1 = hd1[0:len(hd1)-1]
		hd1 += Date_string2 + "\n"
		hd1 += hd2
		with open(HP_file, 'w') as data_out:
			data_out.write(hd1)
		data_out.close()

	data_out = open(HP_file, 'a') 	
	if NS == 0: 
		if NorthSouth == True: 
			data_out.write(Date_Time_string + "    " + FDate_Time_string + HP_str)
		else:
			data_out.write(Date_Time_string + "    " + FDate_Time_string + HP_str + "\n")
	if NS != 0: data_out.write(HP_str +"\n")
	
	data_out.close()
	
	return()
