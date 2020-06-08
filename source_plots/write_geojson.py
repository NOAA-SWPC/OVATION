 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
import pandas as pd
import datetime as dt
import shutil

import geojson 


import geojson


# from geojson import Point, Feature, FeatureCollection, dump

# from df_to_geojson import df_to_geojson



def write_geojson( Gridded_Output_path, mdate, fdate, file_date, global_array):
	

	ofile_json = Gridded_Output_path+'Aurora_'+file_date+'.gjson'
	ofile_latest = Gridded_Output_path+'Aurora_latest.gjson'
	
	cols = ['lon','lat','aur']
		
	for ilon in range (360): 	
		for ilat in range(181):
 			if ilat == 0 and ilon == 0: data_array = np.array([(ilon-180),(ilat-90), (global_array[ilon,ilat])]) 
 			else: data_array = np.vstack((data_array,[(ilon-180),(ilat-90),(global_array[ilon,ilat])]))	

	data_array = data_array.reshape(65160,3)

# 		 These lines create a smaller geojsaon array for debugging
	
# 	for ilon in range (10):
# 	
# 		for ilat in range(5):
#  			if ilat == 0 and ilon == 0: data_array = np.array([(ilat-90),(ilon),(global_array[ilon,ilat])]) 
#  			else: data_array = np.vstack((data_array,[(ilat-90),(ilon),(global_array[ilon,ilat])]))				  

# 	data_array = data_array.reshape(50,3)
	 
	 
# 	Convert array into Panda array
	
	
	
	df = pd.DataFrame(data_array, columns=cols)
# 	df = df.truncate(after=50)
	
	points = []

	df.apply(lambda X: points.append( (int(X["lon"]), int(X["lat"]), int(X["aur"]))), axis=1)
	
	obs_time = mdate.strftime('%Y-%m-%dT%H:%M:%SZ')
	for_time = fdate.strftime('%Y-%m-%dT%H:%M:%SZ')
	
	with open(ofile_json, 'w+') as fp:

		fp.write('{\"Observation Time\": \"' + obs_time + '\", \"Forecast Time\": \"'+ for_time +'\", \"Data Format\": \"[Longitude, Latitude, Aurora]\", ')
		geojson.dump(geojson.MultiPoint(points), fp, sort_keys=True)
	
# 	Remove  sguiggly brackets
		fp.seek(134, os.SEEK_SET)
		fp.write(" ")
# 	
	
	fp.close()
	
	shutil.copyfile(ofile_json, ofile_latest)
	
	return()
