#  ***************************  Plots OVATION Model Output as Aurora Probability *****************

#  This code reads the output text file from the OVATION model and creates
#  northern and southern maps of the aurora with a polar perspective.
#  It reads in three differnt electron precipitation maps and combines them
#  into a single map of the total electron flux.   The fluxes are adjusted to
#  better represent the visual element of the aurora (i.e. descrete aurora is much easier
#  to see than diffuse aurora)

#  Input:  text file with...
#         - geomagnetic latitutde
#         - local time
#         - three electron (ignores the protons) energy depostion rates for
#         - Hemispheric Power

# Output: Two jpeg images (800x800 pixels) of the northern and southern polar regions
#
#  ******************************************************************************************

import os
# from os import path
import numpy as np
import pylab
# import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
# import cartopy.features as cpf

import matplotlib.pyplot as plt
import datetime as dt
# import math
import scipy.interpolate as sciint
from scipy.ndimage import gaussian_filter
import aacgmv2
from matplotlib.colors import LinearSegmentedColormap

from PIL import Image

# from output_geojson import output_geojson
# from write_HP_file import write_HP_file
from set_plot_colors_3_day_forecast import set_plot_colors_3_day_forecast
# from kp_and_g_scale import kp_and_g_scale



#  ***************** Set Run Parameters  *************************
#Home_path = os.environ.get('Home_path','./')
#Input_path = os.environ.get('Input_path','../Output')
#Output_path = os.environ.get('Output_path','../Output/NOWCAST')

Home_path = os.environ.get('home_path','./')
Input_path = os.environ.get('input_path','../output/FORECAST/model_output/')
Output_path = os.environ.get('output_path','../output/FORECAST/')
Logo_path = os.environ.get('logo_path','../logo/')

Image_Output_path = Output_path + 'images/'
# Gridded_Output_path = Output_path  + 'Text/Gridded/'
# HP_Output_path = Output_path + 'Text/Hemispheric_power/'

os.makedirs(Output_path + 'images/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'images/'))
os.makedirs(Output_path + 'images/north/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'images/north/'))
os.makedirs(Output_path + 'images/south/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'images/south/'))
# os.makedirs(Output_path + 'Text/Gridded/', exist_ok=True)
# print("making directory path (if necessary) {}".format(Output_path +'Text/Gridded/'))
# os.makedirs(Output_path + 'Text/Hemispheric_power/', exist_ok=True)
# print("making directory path (if necessary) {}".format(Output_path +'Text/Hemispheric_power/'))

# run_realtime = os.environ.get('run_realtime',True)   #If = True, then run once and use the standard "latest" file from the ovation model
# run_once = os.environ.get('run_once',False)  #If = True, then only plot one map or one pair of maps if plot_south = 1
plot_south = os.environ.get('plot_south',False)    #if = False, then north only, if = 1, then include southern hemispher plots
# output_global_ascii = os.environ.get('output_globa_ascii',True)   # Set this to True to get gridded ASCII output of global aurora distribution

lat_mult = 1.

# *  *************************************************************

	#  Relative weights of the three types of electron aurora (Biased toward discrete aurora)
mult_dif = .75
mult_mono = 1.25
mult_wave = 1.



	#  ****************	Define color map and color table for aurora**************

upper = 4.3		#Upper Limit for Aurora in plots	 in ergs/cm2/sec.   Setting this prevents the color map from wrapping


# ****************   Set the color scale for the aurora   ***************************

cdict1, fg_color, land_color, ocean_color, boarder_color = set_plot_colors_3_day_forecast()

os.system('rm ' + Output_path + 'images/north/*.jpg')
os.system('rm ' + Output_path + 'images/south/*.jpg')	
	
# NS = 1

for NS in range(2):

	lat_mult = 1.
	if NS == 1: lat_mult = -1
	
	file_list = []
	
	ipath = Input_path+"north/"
	
	if NS == 1: ipath = Input_path+"south/"
	
	for file in os.listdir(ipath):
		file_list = np.append(file_list, file)
		nloops = np.size(file_list)

	
	#******************   Create two plots for North and South hemispheres   **********
# 	
	nloops = 2 #  Set small for debugging
	
	for iloop in range(nloops):
		
		in_file = file_list[iloop]
		print ("Input File   ", in_file)
		
# 		out_file = in_file[0:24] + ".png"
		
		
		vlat = 80.
		vlon = -100
		if NS == 1:
			vlat = -90.
			vlon = -140
	
		start_hour = 0
	
		imin = 0
		ihour = 0
	
		input_file = open(ipath + in_file, 'r')		
	
	#  Read Input File and Parse Critical Variables
		
	# 	Parse Measurement Time
	
		fhead = input_file.readline()
		
		mdatei = input_file.readline() .strip()
		mdatei =mdatei.split()
		mdate1 = mdatei[2].split("-")
		mdate2 = mdatei[3].split(":")
	
		myear = int(mdate1[0])
		mmonth = (int(mdate1[1]))
		mday = int(mdate1[2])
		mhour = int(mdate2[0])
		mminute = int(mdate2[1])
		
	
		mdate = dt.datetime(myear, mmonth, mday, mhour, mminute)
		
	# 	Parse Forecast Time
	
		fdatei = input_file.readline() .strip()
		fdatei =fdatei.split()
		fdate1 = fdatei[2].split("-")
		fdate2 = fdatei[3].split(":")
	
		fyear = int(fdate1[0])
		fmonth = (int(fdate1[1]))
		fday = int(fdate1[2])
		fhour = int(fdate2[0])
		fminute = int(fdate2[1])
	
		fdate = dt.datetime(fyear, fmonth, fday, fhour, fminute)
		
	# 	Pares Hemispheric Power
	
		f_hp_head = input_file.readline().strip().split()
	
		HP = float(f_hp_head[2])			#Read Hemispheric Power
		
		f_hp_head = input_file.readline().strip().split()
		
		Kp = float(f_hp_head[2])			#Read Kp
	
		fhead = input_file.readline()
# 		fhead = input_file.readline()
		
		
	# 	**********************    Read Data from File   ****************************
	
	#  Set up lat-lon arrays and fill them from the input data
	
		numlon= 96
		numlat= 80
	
		glat = np.zeros(shape=(numlat,numlon+1))
		glon= np.zeros(shape=(numlat,numlon+1))
		alt=np.zeros(shape=(numlat,numlon+1))
		aur = np.zeros(shape=(numlat,numlon+1))
	
		lat = np.zeros(numlat)
		lon = np.zeros(numlat)
	
		itot = 0
		
		dec_time = float(fhour) + float(fminute)/60.
	
	   
	
		for ilon in range(0,numlon):
			for ilat in range(0,numlat):
				itot=itot+1
				row = input_file.readline().strip().split()
	
				aur[ilat,ilon] = (mult_dif * float(row[3]) + mult_mono * float(row[4]) + mult_wave * float(row[5]))
				lat[ilat] = float(row[1])
	
				#Convert local time into magnetic longitude
				lon[ilat] = (float(row[0])-dec_time)*15. + 71.
	
				#Convert Magnetic Coordinates to Geographic
				glat[(ilat,ilon)], glon[(ilat,ilon)], alt = aacgmv2.convert_latlon_arr(lat[ilat],lon[ilat], 500, fdate, method_code="A2G")
	
		input_file.close()
	
		glon[:,96] = glon[:,0]
		glat[:,96] = glat[:,0]
		aur[:,96] = aur[:,0]
	
		X= glat
		Y= glon
		Z= aur.ravel()
	
	
	
	# ***************	regrid to regular geographic coordinates **************
	
		glatlon = np.zeros(shape = (((numlat)*(numlon+1)),2))
	
		glatlon[:,0] = glon.ravel()
		glatlon[:,1] = glat.ravel()
	
		gx,gy = np.mgrid[-180:180:360j, 0:lat_mult*90:90j]			#Set even lat lon spacing
		geo_aur = sciint.griddata(glatlon,Z,(gx,gy),method='cubic')		#Regrid to even lat lon spacing
		geo_aur[359,:] = geo_aur[358,:]
		geo_aur[0,:] = geo_aur[1,:]
	
		geo_aur = np.nan_to_num(geo_aur)			#replace NANs with Zeros
	
		geo_aur = 1.1 * gaussian_filter(geo_aur,sigma = 2,mode = 'wrap')		#Smooth data to remove artifacts from geomag to geograph transform
	# 		geo_aur = 1.2 * gaussian_filter(geo_aur,sigma = 3,mode = 'wrap')		#Smooth data to remove artifacts from geomag to geograph transform
	
		np.clip(geo_aur, 0.,upper,out=geo_aur)  #clip the data at the upper plotting limit
		
	
		G = 0
		if Kp>=5: G=1
		if Kp>=6: G=2
		if Kp>=7: G=3
		if Kp>=8: G=4
		if Kp>=9: G=5
	
	# 	************   Set up aurora color map  **********************
	
		clevs = np.arange(.2,1.1*upper,0.25)
	
		auro_col = LinearSegmentedColormap('aurora_color',cdict1)
		plt.register_cmap(cmap=auro_col)
		cmap=plt.get_cmap('aurora_color')
	
	
	
	# ***************	Plot the aurora on a map of the globe**************

		F = pylab.gcf()		#Det screen DPI so that the final images are the correct size in pixels
		DPInch = F.get_dpi()
# 		DPInch =75
		font_mult = 75./DPInch

		fig = plt.figure(figsize=[800/DPInch,800/DPInch], facecolor = 'black' )#,dpi = DPInch)
		 
		ax1 = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(vlon,vlat))
	
		plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
		plt.margins(0,0)
 	
		ax1.gridlines()
		ax1.add_feature(Nightshade(fdate, color = 'black',delta=0.25, alpha = 0.4, zorder = 2))
 	
		ax1.add_feature(cfeature.OCEAN, facecolor = ocean_color, edgecolor = boarder_color )
		ax1.add_feature(cfeature.LAND, facecolor = land_color, edgecolor = boarder_color )
		ax1.add_feature(cfeature.BORDERS, edgecolor = boarder_color )
		ax1.add_feature(cfeature.STATES,  edgecolor = boarder_color)
 	
		geo_aur = np.transpose(geo_aur)

		lat_min = 1
		lat_max = 90
		if NS == 1:
 			lat_min = -1
 			lat_max = -90
 	
		im1 = ax1.imshow(geo_aur, vmin= 0, vmax = 5, transform = ccrs.PlateCarree(),
			 extent = (-180, 180, lat_min, lat_max), origin = (0,0), cmap = cmap, zorder = 2 )
 	

	# 		  *********************  Lable plot  ***************************
	
		text_color = (1.,1.,.4)
		background_color = (0,0,0,.5)
		
		now = dt.datetime.utcnow()
		lt = (1/60.)*(fdate - mdate).seconds       #Local Time
 	
		#Upper Right Lable

		ur_text = '                           NOAA Space Weather Prediction Center\n\n\nFor %s (UTC)' % fdate.strftime("%Y-%m-%d %H:%M")
		plt.figtext(0.89,0.91,ur_text, color = text_color, backgroundcolor = background_color, size =16,ha = 'right')	
		plt.figtext(0.89,0.94,'Multi-Day Aurora Forecast', color = text_color, size = 30,weight = 'bold',ha = 'right')

		ur_text2 = 'Forecast Kp: %1i (Range 0 to 9)\nForecast G-Scale: %1i (Range 0 to 5)' %(int(Kp), int(G))
		plt.figtext(0.99,.85, ur_text2, color=text_color, backgroundcolor = background_color, size = 14,weight = 'bold',ha = 'right')

 	
		# ~ Lower Right
		
		lr_text = '\nModel Run at %s (UTC)\nKp Forecast Made at %s (UTC)'%(now.strftime("%Y-%m-%d %H:%M"), mdate.strftime("%Y-%m-%d %H:%M"))
		plt.figtext(0.99,0.02, lr_text,color=text_color,backgroundcolor = background_color, size = 12,ha = 'right')
		plt.figtext(0.99,0.06,'OVATION Aurora Model', color = text_color, size =18, ha = 'right')
# 		plt.figtext(0.99,0.08, 'Model Run at %s (UTC)' % now.strftime("%Y-%m-%d %H:%M"),color=text_color,size = 12,ha = 'right')
		
 	
 	
 	
		#  **********************   Add and Label Colorbar.  ***********************
		plt.figtext(0.2,.055,'Probability of Aurora\n', color =text_color, backgroundcolor = background_color, size =16, ha = "center",weight = 'bold')
		plt.figtext(0.1,.055,'Low                     High', color = text_color, size =14)
		
		cax = plt.axes([0.1,0.01,0.20,0.027])
		cb=plt.colorbar(im1,cax=cax,orientation = "horizontal", ticks= [1,2,3,4,5])
 	
		cb.patch.set_facecolor('darkgray')
 	
 	# 		cb.set_label('Energy Deposition  (erg/cm2)', color=fg_color,size =12)
		cb.ax.xaxis.set_tick_params(color=fg_color)
		cb.ax.set_xticklabels([1,2,3,4],color=fg_color)
		cb.outline.set_edgecolor(fg_color)

	

	#  *************************   Save Graphics and Text *****************************
		
	
	
	#    ****   Create date-time string for file names    ****
	# 	Note that the filenames will have the forecast date-time not the actual observation time
	
		syear = str("%4.4i"%fyear)
		smonth = str("%2.2i"%fmonth)
		sday = str("%2.2i"%fday)
		shour = str("%2.2i"%fhour)
		sminute = str("%2.2i"%fminute)
	
		file_date = 'Text_' + syear + '-' + smonth + '-' + sday + '_' + shour + sminute
		
				
				
	#
	# ***********************   Add NOAA Logo  to image  **************************			
	
		icon = Image.open(Logo_path + "NOAA_logo4.png")
		
		height = icon.size[1]
		icon = np.array(icon).astype(np.float) / 255

		newax = fig.add_axes([0.9, 0.9, 0.1, 0.1], anchor='NE', zorder=-1)
		newax.imshow(icon)
		newax.axis('off')

	#  ******************   Save North and South Image Files   *******************************	
		

		ofl = 'north/aurora_forecast_N_' + file_date + '.jpg'
		if NS == 1:
			ofl = 'south/aurora_forecast_S_' + file_date + '.jpg'

	
	#   *******************************Save Final Image to "Latest Image" File and Final file****************************************
	
		
	
		ofile_image = Image_Output_path+ofl
		plt.savefig(ofile_image ,facecolor = 'black',dpi=DPInch)#,bbox_inches = 'tight', pad_inches = 0)#,dpi=DPInch)
		
		plt.close()

