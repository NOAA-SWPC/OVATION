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

# import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
# import cartopy.features as cpf

import matplotlib.pyplot as plt

import pylab
import datetime as dt

import scipy.interpolate as sciint
from scipy.ndimage import gaussian_filter
import aacgmv2
from matplotlib.colors import LinearSegmentedColormap

from PIL import Image

from write_geojson import write_geojson
from set_plot_colors import set_plot_colors
from kp_and_g_scale import kp_and_g_scale



#  ***************** Set Run Parameters  *************************

#Set mode to be NOWCAST or HISTORIC
#Note:  FORECAST mode plots are made by a seperate program

mode = 'NOWCAST'

Home_path = os.environ.get('Home_path','./')
Input_path = os.environ.get('Input_path','../output')
Output_path = os.environ.get('Output_path','../output/'+mode+'/ovation_products/')
Logo_path = os.environ.get('Logo_path','../logo/')

Image_Output_path = Output_path + 'images/'
Gridded_Output_path = Output_path  + 'gridded_text/'
# HP_Output_path = Output_path + 'Text/Hemispheric_power/'

os.makedirs(Output_path + 'images/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'images/'))
os.makedirs(Output_path + 'images/north/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'images/north/'))
os.makedirs(Output_path + 'images/south/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'images/soutth/'))
os.makedirs(Output_path + 'gridded_text/', exist_ok=True)
print("making directory path (if necessary) {}".format(Output_path +'gridded_text/'))

run_realtime = os.environ.get('run_realtime',True)   #If = True, then run once and use the standard "latest" file from the ovation model
run_once = os.environ.get('run_once',False)  #If = True, then only plot one map or one pair of maps if plot_south = 1
plot_south = os.environ.get('plot_south',True)    #if = False, then north only, if = 1, then include southern hemispher plots
output_global_ascii = os.environ.get('output_globa_ascii',True)   # Set this to True to get gridded ASCII output of global aurora distribution
subset = os.environ.get('subset',True)	#If subset = True and mode = 'HISTORIC' then priocess only a suybset of the data
#  ******************************   To look at a subset of the output from the model  *********************

if mode == 'HISTORIC' and subset == True:
	start_date = os.environ.get('start_date', '09-28-2017 00:00')
	start_date = dt.datetime.strptime(start_date, '%m-%d-%Y %H:%M')
	end_date = os.environ.get('end_date', '09-29-2017 00:00')
	end_date = dt.datetime.strptime(end_date, '%m-%d-%Y %H:%M')
	cadence = os.environ.get('cadence', 60)   #Cadence in Minutes

	print ("Processing only a subset between ",start_date," and ",end_date)
# *  *************************************************************

if mode == 'NOWCAST':
	run_once = True
	plot_south = True
	create_json = True
	nloops = 1
	Input_path = Input_path + '/NOWCAST/model_output/'
	
if mode == 'HISTORIC':
	run_once = False
	plot_south = False
	create_json = False
	Input_path = Input_path + '/HISTORIC/model_output/'

ins = 1
if plot_south == True:
	ins = 2


#  Relative weights of the three types of electron aurora (Biased toward discrete aurora)
mult_dif = .75
mult_mono = 1.25
mult_wave = 1.

	#  ****************	Define color map and color table for aurora**************

upper_limit = 4.3		#Upper Limit for Aurora in plots	 in ergs/cm2/sec.   Setting this prevents the color map from wrapping


# ****************   Set the color scale for the aurora   ***************************

cdict1, fg_color, land_color, ocean_color, boarder_color = set_plot_colors()

#******************   Create two plots for North and South hemispheres   **********

# NS = 1   #One plot only for debuggingh

global_array = np.zeros([360,181])   #Set up global array of aurora for sending to the geojson file.

for NS in range(ins):
	lat_mult = 1
	if NS == 1: lat_mult = -1



	ifl = 	'Ovation_Latest_aurora_N'
	if NS == 1: ifl = 'Ovation_Latest_aurora_S'

	vlat = 80.
	vlon = -100.
	if NS == 1:
		vlat = -90.
		vlon = -120.

	start_hour = 0

	imin = 0
	ihour = 0
	
	
	if mode == 'HISTORIC':
		
		file_list = []
		
		ipath = Input_path + 'north/'
		if NS == 1: ipath = Input_path + 'north/'
		
		for file in os.listdir(ipath):
			file_list = np.append(file_list, file)
			nloops = np.size(file_list)
			
# 			nloops = 1		# debug mode

	
	for iloop in range(nloops):
		

		if mode == 'NOWCAST':
			in_file = Input_path + ifl + '_.txt'
			input_file = open(in_file, 'r')				
			print ("Input File   ", in_file)	
            
		else:		
			in_file = file_list[iloop]
			input_file = open(ipath + in_file, 'r')				
			print ("Input File   ", ipath+in_file)		
		
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
		
		dec_time = float(mhour) + float(mminute)/60.
	
		mdate = dt.datetime(myear, mmonth, mday, mhour, mminute)
		
# 		if mdate >= start_date and mdate <= end_date: 
	# 	lt_date = dt.datetime(myear, mmonth, mday,mhour)
		
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
	
		fhead = input_file.readline()
		fhead = input_file.readline()
		
		
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
	
	
		for ilon in range(0,numlon):
			for ilat in range(0,numlat):
				itot=itot+1
				row = input_file.readline().strip().split()
	
				aur[ilat,ilon] = (mult_dif * float(row[3]) + mult_mono * float(row[4]) + mult_wave * float(row[5]))
				
				lat[ilat] = float(row[1])
	
				#Convert local time into magnetic longitude
	
				lt_hour = float(row[0])			
				lon[ilat] = aacgmv2.convert_mlt(lt_hour, mdate, m2a=True)
				
				#Convert Magnetic Coordinates to Geographic
				glat[(ilat,ilon)], glon[(ilat,ilon)], alt = aacgmv2.convert_latlon_arr(lat[ilat],lon[ilat], 500, fdate, method_code="A2G")
	
# 		input_file.close()
		
	
		glon[:,96] = glon[:,0]
		glat[:,96] = glat[:,0]
		aur[:,96] = aur[:,0]
	
		Z= aur.ravel()
	
	
	
	# ***************	regrid to regular geographic coordinates **************
	
		glatlon = np.zeros(shape = (((numlat)*(numlon+1)),2))
	
		glatlon[:,0] = glon.ravel()
		glatlon[:,1] = glat.ravel()
	
		gx,gy = np.mgrid[-180:180:361j, 0:lat_mult*90:90j]			#Set even lat lon spacing
		geo_aur = sciint.griddata(glatlon,Z,(gx,gy),method='cubic')		#Regrid to even lat lon spacing
	
		geo_aur[360,:] = geo_aur[359,:]
		geo_aur[0,:] = geo_aur[1,:]
		geo_aur[:,89] = geo_aur[:,88]
		geo_aur[:,0] = geo_aur[:,1]	
		
	# 	************************   Clean up the model results   *********************************8888888
	
		geo_aur = np.nan_to_num(geo_aur)			#replace NANs with Zeros
		geo_aur = 1.1 * gaussian_filter(geo_aur,sigma = 2,mode = 'wrap')		#Smooth data to remove artifacts from geomag to geograph transform
		np.clip(geo_aur, 0.,upper_limit,out=geo_aur)  #clip the data at the upper plotting limit
		
		for ilat in range(90):
			for ilon in range (361):
				iilon = int(gx[ilon,ilat])
				iilat = 90 + int(gy[ilon,ilat])
				
				global_array[iilon,iilat] = geo_aur[ilon,ilat]
	
	#  *****************   Call subroutine to Calculate Kp and G scales from Hemispheric Power  ********************
		
		Kp, G = kp_and_g_scale(HP)
	
	
	# 	************   Set up aurora color map  **********************
	
		clevs = np.arange(.2,1.1*upper_limit,0.25)
	
		auro_col = LinearSegmentedColormap('aurora_color',cdict1)
		plt.register_cmap(cmap=auro_col)
		cmap=plt.get_cmap('aurora_color')
	
	
	# ***************	Plot the aurora on a map of the globe**************
		
		F = pylab.gcf()		#Det screen DPI so that the final images are the correct size in pixels
	# 	DPI = F.get_dpi()
		DPInch =75
		font_mult = 75./DPInch
	
	
		fig = plt.figure(figsize=[800/DPInch, 800/DPInch], facecolor = 'black' ,dpi = DPInch)
	
		ax1 = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(vlon,vlat))
		
	#	ax1 = plt.subplot(1, 1, 1, projection=ccrs.Miller(central_longitude=0, globe=None))
	#	ax1 = plt.subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0, globe=None, false_easting=None, false_northing=None))
	
		ax1.gridlines()
		ax1.add_feature(Nightshade(fdate, color = 'black',delta=0.25, alpha = 0.6, zorder = 2))   #Add Day-night terminator
	
		ax1.add_feature(cfeature.OCEAN, facecolor = ocean_color, edgecolor = boarder_color  )
		ax1.add_feature(cfeature.LAND, facecolor = land_color)
		ax1.add_feature(cfeature.BORDERS, edgecolor = boarder_color, zorder = 3 )
		ax1.add_feature(cfeature.STATES,  edgecolor = boarder_color, zorder = 3 )
		ax1.add_feature(cfeature.LAND, facecolor = (0,0,0,0),edgecolor = boarder_color, zorder = 3)
	
		geo_aur = np.transpose(geo_aur)
	
		lat_min = 1
		lat_max = 90
		if NS == 1:
			lat_min = -1
			lat_max = -90
	
		im1 = ax1.imshow(geo_aur, vmin= 0, vmax = 5, transform = ccrs.PlateCarree(),
			 extent = (-180,180, lat_min, lat_max), origin = (0,0), cmap = cmap, zorder = 2 )
	
	
	# 		  *********************  Lable plot  ***************************
	
		now = dt.datetime.now()
		lt = (1/60.)*(fdate - mdate).seconds       #Local Time
		
		text_color = 'yellow'
		background_color = (0,0,0,.5)
		
		font_30 = (30. * font_mult)
		font_16 = (16. * font_mult)
		font_14 = (14. * font_mult)
		font_12 = (12. * font_mult)
	
		#Upper Left lable
	
		ul_text = 'NOAA Space Weather Prediction Center\n\n\nFor %s (UTC)' % fdate.strftime("%Y-%m-%d %H:%M")
		plt.figtext(0.11,0.91,ul_text, color = text_color, backgroundcolor = background_color, size =font_16)
		
		plt.figtext(0.11,0.935,'Aurora Forecast', color = text_color, size = font_30,weight = 'bold')
	
	
		#~ Upper Right Lable
	
		
		
		ur_text = 'Forecast Lead Time: %3i minutes\n\
Estimated Kp: %1i (Range 0 to 9)\n\
Estimated G-Scale: %1i (Range 0 to 5)\n\
HPI: %5.1f GW (Range 5 to 200)'\
			 %(int(lt),int(Kp), int(G),  HP)
	
		plt.figtext(0.99,0.91, ur_text,color = text_color, backgroundcolor=background_color,size = font_16, ha = 'right')
	
		# ~ Lower Right Lable
			 
		lr_text = 'OVATION Aurora Model\nModel Run at %s (UTC)\nL1 Observations at %s (UTC)'\
			  %(now.strftime("%Y-%m-%d %H:%M"), mdate.strftime("%Y-%m-%d %H:%M"))
	
		plt.figtext(0.99,0.01, lr_text,color = text_color, backgroundcolor = background_color,size = font_14,ha = 'right')
	
	
		#  **********************   Add and Lable Colorbar.  ***********************
	
		cax = plt.axes([0.05,0.08,0.30,0.027])
		cb=plt.colorbar(im1,cax=cax,orientation = "horizontal", ticks= [1,2,3,4,5])
	
		cb.patch.set_facecolor('darkgray')
		cb.set_label('Approximate Energy Deposition', color=fg_color)

		cb.ax.set_xticklabels(['1','2','3','4','5'],color=fg_color)
		cb.outline.set_edgecolor(fg_color) 	
	
		plt.figtext(0.2,.12,'   Probability of Aurora   \n', color = 'yellow',backgroundcolor = background_color, size =font_16, ha = "center",weight = 'bold')
		plt.figtext(0.05,.12,'10%            50%            90%', color = 'yellow', size =font_14) 
		cb_text = '       1        2        3        4        >4       \n(erg/cm2)\nApproximate Energy Deposition'
		plt.figtext(0.2,.02, cb_text, color=fg_color,ha = "center",size =font_12,backgroundcolor = background_color)

	
	#  ************************   Getting rid of extra black boarder   **********************
	# 	Without these lines, the output image has excessive black space around it
	
		plt.gca().set_axis_off()
		plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
		plt.margins(0,0)
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		plt.gca().yaxis.set_major_locator(plt.NullLocator())
		

	
	#  *************************   Save Graphics and Text *****************************
		
	
	
	#    ****   Create date-time string for file names    ****
	# 	Note that the filenames will have the observation date-time not the actual foreast time
	
		syear = str("%4.4i"%myear)
		smonth = str("%2.2i"%mmonth)
		sday = str("%2.2i"%mday)
		shour = str("%2.2i"%mhour)
		sminute = str("%2.2i"%mminute)
	
		file_date = syear + '-' + smonth + '-' + sday + '_' + shour + sminute
		
	
	#
	# ***********************   Add NOAA Logo  to image  **************************		
	#
	# ***********************    Save in a temp file and then read it back in as an image   **************************
	# 			This is the only way I could figure out how to add the NOAA logo (a jpeg file) to the aurora map
	#

		
		ofile_image = Image_Output_path+"temp1.png"
		plt.savefig(ofile_image ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0,dpi=DPInch)
		plt.clf()
		im = Image.open(ofile_image)	
		os.remove(ofile_image)
	
		icon = Image.open(Logo_path + "NOAA_logo4.png")
		icon_size = (70,70)
	
		icon2 = icon.resize(icon_size)
		xs,ys = icon2.size
	
	
	# 		#NOAA logo position
		x0 = 0
		y0 = 7
		x1 = x0+xs
		y1 = y0+ys
	
		im.paste(icon2, (x0,y0,x1,y1), icon2)
		im = im.convert('RGB')  #Converting from RGBA to RGB so that it can be saved as a jpg
	
		
	#  ******************   Save North and South Image Files   *******************************	
		
		image_type = '.jpg'
		
		ofl_latest = 'latest_aurora_N'+image_type
		ofl = 'north/aurora_N_' + file_date + image_type
		if NS == 1:
			ofl_latest = 'latest_aurora_S'+image_type
			ofl = 'south/aurora_S_' + file_date + image_type		
	
	#   *******************************Save Final Image to "Latest Image" File and Final file****************************************
	
		ofile_image = Image_Output_path+ofl
		im.save(ofile_image ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0) 
		if mode == 'NOWCAST':
			ofile_image_latest = Image_Output_path+ofl_latest
			im.save(ofile_image_latest ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0)
		
		im.close()
	#   **************************   Call output_geojson subroutine    *************************
		
	if create_json: 
		global_array = global_array * 100./upper_limit    #Scale data to go from 0 to 100
		np.clip(global_array, 0.,100.,out=global_array)  #clip the data at 100
					
		write_geojson(Gridded_Output_path, mdate, fdate,file_date, global_array)	

	input_file.close()