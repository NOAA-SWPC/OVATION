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
import numpy as np

# from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
# import cartopy.features as cpf

import matplotlib.pyplot as plt
import datetime as dt
import math
import scipy.interpolate as sciint
from scipy.ndimage import gaussian_filter
import aacgmv2
from matplotlib.colors import LinearSegmentedColormap

from PIL import Image

# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
# from matplotlib.figure import Figure


#  ***************** Set Run Parameters  *************************

Home_path = os.environ.get('Home_path','./')
Input_path = os.environ.get('Input_path','../Output')
Output_path = os.environ.get('Output_path','../Output/NOWCAST')
Logo_path = os.environ.get('Home_path') + 'logo/'
Image_Output_path = os.environ.get('Output_path','Final_Output/Images/')
Text_Output_path = os.environ.get('Output_path','Final_Output/Text/')

run_realtime = os.environ.get('run_realtime',True)   #If = True, then run once and use the standard "latest" file from the ovation model
run_once = os.environ.get('run_once',False)  #If = True, then only plot one map or one pair of maps if plot_south = 1
plot_south = os.environ.get('plot_south',False)    #if = False, then north only, if = 1, then include southern hemispher plots
output_global_ascii = os.environ.get('output_globa_ascii',True)   # Set this to True to get gridded ASCII output of global aurora distribution

#  If not running realtime, set the starttime, stoptime and cadence.
#         Note, the filenames will be generated based on these parameters and filenames must match existing files.
# start_date = os.environ.get('start_date',dt.datetime(2017,9,28,17,00))   	#YYYY, MM, DD, HH, MM
# end_date =  os.environ.get('end_date',dt.datetime(2017,9,29,0,0))   	#YYYY, MM, DD, HH, MM
# cadence = os.environ.get('cadence',60)		#minutes
start_date = os.environ.get('start_date', '09-27-2017 00:00')
start_date = dt.datetime.strptime(start_date, '%m-%d-%Y %H:%M')
end_date = os.environ.get('end_date', '09-30-2017 00:00')
end_date = dt.datetime.strptime(end_date, '%m-%d-%Y %H:%M')
cadence = os.environ.get('cadence', 60)

# *  *************************************************************

if run_realtime == True:
	run_once = True
	plot_south = True
	Input_path = Input_path + '/NOWCAST/Text/'
else:
	Input_path = Input_path + '/Historic/Text/'

ins = 1
if plot_south == True:
	ins = 2


	# ******************	 Variable Definition  **************

if run_once == True:end_date = start_date

run_date = start_date

imin = 0
ihour = 0


#  Relative weights of the three types of electron aurora (Biased toward discrete aurora)
mult_dif = .75
mult_mono = 1.25
mult_wave = 1.

#  Parameters used to Calculate Kp from Hp
y0 = 6.9498
A1 = 2.04984
t1 = 1.8400


	#  ****************	Define color map and color table for aurora**************

upper = 4.3		#Upper Limit for Aurora in plots	 in ergs/cm2/sec.   Setting this prevents the color map from wrapping

fg_color = 'yellow'	  #Set forground color

# ****************   Set the color scale for the aurora   ***************************
#     Note, the Alpha term is for transparancy

cdict1 = {'red':  ((0.0, 0.0, 0.0),
				   (0.1, 0.0, 0.0),
				   (0.3, 0.2, 0.2),
				   (0.5, 1.0, 1.0),
				   (0.8, 1.0, 1.0),
				   (0.9, 1.0, 1.0),
			       (1.0, 0.8, 0.8)),

		 'green': ((0.0, 0.5, 0.5),
				   (0.1, 0.9, 0.9),
				   (0.3, 1., 1.),
				   (0.5, 1., 1.),
				   (0.8, .5, .5),
				   (0.9, 0.0, 0.0),
			       (1.0, 0.0, 0.0)),

		 'blue':  ((0.0, 0.0, 0.0),
					(0.1, 0.0, 0.0),
				   (0.3, 0.0, 0.0),
				   (0.5, 0.0, 0.0),
				   (0.8, 0.0, 0.0),
				   (0.9, 0.0, 0.0),
				   (1.0, 0.0, 0.0)),

		'alpha': ((0.0, 0.0, .0),
 				   (0.1, .9, .9),
 				   (0.3, 1.0, 1.0),
 				   (0.5, 1.0, 1.0),
 				   (0.8, 1.0, 1.0),
 				   (0.9, 1.0, 1.0),
 				   (1.0, 1.0, 1.0))
 			}

# *****************************   Start of Loop  *****************************
#       If in real-time mode or run-once mode, the code will run staight through and not loop

# ins = 1   #Remove before launch

while run_date <= end_date:
	for NS in range(ins):
		lat_mult = 1
		if NS == 1: lat_mult = -1
		smin = str("%2.2i%2.2i"%(ihour,imin))

		# 	ifl = 'aurora_N1_2017-09-06_'+smin+'.txt'

		year = str("%4.4i"%run_date.year)
		month = str("%2.2i"%run_date.month)
		day = str("%2.2i"%run_date.day)
		hour = str("%2.2i"%run_date.hour)
		minute = str("%2.2i"%run_date.minute)

# 	 Generate input/output file names

		file_date = 'Text_' + year + '-' + month + '-' + day + '_' + hour + minute

		ifl = 'North/aurora_N_' + year + '-' + month + '-' + day + '_' + hour + minute
		if NS == 1:
			ifl = 'South/aurora_S_' + year + '-' + month + '-' + day + '_' + hour + minute

		ofl = ifl

		if run_realtime == True:
			ifl = 	'Ovation_Latest_aurora_N'
			if NS == 1: ifl = 'Ovation_Latest_aurora_S'

		vlat = 80.
		vlon = -100
		if NS == 1:
			vlat = -90.
			vlon = -90

		start_hour = 0

		imin = 0
		ihour = 0


		input_file = open(Input_path + ifl + '.txt', 'r')

#  Read Input File and Parse Critical Variables

		fhead = input_file.readline()

		fdatei = input_file.readline() .strip()

		fdatei =fdatei.split()
# 		print ('data date  ',fdatei)
		fdate1 = fdatei[2].split("-")
		fdate2 = fdatei[3].split(":")

		year = int(fdate1[0])
		month = (int(fdate1[1]))
		day = int(fdate1[2])

		hour = int(fdate2[0])
		min = int(fdate2[1])
		dec_time = float(hour) + float(min)/60.

		mdate = dt.datetime(year, month, day, hour, min)

		fdatei = input_file.readline() .strip()
		fdatei =fdatei.split()
		fdate1 = fdatei[2].split("-")
		fdate2 = fdatei[3].split(":")

		year = int(fdate1[0])
		month = (int(fdate1[1]))
		day = int(fdate1[2])

		hour = int(fdate2[0])
		min = int(fdate2[1])

		fdate = dt.datetime(year, month, day, hour, min)

		f_hp_head = input_file.readline().strip().split()

		HP = float(f_hp_head[2])			#Read Hemispheric Power

		fhead = input_file.readline()
		fhead = input_file.readline()

#  Set up lat-lon arrays and fill them

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

		np.clip(geo_aur, 0.,upper,out=geo_aur)





#  Calculate Kp and G scales from Hemispheric Power


		lnval = ((HP/A1)-y0)
		if lnval <= 0.0: lnval = 1e-5
		Kp= t1*math.log(lnval)
		if Kp < 0.0: Kp = 0.

# Create G Scale from Kp

		G = 0
		if Kp>=5: G=1
		if Kp>=6: G=2
		if Kp>=7: G=3
		if Kp>=8: G=4
		if Kp>=9: G=5



		now = dt.datetime.now()

	# 	************   Set up aurora color map  **********************

		clevs = np.arange(.2,1.1*upper,0.25)

		auro_col = LinearSegmentedColormap('aurora_color',cdict1)
		plt.register_cmap(cmap=auro_col)
		cmap=plt.get_cmap('aurora_color')



	# ***************	Plot the aurora on a map of the globe**************



		landcolor = (.5,.5,.1)
		oceancolor = (0.5, 0.5, 0.9)
		boardercolor = (.7,.7,.7)

		fig = plt.figure(figsize=[10, 10],facecolor = 'black' )

		ax1 = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(vlon,vlat))

		ax1.gridlines()
		ax1.add_feature(Nightshade(fdate, color = 'black',delta=0.25, alpha = 0.4, zorder = 2))

		ax1.add_feature(cfeature.OCEAN, facecolor = oceancolor, edgecolor = boardercolor )
		ax1.add_feature(cfeature.LAND, facecolor = landcolor, edgecolor = boardercolor )
		ax1.add_feature(cfeature.BORDERS, edgecolor = boardercolor )
		ax1.add_feature(cfeature.STATES,  edgecolor = boardercolor)


		geo_aur = np.transpose(geo_aur)

		lat_min = 1
		lat_max = 90
		if NS == 1:
			lat_min = -1
			lat_max = -90

		im1 = ax1.imshow(geo_aur, vmin= 0, vmax = 5, transform = ccrs.PlateCarree(),
			 extent = (-180, 180, lat_min, lat_max), origin = (0,0), cmap = cmap, zorder = 2 )


# 		  *********************  Lable plot  ***************************

		lt = (1/60.)*(fdate - mdate).seconds       #Local Time

		#Upper Left

		plt.figtext(0.11,0.976,'NOAA Space Weather Prediction Center', color = 'yellow', size =16)
		plt.figtext(0.11,0.94,'Aurora Forecast', color = 'yellow', size = 30,weight = 'bold')
		plt.figtext(0.11,0.915, 'For %s (UTC)' % fdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 16)

		#~ Upper Right

		plt.figtext(0.99,0.92,'Hemispheric Power: %5.1f GW (Range 5 to 200)' %HP,color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.94,'Estimated Kp: %1i (Range 0 to 9)' %int(Kp),color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.96,'Estimated G-Scale: %1i (Range 0 to 5)' %int(G),color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.98,'Forecast Lead Time: %3i minutes' %int(lt), color = 'yellow', size =12, ha = 'right')

		# ~ Lower Right

		plt.figtext(0.99,0.08,'OVATION Aurora Model', color = 'yellow', size =18, ha = 'right')
		plt.figtext(0.99,0.06, 'Model Run at %s (UTC)' % now.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.04, 'L1 Observations at %s (UTC)' % mdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 12,ha = 'right')

# 		#  ************************   Getting rid of black boarder   **********************
		plt.gca().set_axis_off()
		plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
		            hspace = 0, wspace = 0)
		plt.margins(0,0)
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		plt.gca().yaxis.set_major_locator(plt.NullLocator())


		#  **********************   Add and Label Colorbar.  ***********************

		cax = plt.axes([0.05,0.08,0.30,0.027])
		cb=plt.colorbar(im1,cax=cax,orientation = "horizontal", ticks= [1,2,3,4,5])

		cb.patch.set_facecolor('darkgray')

# 		cb.set_label('Energy Deposition  (erg/cm2)', color=fg_color,size =12)
		cb.ax.xaxis.set_tick_params(color=fg_color)
		cb.ax.set_xticklabels([1,2,3,4],color=fg_color)
		cb.outline.set_edgecolor(fg_color)

		cb.set_label('Approximate Energy Deposition', color=fg_color)

		plt.figtext(0.18,.02,'(erg/cm2)', color=fg_color,ha = "center",size =12)
		plt.figtext(0.32,.05,'>4', color=fg_color,size =12)
		plt.figtext(0.05,.12,'10%        50%       90%', color = 'yellow', size =14)
		plt.figtext(0.2,.14,'Probability of Aurora', color = 'yellow', size =16, ha = "center",weight = 'bold')


#
		# ***********************    Save in a temp file and then read it back in as an image   **************************
#

		ofile_image = Image_Output_path+"temp1.png"

		plt.savefig(ofile_image ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0)


		im = Image.open(ofile_image)
#
		# ***********************   Add NOAA Logo  to image  **************************


		icon = Image.open(Logo_path + "NOAA_logo4.png")
		icon_size = (70,70)

		icon2 = icon.resize(icon_size)
		xs,ys = icon2.size

#
# # 		#NOAA logo position
		x0 = 0
		y0 = 7
		x1 = x0+xs
		y1 = y0+ys

		im.paste(icon2, (x0,y0,x1,y1), icon2)




#  *************************   Save Graphics and Text *****************************

#    ******************   Output Text Data to ASCII Array  *********************
		if output_global_ascii == True and plot_south == True:
			if NS ==0:
				global_array = np.zeros([360,181])
				global_array[0:360,91:181] = np.transpose(geo_aur)
			if NS ==1:
				gaf = np.flip(np.transpose(geo_aur), 1)
				global_array[0:360,0:90] = gaf
				global_array = np.flip(global_array, 0)   #(np.transpose(global_array),0)
				global_array = global_array *100.


			ofile_text = Text_Output_path+'OVATION_Text_Latest.txt'
			np.savetxt(ofile_text, global_array, delimiter=' ', fmt='%5.1f')
			ofile_text = Text_Output_path+'OVATION_'+file_date+'.txt'
			np.savetxt(ofile_text, global_array, delimiter=' ', fmt='%5.1f')

#   ***********************************************************************

		#Save Final Image to "Latest Image" File
		if run_realtime == True:
			ofile_image = Image_Output_path+ifl+'.png'
			im.save(ofile_image ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0)


		#Save Final Image to File
		ofile_image = Image_Output_path+ofl+'.png'
		plt.savefig(ofile_image ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0)


# 	Increment run_date and return to top of loop
	run_date = run_date + dt.timedelta(minutes = cadence)


	# im2.show()








# 	run_date = run_date + dt.timedelta(minutes = cadence)


	# im2.show()
