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


import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import datetime as dt
import math 
import scipy.interpolate as sciint
from scipy.ndimage import gaussian_filter
import aacgmv2
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image


#  ***************** Run Parameters  *************************
    

Input_path = '/home/rodney.viereck/python/ovation/Output/Historic/Text/'
Output_path = '/home/rodney.viereck/python/ovation/Final_Output/Final_Images/'

run_once =  1   #If = 1, then only plot one map

plot_south = 0    #if = 0, then north only, if = 1, then include southern hemispher plots
ins = 1 + plot_south

run_realtime = 1   #If = 1, then run once and use the standard "latest" file from the ovation model
if run_realtime == 1:
	run_once = 1
	plot_south = 1
	
output_global_ascii = 1   # Set this to 1 to get gridded ASCII output of global aurora distribution


#  The following start and stop times along with the cadence must match the input data available.   

start_date = dt.datetime(2017,9,28,17,00)   	#YYYY, MM, DD, HH, MM
end_date =  dt.datetime(2017,9,29,0,0)   	#YYYY, MM, DD, HH, MM    
cadence = 60		#minutes
		
	# ******************	 Variable Definition  **************

if run_once == 1:end_date = start_date

run_date = start_date

imin = 0
ihour = 0


#  Relative weights of the three types of electron aurora (Biased toward discrete aurora)
mult_dif = .75
mult_mono = 1.25
mult_wave = 1.

#  Parameter used to Calculate Kp  from Hp
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
					   
		'alpha': ((0.0, 0.6, .6),
 				   (0.1, .9, .9),		
 				   (0.3, 1.0, 1.0),
 				   (0.5, 1.0, 1.0),
 				   (0.8, 1.0, 1.0),
 				   (0.9, 1.0, 1.0),
 				   (1.0, 1.0, 1.0))
 			}

# *****************************   Start of Loop  *****************************
#       If in real-time mode or run-once mode, the code will run staight through and not loop

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
			
		ifl = 'North/aurora_N_' + year + '-' + month + '-' + day + '_' + hour + minute 
		vlat = 80.
		vlon = -100
		if NS == 1: 
			ifl = 'South/aurora_S_' + year + '-' + month + '-' + day + '_' + hour + minute 
			vlat = -90.
			vlon = -90
			
		start_hour = 0
		
		imin = 0
		ihour = 0
	
# 		print (ifl)
		
		
		input_file = open(Input_path + ifl + '.txt', 'r')
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
		
# 		numarray = 7680
			
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
		


#    ******************   Output Data to ASCII Array  *********************
		if output_global_ascii == 1:
			if NS ==0:
				global_array = np.zeros([360,181])
				global_array[0:360,91:181] = geo_aur
			if NS ==1:
				gaf = np.flip(geo_aur, 1)
				global_array[0:360,0:90] = gaf
				global_array = np.flip(np.transpose(global_array),0)
				global_array = global_array *100.
				ofilet = Output_path+"temp.dat"
				# open (ofilet, mode = "w")

			ofilet = Output_path+ifl+'.txt'
			np.savetxt(ofilet, global_array, delimiter=' ', fmt='%5.1f')

#   ***********************************************************************
				
		
			# ***************************************
		
		
		#  Calculate Kp and G scales from Hemispheric Power
		
	
		lnval = ((HP/A1)-y0)
		if lnval <= 0.0: lnval = 1e-5
		Kp= t1*math.log(lnval)
		if Kp < 0.0: Kp = 0.
		
		
		G = 0
		if Kp>=5: G=1
		if Kp>=6: G=2
		if Kp>=7: G=3
		if Kp>=8: G=4
		if Kp>=9: G=5
		
		input_file.close()
		
		now = dt.datetime.now()
		
	# 	************   Set up aurora color map  **********************		
						
		clevs = np.arange(.2,1.1*upper,0.25)
					
		auro_col = LinearSegmentedColormap('aurora_color',cdict1)
		plt.register_cmap(cmap=auro_col)
		cmap=plt.get_cmap('aurora_color')
		
	# ***************	Plot the aurora on a map of the globe**************

		
		height = 50000. #View Height in km
					
		lgb = plt.figure(figsize=(10,10), dpi = 80, facecolor='black')
		
		map = Basemap(projection='nsper',lon_0=vlon,lat_0=vlat,
 			satellite_height=height*1000.,resolution='l')

	
		
		cs=map.nightshade(fdate, color = 'black',delta=0.25, alpha = 0.6, zorder = 2)
		
		landcolor = (.5,.5,.1)
		oceancolor = (0.4, 0.4, 0.8)
		boardercolor = (.7,.7,.7)
		
		map.drawmapboundary(fill_color = oceancolor)
		map.drawcoastlines(color=boardercolor,zorder=5,linewidth = 0.5)
		map.fillcontinents(color=landcolor,lake_color=oceancolor)
		map.drawcountries(color=boardercolor,linewidth=0.5)
		map.drawstates(color=boardercolor,linewidth=0.5)
		map.drawparallels(np.arange(-90.,120.,30.),color=boardercolor)
		map.drawmeridians(np.arange(0.,420.,60.),color=boardercolor)
		
		x,y = map(gx, gy)
		
		#~ clevs = np.arange(1,15,0.5)
		
		cs = map.contourf(x,y,geo_aur,clevs, cmap=cmap , zorder=4)
		

		
# 		  *********************  Lable plot  ***************************
		
		lt = (1/60.)*(fdate - mdate).seconds
		
		#Upper Left
		
		plt.figtext(0.04,0.99,'NOAA Space Weather Prediction Center', color = 'yellow', size =16)
		plt.figtext(0.04,0.953,'Aurora Forecast', color = 'yellow', size = 30,weight = 'bold')
		plt.figtext(0.04,0.925, 'For %s (UTC)' % fdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 16)
		
		#~ Upper Right
		
		plt.figtext(0.99,0.93,'Hemispheric Power: %5.1f GW (Range 5 to 200)' %HP,color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.95,'Estimated Kp: %1i (Range 0 to 9)' %int(Kp),color='yellow',size = 12,ha = 'right')	
		plt.figtext(0.99,0.97,'Estimated G-Scale: %1i (Range 0 to 5)' %int(G),color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.99,'Forecast Lead Time: %3i minutes' %int(lt), color = 'yellow', size =12, ha = 'right')
		
		# ~ Lower Right
		
		plt.figtext(0.99,0.04,'OVATION Aurora Model', color = 'yellow', size =18, ha = 'right')
		plt.figtext(0.99,0.02, 'Model Run at %s (UTC)' % now.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 12,ha = 'right')
		plt.figtext(0.99,0.00, 'L1 Observations at %s (UTC)' % mdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 12,ha = 'right')

# 		#  ************************   Getting rid of boarder   **********************
		plt.gca().set_axis_off()
		plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
		            hspace = 0, wspace = 0)
		plt.margins(0,0)
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		plt.gca().yaxis.set_major_locator(plt.NullLocator())

		
		#  **********************   Add and Label Colorbar.  ***********************

		cax = plt.axes([0.05,0.06,0.30,0.025])
		cb=plt.colorbar(cs,cax=cax,orientation = "horizontal", ticks= [1,2,3,4,5])
	
	
		cb.set_label('Energy Deposition  (erg/cm2)', color=fg_color,size =12)
		cb.ax.xaxis.set_tick_params(color=fg_color)
		cb.ax.set_xticklabels([1,2,3,4],color=fg_color)
		cb.outline.set_edgecolor(fg_color)
		
		cb.set_label('Approximate Energy Deposition', color=fg_color)
 		
		plt.figtext(0.18,.00,'(erg/cm2)', color=fg_color,ha = "center",size =12)
		plt.figtext(0.32,.04,'>4', color=fg_color,size =12)		
		plt.figtext(0.05,.10,'10%        50%       90%', color = 'yellow', size =14)
		plt.figtext(0.2,.13,'Probability of Aurora', color = 'yellow', size =16, ha = "center",weight = 'bold')
				
 			
# 		#**************************  Save image to a temp file  **************************
 		

		ofile1 = Output_path+"temp.png"
 			
		plt.savefig(ofile1 ,facecolor = 'black', edgecolor = 'black',bbox_inches = 'tight', pad_inches = 0)

		
		im = Image.open(ofile1)	
# 
		# ***********************   Add NOAA Logo   **************************	
		
		icon = Image.open("/home/rodney.viereck/python/ovation/configuration/NOAA_logo4.png")
		icon_size = (90,90)
		
		icon2 = icon.resize(icon_size)
		x,y = icon2.size
		
# 		
# 		#NOAA logo position
		px = 20
		py = 75
		im.paste(icon2, (px, py, px + x, py + y), icon2)
		im = im.convert('RGB')
# 		
		ofile = Output_path+ifl+'.png'
		print (ofile)
# 			
		im.save(ofile)
	
	run_date = run_date + dt.timedelta(minutes = cadence)
	
	
	# im2.show()
	
	








