#file:///C:/ProgramData/Anaconda3/pkgs/proj4-5.2.0-ha925a31_1/Library/share/epsg

import os
os.environ["PROJ_LIB"] = "C:/ProgramData/Anaconda3/pkgs/proj4-5.2.0-ha925a31_1/Library/share"

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import datetime as dt
import math as math
import scipy.interpolate as sciint
from scipy.ndimage import gaussian_filter
from aacgmv2 import convert
from matplotlib.colors import LinearSegmentedColormap
#~ from mag2geo import mag2geo
from PIL import Image
from rem_out import rem_out

# ******************	Open and read the file created by the Ovation Model**************

run_once = 1  #only plot one map
override = 1  #to specifify an input file
#if override == 1: run_once = 1

override_file = "HF_Kp_4"

#Kp_array = np.zeros((25,360,90))

Input_path = 'C:/Docs/Python/Ovation_New_realtime_2019/Output/Text/'
Output_path = 'C:/Docs/Python/Ovation_New_realtime_2019/Output/HF_Final_Images/'

start_date = dt.datetime(2019,9,28,0,0)   	#YYYY, MM, DD, HH, MM
end_date =  dt.datetime(2019,9,30,0,0)   	#YYYY, MM, DD, HH, MM

if run_once == 1:end_date = start_date


run_date = start_date

imin = 0
ihour = 0

cadence = 60		#minutes

icad = 0

while run_date <= end_date:
	
	delta_t = icad * cadence * 60  #Step size in seconds
		
	smin = str("%2.2i%2.2i"%(ihour,imin))
	
	ifl = 'aurora_N1_2017-09-06_'+smin+'.txt'
	
	year = str("%4.4i"%run_date.year)
	month = str("%2.2i"%run_date.month)
	day = str("%2.2i"%run_date.day)
	hour = str("%2.2i"%run_date.hour)
#	hour = str("%2.2i"%(run_date.hour + delta_t))
	minute = str("%2.2i"%run_date.minute)
	
	ifl = 'aurora_N1_' + year + '-' + month + '-' + day + '_' + hour + minute 
#	ifl = 'aurora_N1_2017-09-07_0000.txt'
	#~ ifl = 'aurora_N_2017-09-07'
	#~ Input_file = ifl +'.txt'
	
	if override == 1: ifl = override_file
	
	print (ifl)

	input_file = open(Input_path+ifl+'.txt', 'r')

	fhead = input_file.readline()

	fdatei = input_file.readline() .strip()
	
#	fdatei =fdatei.split()
#	print ('data date  ',fdatei)
#	fdate1 = fdatei[2].split("-")
#	fdate2 = fdatei[3].split(":")
#
#	year = int(fdate1[0])
#	month = (int(fdate1[1]))
#	day = int(fdate1[2])
#
#	hour = int(fdate2[0])
#	min = int(fdate2[1])
	dec_time = float(hour) + float(min)/60.

	mdate = dt.datetime(int(year), int(month), int(day), int(hour), int(min))
	fdate = mdate

	fdatei = input_file.readline() .strip()
#	fdatei =fdatei.split()
#	fdate1 = fdatei[2].split("-")
#	fdate2 = fdatei[3].split(":")
#
#	year = int(fdate1[0])
#	month = (int(fdate1[1]))
#	day = int(fdate1[2])
#
#	hour = int(fdate2[0])
#	min = int(fdate2[1])

#	fdate = dt.datetime(year, month, day, hour, min)

	f_hp_head = input_file.readline().strip().split()
	HP = float(f_hp_head[2])
	
	f_hp_head = input_file.readline().strip().split()
	Kp = int(float(f_hp_head[2]))

#	fhead = input_file.readline()
	fhead = input_file.readline()

	numarray = 7680

	numlon= 96
	numlat= 80

	glat = np.zeros(shape=(numlat,numlon+1))
	glon= np.zeros(shape=(numlat,numlon+1))
	alt=np.zeros(shape=(numlat,numlon+1))
	aur_diff = np.zeros(shape=(numlat,numlon+1))
	aur_mono = np.zeros(shape=(numlat,numlon+1))
	aur_wave = np.zeros(shape=(numlat,numlon+1))

	lat = np.zeros(numlat)
	lon = np.zeros(numlat)

	itot = 0
	
	offset = 0.
	scale_diff = .01  #(Green)     Diffuse Aurora  (Lower Energy)
	scale_mono = .2  #(Yellow)   Monoenergetic Aurora (Higher Energy but lower Energy Flux)
	scale_wave = .4   #(Red)     Alfven Wave driven aurora (High Energy and Higher Energy Flux)

	for ilon in range(0,numlon):
		for ilat in range(0,numlat):
			itot=itot+1
			row = input_file.readline().strip().split()
			aur_diff[ilat,ilon] = offset + scale_diff * float(row[3])
			aur_mono[ilat,ilon] = offset + scale_mono * float(row[4]) 
			aur_wave[ilat,ilon] = offset + scale_wave * float(row[5]) 
			lat[ilat] = (float(row[1]))
			#Convert local time into longitude
			lon[ilat] = (float(row[0])-dec_time)*15. + 71.
			#~ itot = itot+1
			#~ glat[ilat,ilon] = lat[ilat]
			#~ glon[ilat,ilon] = lon[ilat]
			glat[(ilat,ilon)], glon[(ilat,ilon)] = convert(lat[ilat], lon[ilat],1000, fdate,a2g = True)
	
	input_file.close()
	
	glon[:,96] = glon[:,0]	
	glat[:,96] = glat[:,0]	
	aur_diff[:,96] = aur_diff[:,0]	
	aur_mono[:,96] = aur_mono[:,0]	
	aur_wave[:,96] = aur_wave[:,0]	

	
	aur_diff = rem_out(aur_diff)		#Remove Outliers
	aur_mono = rem_out(aur_mono)
	aur_wave = rem_out(aur_wave)
	
	X= glat
	Y= glon
	
	Z_sum = aur_diff.ravel() + aur_mono.ravel() + aur_wave.ravel()
	
	Z_sum = 3. * np.log(1.0 + (Z_sum - np.min(Z_sum)))
	
	
#	Z_wave = Z_wave - np.min(Z_wave)
#	Z_wave = 8. * (Z_wave / np.max(Z_wave))

	#~ Z=np.where(Z>.01,Z, Z+2)


	# ***************	regrid to regular geographic coordinates **************

#	for ia in range (3):

	glatlon = np.zeros(shape = (((numlat)*(numlon+1)),2))

	glatlon[:,0] = glon.ravel()
	glatlon[:,1] = glat.ravel()

	gx,gy = np.mgrid[-180:180:360j, 0:90:90j]			#Set even lat lon spacing
#		if ia == 0: geo_aur = sciint.griddata(glatlon,Z_diff,(gx,gy),method='cubic')		#Regrid to even lat lon spacing
#		if ia == 1: geo_aur = sciint.griddata(glatlon,Z_mono,(gx,gy),method='cubic')
	geo_aur = sciint.griddata(glatlon,Z_sum,(gx,gy),method='cubic')
	geo_aur[359,:] = geo_aur[358,:]
	geo_aur[0,:] = geo_aur[1,:]
	geo_aur = np.nan_to_num(geo_aur)			#replace NANs with Zeros
	
	geo_aur[geo_aur < 1.5] = 0.	
	geo_aur = np.clip(geo_aur,0.,10.)
	
	geo_aur = 1.2 * gaussian_filter(geo_aur,sigma = 3,mode = 'wrap')		#Smooth data to remove artifacts from geomag to geograph transform
	geo_aur = 1.2 * gaussian_filter(geo_aur,sigma = 3,mode = 'wrap')		#Smooth data to remove artifacts from geomag to geograph transform


#	Kp_array[icad,:,:] = geo_aur        #Set array for longitude plotting
	
#	if ia == 0: geo_aur_diff = geo_aur
#	if ia == 1: geo_aur_mono = geo_aur
#	if ia == 2: geo_aur_wave = geo_aur


	# ***************************************


	#  ***********  Calculate Kp and G  from Hp ************************

#	y0 = 6.9498
#	A1 = 2.04984
#	t1 = 1.500
#
#	lnval = ((HP/A1)-y0)
#	if lnval <= 0.0: lnval = 1e-5
#	Kp= t1*math.log(lnval)
#	
#	if Kp < 0.0: Kp = 0.
#	if Kp > 9.5: Kp = 9.0


	G = 0
	if Kp>=5: G=1
	if Kp>=6: G=2
	if Kp>=7: G=3
	if Kp>=8: G=4
	if Kp>=9: G=5

	now = dt.datetime.now()

	#  ****************	Define color map and color table for HF Absorption**************
			 
	cdict1 = {'red':  ((0.0, 0.0, 0.0),
					   (0.25, 0.,0.),
					   (0.5, .5, .5),
					   (0.75, 1., 1.),
					   (1.0, 1., 1.)),

			 'green': ((0.0, 0.25, 0.25),
					   (0.25, .5, .5),
					   (0.5, 1., 1.),
					   (0.75, 1., 1.),
					   (1.0, 0., 0.)),

			 'blue':  ((0.0, 0.0, 0.0),
					   (0.25, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.75, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),
					   
			'alpha': ((0.0, 0.5, .75),		
					   (0.25, .75, 1.),
					   (0.5, 1.0, 1.0),
					   (0.75, 1.0, 1.0),
					   (1.0, 1.0, 1.0))
				}

	#  ****************	Define color map and color table for electron bands**************
			 
	cdictr = {'red':  ((0.0, 0.4, 0.4),
					   (0.1, 0.8, 0.8),
					   (0.25, 1., 1.),
					   (0.5,  1., 1.),
					   (0.75, 1., 1.),
					   (1.0,  1., 1.)),

			 'green': ((0.0, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.25, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.75, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),

			 'blue':  ((0.0, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.25, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.75, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),
					   
			'alpha': ((0.0, 0.0, .25),
					   (0.1, 0.25, 0.5),			
					   (0.25, .5, .6),
					   (0.5, .6, .7),
					   (0.75, .8, .9),
					   (1.0, 1.0, 1.0))
				}
	
	cdictg = {'red':  ((0.0, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.25, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.75, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),

			 'green': ((0.0, 0.0, 0.2),
					   (0.1, 0.2, 0.3),
					   (0.25, 3., 4.),
					   (0.5,  4., 6.),
					   (0.75, 6., 8.),
					   (1.0,  8., 1.)),

			 'blue':  ((0.0, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.25, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.75, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),
					   
			'alpha': ((0.0, 0.0, .25),
					   (0.1, 0.25, 0.5),			
					   (0.25, .5, .6),
					   (0.5, .6, .7),
					   (0.75, .8, .9),
					   (1.0, 1.0, 1.0))
				}
		
	cdictp = {'red':  ((0.0, 0.0, 0.2),
					   (0.1, 0.2, 0.3),
					   (0.25, 3., 4.),
					   (0.5,  4., 6.),
					   (0.75, 6., 8.),
					   (1.0,  8., 1.)),

			 'green': ((0.0, 0.0, 0.2),
					   (0.1, 0.2, 0.3),
					   (0.25, 3., 4.),
					   (0.5,  4., 6.),
					   (0.75, 6., 8.),
					   (1.0,  8., 1.)),

			 'blue':  ((0.0, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.25, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.75, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),
					   
			'alpha': ((0.0, 0.0, .25),
					   (0.1, 0.25, 0.5),			
					   (0.25, .5, .6),
					   (0.5, .6, .7),
					   (0.75, .8, .9),
					   (1.0, 1.0, 1.0))
				}
				
		
#	col_diff = LinearSegmentedColormap('color_diff',cdictg)
#	plt.register_cmap(cmap=col_diff)
#	cmap_diff=plt.get_cmap('color_diff')
#	
#	col_mono = LinearSegmentedColormap('color_mono',cdictp)
#	plt.register_cmap(cmap=col_mono)
#	cmap_mono=plt.get_cmap('color_mono')

	col_wave = LinearSegmentedColormap('color_HF',cdict1)
	plt.register_cmap(cmap=col_wave)
#	cmap_wave=plt.get_cmap('color_wave')
	cmap_wave=plt.get_cmap('color_HF')	

	# ***************	Plot the aurora on a map of the globe**************

	height = 50000. #View Height in km
			
	lgb = plt.figure(figsize=(10,10),facecolor='black')

	map = Basemap(projection='nsper',lon_0=-105,lat_0=80,
		satellite_height=height*1000.,resolution='l')


	cs=map.nightshade(fdate, color = 'black',delta=0.25, alpha = 0.6, zorder = 2)

	landcolor = (.3,.3,.05)
	oceancolor = (0.2, 0.2, 0.6)

	map.drawmapboundary(fill_color = oceancolor)
	map.drawcoastlines(color='gray',zorder=5,linewidth = 0.5)
	map.fillcontinents(color=landcolor,lake_color='darkblue')
	map.drawcountries(color='gray',zorder=5,linewidth=0.5)
	map.drawstates(color='gray',zorder=5,linewidth=0.5)
	map.drawparallels(np.arange(-90.,120.,30.),zorder=5,color='gray')
	map.drawmeridians(np.arange(0.,420.,60.),zorder=5,color='gray')

	x,y = map(gx, gy)

	clevs = np.arange(1,11,1.)

#	print ('Max Value ',np.max(geo_aur))
#	geo_aur_diff = (geo_aur_diff + 0.5) * 0.90  # Amplify lower values for plotting
#	geo_aur_diff[geo_aur_diff < 0.5] = 0.	
#	geo_aur_diff = np.clip(geo_aur_diff,0,3.75)
	
#	geo_aur_mono = (geo_aur_mono + 0.5) * 0.90  # Amplify lower values for plotting
#	geo_aur_mono[geo_aur_mono < 0.5] = 0.	
#	geo_aur_mono = np.clip(geo_aur_mono,0,3.75)
	
#	geo_aur_wave = (geo_aur_wave + 0.5) * 0.90  # Amplify lower values for plotting
#	geo_aur_wave[geo_aur_wave < 1.5] = 0.	
#	geo_aur_wave = np.clip(geo_aur_wave,0.,10.)
#	print ('Max Value ',np.max(geo_aur))
	
	
#	cs = map.contourf(x,y,geo_aur_wave,clevs, cmap=cmap , zorder=4)
	
#	cs = map.contourf(x,y,geo_aur_diff,clevs, cmap=cmap_diff , zorder=4)
#	cs = map.contourf(x,y,geo_aur_mono,clevs, cmap=cmap_mono , zorder=4)
	
	cs = map.contourf(x,y,geo_aur,clevs, cmap=cmap_wave , zorder=4)


	#  **********************   Add and Label Colorbar.  ***********************

	cax = plt.axes([0.15,0.16,0.30,0.025])
	cb=plt.colorbar(cs,cax=cax,orientation = "horizontal", ticks= [1,2,3,4,5,6,7,8,9,10])

	fg_color = 'yellow'

#	cb.set_label('HF Radio Impact', color=fg_color)
#	cb.ax.xaxis.set_tick_params(color=fg_color)
	cb.ax.set_xticklabels([1,2,3,4,5,6,7,8,9,10],color='yellow')
	cb.outline.set_edgecolor(fg_color)
#	plt.figtext(0.30,.19,'10%                      50%                     90%', color = 'yellow', size =10, ha = "center")
	plt.figtext(0.30,.19,'None          Minor         Moderate           Severe', color = 'yellow', size =10, ha = "center")
	plt.figtext(0.30,.21,'HF Radio Impact', color = 'yellow', size =12, ha = "center")


	#  *********************  Lable plot  ***************************

	#Upper Left

	plt.figtext(0.14,0.85,'HF Radio Communications', color = 'yellow', size =22, weight = 'bold')
	plt.figtext(0.14,0.825,'Model: OVATION-Py 2013', color = 'yellow', size =14)
	plt.figtext(0.14,0.80,'15 - 30 Minutes Forecast Lead Time', color = 'yellow', size =14)
	

#	plt.figtext(0.05,0.50,'Differential', color = 'Green', size =12)
#	plt.figtext(0.05,0.48,'Monoenergetic', color = 'Yellow', size =12)
#	plt.figtext(0.05,0.46,'Wave', color = 'Red', size =12)	
	
	#~ Upper Right

	plt.figtext(0.90,0.83, 'Forecast for %s (UTC)' % fdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 12,ha = 'right')
#	plt.figtext(0.90,0.84,'Hemispheric Power: %5.1f GW (Range 5 to 200)' %HP,color='yellow',size = 12,ha = 'right')
	plt.figtext(0.90,0.81,'Estimated Kp: %1i (Range 0 to 9)' %int(Kp),color='yellow',size = 12,ha = 'right')	
	plt.figtext(0.90,0.79,'Estimated G-Scale: %1i (Range 0 to 5)' %int(G),color='yellow',size = 12,ha = 'right')
	plt.figtext(0.90,0.77,'Hemispheric Power: %5.1f GW (Range 5 to 200)' %HP,color='yellow',size = 12,ha = 'right')

	# ~ Lower Right

	plt.figtext(0.90,0.14, 'Model Run at %s (UTC)' % now.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 10,ha = 'right')
	plt.figtext(0.90,0.12, 'L1 Observations at %s (UTC)' % mdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 10,ha = 'right')

	#~ plt.show()
	# ***********************   Crop and Save File   **************************	

	#~ ofile1 = Output_path+"temp"+smin+".jpg"
	ofile1 = Output_path+ifl+'_'+str("%2.2i" %icad)+'.jpg'
	print (ofile1)
	
	plt.savefig(ofile1,facecolor = 'black', edgecolor = 'black')
	
	plt.close()
	
	run_date = run_date + dt.timedelta(minutes = cadence)
	
	icad = icad + 1    

#Crop and resave image for web page
	
#	im = Image.open(ofile1)
#	im2 = im.crop((115,110,915,910))
#	im2 = im2.convert('RGB')
#	ofile = Output_path+ifl+'_resz.jpg'
#	print (ofile)
	
#	im2.save(ofile)












