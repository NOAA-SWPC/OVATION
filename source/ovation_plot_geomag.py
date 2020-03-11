import numpy as np
#from mpl_toolkits.basemap import Basemap, cm

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import matplotlib.pyplot as plt
import datetime as dt
import math as math
import scipy.interpolate as sciint

from scipy.ndimage import gaussian_filter
# from aacgmv2 import convert
from matplotlib.colors import LinearSegmentedColormap
from mag2geo import mag2geo
#from PIL import Image
import pdb
import io, json


def cart2polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def polar2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


def ovation_plot_geomag (ipath,ifile,ofile):

# ******************	Open and read the file created by the Ovation Model**************

	#~ Input_path = 'C:/Users/SWPC1/Documents/Python/Ovation_Realtime/Output/Text/'
	#~ Input_file = 'aurora_N1_2017-09-06_0000.txt'

	#~ input_file = open(Input_path + Input_file, 'r')

	if ifile.find('aurora_N') != -1:
		input_file = open(ipath+ 'North/'+ ifile+'.txt','r')
	else:
		input_file = open(ipath+ 'South/'+ ifile+'.txt','r')

	fhead = input_file.readline()

	# *************  Read Observation datetime *******************
	#print("ovation_plot_geomag with ipath: {}, ifile: {}, ofile: {}".format(ipath, ifile, ofile))

	fdatei = input_file.readline() .strip()
	fdatei =fdatei.split()
	fdate1 = fdatei[2].split("-")
	fdate2 = fdatei[3].split(":")

	year = int(fdate1[0])
	month = (int(fdate1[1]))
	day = int(fdate1[2])

	hour = int(fdate2[0])
	min = int(fdate2[1])
	dec_time = float(hour) + float(min)/60.

	mdate = dt.datetime(year, month, day, hour, min)

	# ******************** Read Forecast datetime  ****************

	fdatei = input_file.readline() .strip()
	fdatei =fdatei.split()
	fdate1 = fdatei[2].split("-")
	fdate2 = fdatei[3].split(":")

	year = int(fdate1[0])
	month = (int(fdate1[1]))
	day = int(fdate1[2])

	hour = int(fdate2[0])
	min = int(fdate2[1])
	dec_time = float(hour) + float(min)/60.

	fdate = dt.datetime(year, month, day, hour, min)

	# ******************  Read Hemispheric Power  *****************

	f_hp_head = input_file.readline().strip().split()

	HP = float(f_hp_head[2])

	# Skip the rest of the header

	fhead = input_file.readline()
	fhead = input_file.readline()

	#~ ********************* Read Data ***************
	numarray = 7680

	numlon= 96
	numlat= 80

	mlat = np.zeros(shape=(numlat,numlon+1))
	mlon= np.zeros(shape=(numlat,numlon+1))
	glat = np.zeros(shape=(numlat,numlon+1))
	glon= np.zeros(shape=(numlat,numlon+1))
	alt=np.zeros(shape=(numlat,numlon+1))
# 	aur = np.zeros(shape=(numlat,numlon+1))
	je_d = np.zeros(shape=(numlat,numlon+1))
	je_m = np.zeros(shape=(numlat,numlon+1))
	je_w= np.zeros(shape=(numlat,numlon+1))
	je_i = np.zeros(shape=(numlat,numlon+1))

	lat = np.zeros(numlat)
	lon = np.zeros(numlat)

	itot = 0

	for ilon in range(0,numlon):
		for ilat in range(0,numlat):
			row = input_file.readline().strip().split()
# 			aur[ilat,ilon] = (float(row[2]))
			je_d[ilat,ilon] = (float(row[2]))
			je_m[ilat,ilon] = (float(row[3]))
			je_w[ilat,ilon] = (float(row[4]))
			je_i[ilat,ilon] = (float(row[5]))

			mlat[ilat,ilon] = (float(row[1]))
			if mlat[ilat,ilon] < 0: mlat[ilat,ilon] = -mlat[ilat,ilon]
			mlon[ilat,ilon] = -90. + (float(row[0]))* 15 #-dec_time)*15. + 71.#Convert local time into longitude but keep local noon at time

			itot = itot+1


	for ilat in range(0,numlat):
		mlat[(ilat,96)] = mlat[(ilat,0)]
		mlon[(ilat,96)] = mlon[(ilat,0)]
# 		aur[(ilat,96)] = aur[(ilat,0)]

	#***************   Remove Outliers  *******************************
	#~ mv = np.median(je_w)
	#~ je_w = np.where(je_w < 1000*mv, je_w, mv)		#Remove really big spikes
	#~ temp = 1.3* gaussian_filter(je_w,sigma = 5,mode = 'wrap')	  #Establish a smooth background
	#~ mv = np.median(temp)
	#~ je_w = np.where((je_w - temp) <4*mv,  je_w, temp)			#Replace with smoothed where outlier
	#~ je_w = 1.1* gaussian_filter(je_w,sigma = 1.3,mode = 'wrap')    #Smooth
	#***************************************************************

	#~ plt.imshow(je_w)
	#~ plt.show()
	#~ quit()

	#~ je_w = np.where(je_w < 1000*mv, je_w, mv)

	#~ mv = np.median(je_w)
	#~ print mv

	#~ je_w = np.where(je_w < 300*mv, je_w, mv)

	#~ plt.imshow(je_w)
	#~ plt.show()
	#~ quit()


	#~ print (np.amax(je_d))
	#~ print (np.amax(je_m))
	#~ print (np.amax(je_w))
	#~ print (np.amax(je_i))
	#~ print (np.amax(aur))

	#~ plt.plot(je_d[:,10])
	#~ plt.show()
	#~ quit()

	#~ x, y = np.meshgrid(mlon,mlat)

	theta= (2*np.pi)*mlon/356.
	r = 90-mlat
	X,Y = polar2cart(r,theta)
	#~ X = theta
	#~ Y= r



	#***************   Geomag Plots  *****************************
	cmap = plt.get_cmap('plasma')

# 	aur_col=[]
	je_dif=[]
	je_mon=[]
	je_wav=[]
	je_ion=[]

	patches =[]

	for ix in range (0,numlat-1):
		for iy in range (0,numlon-1):
			vert1 = [X[ix,iy],Y[ix,iy]]
			vert2 = [X[ix+1,iy],Y[ix+1,iy]]
			vert4 = [X[ix,iy+1],Y[ix,iy+1]]
			vert3 = [X[ix+1,iy+1],Y[ix+1,iy+1]]
# 			aur_col.append(aur[ix,iy])
			je_dif.append(je_d[ix,iy])
			je_mon.append(je_m[ix,iy])
			je_wav.append(je_w[ix,iy])
			je_ion.append(je_i[ix,iy])

			polygon = Polygon([vert1,vert2,vert3,vert4],closed=True)
			patches.append(polygon)


	#  ****************	Define color map and color table **************

	cdict1 = {'red':  ((0.0, 0.0, 0.0),
						   (0.1, 0.0, 0.0),
						   (0.25, 0.0, 0.0),
						   (0.5, 0.5, 0.5),
						   (0.75, 1., 1.1),
						   (1.0, 1.0, 1.0)),

				 'green': ((0.0, 0.1, 0.1),
						   (0.1, 0.3, 0.3),
						   (0.25, .6, .6),
						   (0.5, .8, .8),
						   (0.75, 1.0, 1.0),
						   (1.0, 0.0, 0.0)),

				 'blue':  ((0.0, 0.0, 0.0),
							(0.1, 0.0, 0.0),
						   (0.25, 0.0, 0.0),
						   (0.5, 0.0, 0.0),
						   (0.75, 0.0, 0.0),
						   (1.0, 0.0, 0.0)),
					}

				#~ 'alpha': ((0.0, 1.0, 1.0),
						   #~ (0.1, 1.0, 1.0),
						   #~ (0.25, 1.0, 1.0),
						   #~ (0.5, 1.0, 1.0),
						   #~ (0.75, 1.0, 1.0),
						   #~ (1.0, 1.0, 1.0))
					#~ }

	auro_col = LinearSegmentedColormap('aurora_color',cdict1)
	plt.register_cmap(cmap=auro_col)
	colormap=plt.get_cmap('aurora_color')

	plt_col=(.6,.6,.6)
	back_col=(.2,.2,.2)

	#  Set Limits

	pmin = 0.
	pmax_array = [np.max(je_d), np.max(je_m), np.max(je_w), np.max(je_i)]
#	print (pmax_array)
	pmax = 1.1 * np.max(pmax_array)


	fig1 = plt.figure(figsize=[10,8])
	rect = fig1.patch
	rect.set_facecolor(plt_col)

	Z=np.asarray(je_d)
	Z=np.where(Z > 0., Z, 0)
	ax1 = plt.subplot(2,2,1)
	ax1.set_facecolor(back_col)
	ax1 = plt.pcolormesh(X,Y,Z, vmin = pmin, vmax = pmax, cmap=colormap)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.axes.get_yaxis().set_visible(False)
	ax1.axes.labelcolor = 'white'
	circle = plt.Circle((10,10),100, color='white')
	plt.colorbar()

	Z=np.asarray(je_m)
	Z=np.where(Z > 0., Z, 0)
	ax1 = plt.subplot(2,2,2)
	ax1.set_facecolor(back_col)
	ax1 = plt.pcolormesh(X,Y,Z,vmin = pmin, vmax = pmax, cmap=colormap)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.axes.get_yaxis().set_visible(False)
	plt.colorbar()

	Z=np.asarray(je_w)
	Z=np.where(Z > 0., Z, 0)
	ax1 = plt.subplot(2,2,3)
	ax1.set_facecolor(back_col)
	ax1 = plt.pcolormesh(X,Y,Z,vmin = pmin, vmax = pmax, cmap=colormap)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.axes.get_yaxis().set_visible(False)
	plt.colorbar()

	Z=np.asarray(je_i)
	Z=np.where(Z > 0., Z, 0)
	ax1 = plt.subplot(2,2,4)
	ax1.set_facecolor(back_col)
	ax1 = plt.pcolormesh(X,Y,Z,vmin = pmin, vmax = pmax, cmap=colormap)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.axes.get_yaxis().set_visible(False)
	plt.colorbar()

	plt.figtext(0.50,0.97,'OVATION Aurora Model', color='yellow',size = 18,ha = 'center')
	plt.figtext(0.50,0.94, 'L1 Observations at %s (UTC)' % mdate.strftime("%Y-%m-%d %H:%M"),color='yellow',size = 10,ha = 'center')
	plt.figtext(0.50,0.92, 'Hemespheric Power %s' %str(HP),color='yellow',size = 10,ha = 'center')
	plt.figtext(0.50,0.02, 'Plots are in geomagnetic coordinates',color='yellow',size = 10,ha = 'center')

	plt.figtext(.48,.78,'Energy Flux (ergs/$cm^2$)',color='white', size = 10, ha='center',rotation = 270)
	plt.figtext(.905,.78,'Energy Flux (ergs/$cm^2$)',color='white', size = 10, ha='center',rotation = 270)
	plt.figtext(.48,.38,'Energy Flux (ergs/$cm^2$)',color='white', size = 10, ha='center',rotation = 270)
	plt.figtext(.905,.38,'Energy Flux (ergs/$cm^2$)',color='white', size = 10, ha='center',rotation = 270)

	plt.figtext(.27,.9,'Diffuse',color='yellow', size = 12, ha='center')
	plt.figtext(.70,.9,'Monoenergetic',color='yellow', size = 12, ha='center')
	plt.figtext(.27,.48,'Wave',color='yellow', size = 12, ha='center')
	plt.figtext(.70,.48,'Ions',color='yellow', size = 12, ha='center')

	dx = 0.15
	dy = 0.12
	fig_cent = np.array([[.275,.27],[.275,.69],[.70,.27],[.70,.69]])

	#~ i1 = 0
	#~ circle = plt.Circle((fig_cent[i1,:]),1, color='white')

	for ix in range(4):
		yy = fig_cent[ix,0] + dx
		xx = fig_cent[ix,1]

		plt.figtext(xx,yy,'Noon',color='white', size = 8, ha='center')
		yy = fig_cent[ix,0] - dx

		plt.figtext(xx,yy,'Midnight',color='white', size = 8, ha='center')
		xx = fig_cent[ix,1] - dy
		yy = fig_cent[ix,0]

		plt.figtext(xx,yy,'Dusk',color='white', size = 8, ha='center')
		xx = fig_cent[ix,1] + dy

		plt.figtext(xx,yy,'Dawn',color='white', size = 8, ha='center')

	#print("plot ofile: {}".format(ofile))
	plt.savefig(ofile,facecolor = 'black', edgecolor = 'black')
	plt.close()

#	plt.show()
