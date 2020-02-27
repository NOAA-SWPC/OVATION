	########################################################
	# writing data to text file
	########################################################

import numpy as np
import datetime as dt 

#from ovation_plot_geomag_new import ovation_plot_geomag
#from ovation_plot_new import ovation_plot

def write_ascii_file(mode,NS,Output_Path,time_sw,time_for,time_lab,mlt_array,mlat_array,je_d,je_m,je_w,je_i,power_hemi,Kp_2):

	header0 = 'OVATION Prime (2013 Version) \n'
	header1 = 'Observation Time: ' + str(time_sw) + '\nForecast Time:    '+  str(time_for) + '\n'
	header2 = 'Hemispheric Power:   ' + str('%6.1f' %power_hemi) + '\nForecast Kp:   ' + str(Kp_2)
	header3 = '\n Local_time   Lat_mag       je_diff     je_mono     je_wav      je_ions \n'
	
	header_all = header0 + header1 + header2 + header3

	
	#~ print (je_s)
	zipped = zip(mlt_array, mlat_array, je_d,je_m,je_w,je_i) 
	
	nrows = len(je_d)
	
	
	if1 = 'aurora_S_'
	ipath = 'South/'		
	if NS == 0: 
		if1 = 'aurora_N_'
		ipath = 'North/'
	
	ofile_name = if1 + str(time_lab)[0:10] + '_' + str(time_lab)[11:13] + str(time_lab)[14:16]
	ofile_name_2 = "Ovation_Latest_"



	print ('ASCII Output File ', Output_Path + ipath + ofile_name +'.txt')
	
	je_w = np.clip(je_w,0,np.amax(je_w))
	je_d = np.clip(je_d,0,np.amax(je_d))
	je_m = np.clip(je_m,0,np.amax(je_m))
	je_i = np.clip(je_i,0,np.amax(je_i))
	
	
	ofile = open(Output_Path + ipath + ofile_name +'.txt',"w")
	ofile.write(header_all)

	for iloop in range(nrows): ofile.write(' %8.4e   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n' %(mlt_array[iloop], mlat_array[iloop], je_d[iloop],je_m[iloop],je_w[iloop],je_i[iloop]))
	
	ofile.close()
	
	if mode == 'NOWCAST':
		ofile2 = open(Output_Path + ofile_name_2 + if1 + '.txt',"w")
		ofile2.write(header_all)
	
		for iloop in range(nrows): ofile2.write(' %8.4e   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n' %(mlt_array[iloop], mlat_array[iloop], je_d[iloop],je_m[iloop],je_w[iloop],je_i[iloop]))
		
		ofile2.close()
	
#	ovation_plot_geomag(Output_Path,ofile_name)
	#~ ovation_plot(Output_Path,ofile_name)
	
	#~ else:
		#~ np.savetxt(Output_Path + 'aurora_S_'+str(time_sw)[0:10] + '_' + str(time_sw)[11:13] + str(time_sw)[14:16]+'.txt', zipped, fmt = '%10.4f', header = header_all)
	
	return(ofile_name)