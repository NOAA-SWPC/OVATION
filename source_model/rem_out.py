#  rem_out  removes outliers from a 2 D array
import numpy as np
from scipy.ndimage.filters import uniform_filter, gaussian_filter, maximum_filter, generic_filter

def rem_out (array_in):

	#***************   Remove Outliers  *******************************
	#~ ax1 = plt.subplot(2,1,1)
	#~ ax1 = plt.imshow(je[i1,:,:])
	#~ print(type(je))
	mv = np.median(array_in)
	array_out = array_in
	#~ array_out = np.where(array_in < 1000*mv, array_in, mv)		#Remove really big spikes		
	temp = 1.2* gaussian_filter(array_out,sigma = 3,mode = 'wrap')	  #Establish a smooth background
	temp = 1.2* gaussian_filter(temp,sigma = 3,mode = 'wrap')
	mv = np.median(temp)
	array_out = np.where((array_in - temp) <4*mv,  array_out, temp)			#Replace with smoothed where outlier
	array_out = 1.2* gaussian_filter(array_out,sigma = 3,mode = 'wrap')    #Smooth the final answer
	#~ ax1 = plt.subplot(2,1,2)
	#~ plt.imshow(je[i1,:,:])
	#~ plt.show()
	
	return array_out