import numpy as np

from scipy.ndimage import gaussian_filter

def spike_removal(arr):

	temp1 = 1.2 * gaussian_filter(arr,sigma = 3,mode = 'wrap')
	#~ temp2 = arr - temp1
	mean_val = np.mean(arr)
	np.where(arr > 10* mean_val,arr,temp1)
	
	return arr
	