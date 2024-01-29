# Smoothing routine base on IDL code by Patrick Newell, May 2013.  Assumes no data gaps (so any
# interpolation across DMSP's 0230 MLT wedge must be done separately, preferably first)

from utilities_functions_season_epoch import mlt_plus, mlat_plus
import numpy as np
from scipy.ndimage import gaussian_filter

def simple_smooth_j(j_call):

    return

    j_smootha = j_call

    #~ j_smootha = 1.2* gaussian_filter(j_call,sigma = 3,mode = 'wrap')

    return j_smootha

    #~ nmlt = 96
    #~ nmlat = 160
    #~ j_smootha = np.zeros((nmlt,nmlat))

    #~ for i in range(nmlt):
        #~ for j in range(nmlat):
            #~ ip = int(mlt_plus(i,1))
            #~ im = int(mlt_plus(i,-1))
            #~ jp = int(mlat_plus(j,1))
            #~ jm = int(mlat_plus(j,-1))

            #~ jmm = j_call[im,jm]
            #~ jm0 = j_call[im,j]
            #~ jmp = j_call[im,jp]

            #~ j0m = j_call[i,jm]
            #~ j00 = j_call[i,j]
            #~ j0p = j_call[i,jp]

            #~ jpm = j_call[ip,jm]
            #~ jp0 = j_call[ip,j]
            #~ jpp = j_call[ip,jp]

            #~ w00 = 8.            #centerpoint
            #~ w01 = 2.            #one over in latitude (pos or neg)
            #~ w10 = 3.            #one over in mlt
            #~ w11 = 1.            #one in both mlt and mlat

            #~ wsum = 1*w00 + 2*w01 + 2*w10 + 4*w11

            #~ j_s = w00*j00 + w01*(j0p + j0m) + w10*(jp0 + jm0) + w11*(jmm + jmp + jpm + jpp)
            #~ j_s = j_s/wsum
            #~ j_smootha[i,j] = j_s

    #~ return j_smootha
