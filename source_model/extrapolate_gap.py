# extrapolate across the DMSP gap from 55-68 MLAT and 0015-0330 MLT
# based on IDL code by Patrick Newell, June 2012
#
# Author: Diana Morosan (morosand@tcd.ie)

from utilities_functions_season_epoch import mlt_bin, mlat_bin
import numpy as np

def extrapolate_gap(a_in):

    #~ print ('Entering extrapolate_gap')
    
    nmlt = 96 # 24*4=96
    nmlat = 160
    delta_mlat = (90.-50.)/float(nmlat)
    delta_mlt = 24./float(nmlt)

    mlat_start = 58.1
    mlat_end = 75.
    i00 = mlt_bin(23.8) #last two good bins
    i01 = mlt_bin(0.1)
    i10 = mlt_bin(3.6)   #first good bin
    i11 = mlt_bin(3.85)  #next good bin
    nsteps = (i10 - i01) - 1    #number of bins to fill in

    mlatc = mlat_start

    while mlatc < mlat_end:
        
        jn_mlat = mlat_bin(mlatc) # north
        js_mlat = mlat_bin(-mlatc) # south
        #linearly interpolate from before (0) to after (1) the gap
        a_start_n00 = a_in[i00,jn_mlat]
        a_start_s00 = a_in[i00,js_mlat]
        a_start_n01 = a_in[i01,jn_mlat]
        a_start_s01 = a_in[i01,js_mlat]
                        
        a_end_n10 = a_in[i10,jn_mlat]
        a_end_s10 = a_in[i10,js_mlat]
        a_end_n11 = a_in[i11,jn_mlat]
        a_end_s11 = a_in[i11,js_mlat]

        a_start_00 = (a_start_n00 + a_start_s00)/2.
        if a_start_n00*a_start_s00 == 0.:
            a_start_00=a_start_n00+a_start_s00
        a_start_01 = (a_start_n01 + a_start_s01)/2.
        if a_start_n01*a_start_s01 == 0.:
            a_start_01=a_start_n01+a_start_s01
        a_start = (a_start_00 + a_start_01)/2.
        if a_start_00*a_start_01 == 0.:
            a_start_00 = a_start_00 + a_start_01

        a_end_10 = (a_end_n10 + a_end_s10)/2.
        if a_end_n10*a_end_s10 == 0.:
            a_end_10 = a_end_n10 + a_end_s10
        a_end_11 = (a_end_n11 + a_end_s11)/2.
        if a_end_n11*a_end_s11 == 0.:
            a_end_11 = a_end_n11 + a_end_s11
        a_end = (a_end_10 + a_end_11)/2.
        if a_end_10*a_end_11 == 0.:
            a_end = a_end_10 + a_end_11

        for i in range(i01 + 1, i10 -1):
            frac1 = float(i-i01)/float(i10 - i01)
            frac0 = 1. - frac1
            a_inter = frac0*a_start + frac1*a_end
            a_in[i,jn_mlat] = a_inter
            a_in[i,js_mlat] = a_inter

        mlatc = mlatc + delta_mlat

    return

