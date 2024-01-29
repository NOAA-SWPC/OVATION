# Salt and pepper noise removal. Assumes no data gaps (so any
# interpolation across DMSP's 0230 MLT wedge must be done separately)
# should work directly on correlation coefficients.
# Noise is determined from b1 and b2 (auroral coefficients).
# Probability is adjusted in whatever way b1 and b2 are
# includes number flux (which passively follows energy flux smoothing)
#
# Based on IDL code by Patrick Newell, November 2012
#
# Author: Diana Morosan
# morosand@tcd.ie

import numpy as np
from utilities_functions_season_epoch import mlt_plus, mlat_plus

def s_and_p(b1,b2,b1p,b2p,b1n,b2n):
    
    print ('Enter s_and_p   ')
    
    # model grid size - always fixed
    nmlt = 96
    nmlat = 160
    dfdt = 4300.            #reasonableness is evaluated at mean driving

    b1_out = np.zeros((nmlt,nmlat))
    b2_out = np.zeros((nmlt,nmlat))
    b1p_out = np.zeros((nmlt,nmlat))
    b2p_out = np.zeros((nmlt,nmlat))
    b1n_out = np.zeros((nmlt,nmlat))
    b2n_out = np.zeros((nmlt,nmlat))

    for i in range(nmlt):
        for j in range(nmlat):
            ip = int(mlt_plus(i,1))
            im = int(mlt_plus(i,-1))
            jp = int(mlat_plus(j,1))    #Not sure why it changed to float around loop 80
            jm = int(mlat_plus(j,-1))

            je_m0 = b2[im,j]*dfdt + b1[im,j]
            je_00 = b2[i,j]*dfdt + b1[i,j]
            je_p0 = b2[ip,j]*dfdt + b1[ip,j]

            je_0m = b2[i,jm]*dfdt + b1[i,jm]
            je_0p = b2[i,jp]*dfdt + b1[i,jp]

            je_mm = b2[im,jm]*dfdt + b1[im,jm]
            je_mp = b2[im,jp]*dfdt + b1[im,jp]
            je_pm = b2[ip,jm]*dfdt + b1[ip,jm]
            je_pp = b2[ip,jp]*dfdt + b1[ip,jp]

            je_all = np.asarray([je_m0,je_p0,je_0m,je_0p,je_mm,je_mp,je_pm,je_pp])

            i_nz = np.where( je_all != 0. )[0]
            je_all = je_all[i_nz]   #exactly zero means no data
            i_neg = np.where ( je_all < 0. )[0]
            je_all[i_neg] = 0.      #negative values have data (but ap >=0)
            n_left = len(je_all)

            if n_left >= 4: # at least 4 neighbor points with data
                je_mean = np.mean(je_all)
                je_sd = np.var(je_all)
                var = np.absolute(je_00 - je_mean)
                if (var < 3.*je_sd) or (var < 0.05):
                    b1_out[i,j] = b1[i,j]
                    b2_out[i,j] = b2[i,j]
                    b1p_out[i,j] = b1p[i,j]
                    b2p_out[i,j] = b2p[i,j]
                    b1n_out[i,j] = b1n[i,j]
                    b2n_out[i,j] = b2n[i,j]
                else:
                    b1_all = [ b1[i,jm],b1[i,jp],b1[im,j],b1[ip,j],b1[im,jm],b1[im,jp],b1[ip,jm],b1[ip,jp] ]
                    b2_all = [ b2[i,jm],b2[i,jp],b2[im,j],b2[ip,j],b2[im,jm],b2[im,jp],b2[ip,jm],b2[ip,jp] ]
                    b1p_all = [b1p[i,jm],b1p[i,jp],b1p[im,j],b1p[ip,j],b1p[im,jm],b1p[im,jp],b1p[ip,jm],b1p[ip,jp] ]
                    b2p_all = [b2p[i,jm],b2p[i,jp],b2p[im,j],b2p[ip,j],b2p[im,jm],b2p[im,jp],b2p[ip,jm],b2p[ip,jp] ]
                    b1n_all = [ b1n[i,jm],b1n[i,jp],b1n[im,j],b1n[ip,j],b1n[im,jm],b1n[im,jp],b1n[ip,jm],b1n[ip,jp] ]
                    b2n_all = [ b2n[i,jm],b2n[i,jp],b2n[im,j],b2n[ip,j],b2n[im,jm],b2n[im,jp],b2n[ip,jm],b2n[ip,jp] ]

                    b1_out[i,j] = np.mean( b1_all )
                    b2_out[i,j] = np.mean( b2_all )
                    b1p_out[i,j] = np.mean(b1p_all)
                    b2p_out[i,j] = np.mean(b2p_all)
                    b1n_out[i,j] = np.mean(b1n_all)
                    b2n_out[i,j] = np.mean(b2n_all)

    # save in temp arrays, since we are creating a new out set
    b1t = b1_out
    b2t = b2_out
    b1pt = b1p_out
    b2pt = b2p_out
    b1nt = b1n_out
    b2nt = b2n_out

    # remove any latitudinally narrow spike (1 or 2 high bins surrounded
    # by at least 3 low bins on either side
    for i in range(nmlt):
        for j in range(3,nmlat-3):
 
            jp1 = int(mlat_plus(j,1))
            jm1 = int(mlat_plus(j,-1))
            jp2 = int(mlat_plus(j,2))
            jm2 = int(mlat_plus(j,-2))
            jp3 = int(mlat_plus(j,3))
            jm3 = int(mlat_plus(j,-3))
            jp4 = int(mlat_plus(j,4))
            jm4 = int(mlat_plus(j,-4))

            je_0 = b2t[i,j]*dfdt + b1t[i,j]

            je_m1 = b2t[i,jm1]*dfdt + b1t[i,jm1]
            je_p1 = b2t[i,jp1]*dfdt + b1t[i,jp1]
            je_m2 = b2t[i,jm2]*dfdt + b1t[i,jm2]
            je_p2 = b2t[i,jp2]*dfdt + b1t[i,jp2]
            je_m3 = b2t[i,jm3]*dfdt + b1t[i,jm3]
            je_p3 = b2t[i,jp3]*dfdt + b1t[i,jp3]
            je_m4 = b2t[i,jm4]*dfdt + b1t[i,jm4]
            je_p4 = b2t[i,jp4]*dfdt + b1t[i,jp4]

            je_low_all = np.array([je_m4,je_m3,je_m2])
            je_high_all = np.array([je_p2,je_p3,je_p4])

            i_nz_low = np.where( je_low_all != 0. )
            je_low_all = je_low_all[i_nz_low]         #exactly zero means no data
            i_neg = np.where( je_low_all < 0. )

            if i_neg[0].all() >= 0:
               je_low_all[i_neg] = 0.
            i_nz_high = np.where( je_high_all != 0. )
            je_high_all = je_high_all[i_nz_high]
            i_neg = np.where( je_high_all < 0. )

            if i_neg[0].all() >= 0:
               je_high_all[i_neg] = 0.

            n_low_left = len(je_low_all)
            n_high_left = len(je_high_all)

            if (n_low_left >= 2) and (n_high_left >= 2):
                je_low_mean = np.mean(je_low_all)
                je_high_mean = np.mean(je_high_all)
                if (je_low_mean < 0.7*je_0) and (je_high_mean < 0.7*je_0):
                    b1_all = [b1t[i,jm2],b1t[i,jm3],b1t[i,jm4],b1t[i,jp2],b1t[i,jp3], b1t[i,jp4] ]
                    b2_all = [b2t[i,jm2],b2t[i,jm3],b2t[i,jm4],b2t[i,jp2],b2t[i,jp3], b2t[i,jp4] ]
                    b1p_all= [b1pt[i,jm2],b1pt[i,jm3],b1pt[i,jm4],b1pt[i,jp2], b1pt[i,jp3], b1pt[i,jp4] ]
                    b2p_all= [b2pt[i,jm2],b2pt[i,jm3],b2pt[i,jm4],b2pt[i,jp2],b2pt[i,jp3], b2pt[i,jp4] ]
                    b1n_all = [b1nt[i,jm2],b1nt[i,jm3],b1nt[i,jm4],b1nt[i,jp2],b1nt[i,jp3], b1nt[i,jp4] ]
                    b2n_all = [b2nt[i,jm2],b2nt[i,jm3],b2nt[i,jm4],b2nt[i,jp2],b2nt[i,jp3], b2nt[i,jp4] ]
                        
                    b1_out[i,j] = np.mean( b1_all )
                    b2_out[i,j] = np.mean( b2_all )
                    b1p_out[i,j] = np.mean(b1p_all)
                    b2p_out[i,j] = np.mean(b2p_all)
                    b1n_out[i,j] = np.mean(b1n_all)
                    b2n_out[i,j] = np.mean(b2n_all)

    return b1_out,b2_out,b1p_out,b2p_out,b1n_out,b2n_out











