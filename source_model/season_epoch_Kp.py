# season_epoch.py:  Python program to make a precipitation map appropriate to a
# specified date and time -- with seasonal variations incorporated
# The focus is on the combined (two hemisphere) result, which then
# applies to the Northern Hemisphere.  The southern hemisphere can be
# taken to be 182 days out of phase with the called time
#
# Author: Diana Morosan
# morosand@tcd.ie
#
# Based on IDL code by Patrick Newell, 2009 and
# Ovation Model Operational IDL code, 2013


import datetime
import numpy as np
#import matplotlib.pyplot as plt
import pickle

from swpc_get_solar_wind_realtime_data_json import swpc_get_solar_wind_realtime_data_json
from solar_wind_coupling import sol_coup
from utilities_functions_season_epoch import prob_estimate, season_weights, mlt_bin, mlat_bin
from extrapolate_gap import extrapolate_gap
from simple_smooth_j import simple_smooth_j
from s_and_p import s_and_p
from spike_removal import spike_removal
from swpc_get_Kp_NOWCAST_data_json import swpc_get_Kp_NOWCAST_data_json


import cProfile
 
# pr = cProfile.Profile()
# pr.enable()
 
# call_function()
 
# pr.disable()
 
# pr.print_stats(sort='time')

## want to extract this to a grid_def.py or something like that ##
nmlt = 96                    #15 min MLT grid
nmlat = 160                #50-90 south, 50-90 north, step 0.5
ndF = 12                      #dFdt bins in Prob
hmlat = int((nmlat-1)/2)

# paths needed to access the regression fits: premodel
# premodel_path = 'premodel/'

# Interpolated solar wind hourly averages
# This auroral power version uses a weighted average over last 4 hours
# Takes as input solar wind structure from realtime data
#
# input: sw_avg, type:'dict'
#
# output:
# Author: Diana Morosan
# morosand@tcd.ie
#
# Based on Ovation Model Operational IDL code

def ap_inter_sol_realtime(sw_avg):
    
    #~ print ('Entering ap_inter_sol_realtime')

    nh = 4     # hours previous to integrate over
    wh = 0.65   # reduce weighting by factor of wh each hour back
    ncoup = 33   # number of coupling functions

    # Compute coupling functions
    Ec_arr = np.zeros((nh,ncoup)) # coupling functions
    w_arr = [] # weight arrays

    for i in range(nh):
        # finding coupling function for each hourly average and appending it to Ec array
        Ec_arr[i,:] = ( sol_coup(sw_avg['Bx'][i], sw_avg['By'][i], sw_avg['Bz'][i],sw_avg['v'][i],sw_avg['ni'][i]) )
        w_arr.append(wh**i)

    #~ print('Rewturning from sol_coup')
    
    wt = np.nansum( w_arr ) # excluding nans
    w_arr = w_arr/wt


    # Ouput Parameters
    Bmag = np.nansum(w_arr*sw_avg['B_average'])
    Bx = np.nansum(w_arr*sw_avg['Bx'])
    By = np.nansum(w_arr*sw_avg['By'])
    Bz = np.nansum(w_arr*sw_avg['Bz'])
    v = np.nansum(w_arr*sw_avg['v'])
    ni = np.nansum(w_arr*sw_avg['ni'])

    Ec = []
    for i in range(ncoup):
        Ec.append( np.nansum( w_arr*Ec_arr[:,i] ))

    return Bmag, Bx, By, Bz, v, ni, Ec


# get_dmsp_smooth.py: returns jea (DMSP) appropriate to a
# specified date and time -- with seasonal variations incorporated
# The array returned applies to the Northern Hemisphere.
# For southern hemisphere call with (doy+182) mod 365
# This smooth version applies a salt and pepper noise removal and
# nearest neighbor smoothing
# Based on IDL code by Patrick Newell. Last updated May 2013

def get_dmsp_data(dmsp_path):

    my_pickle = '{}/dmsp.pkl'.format(dmsp_path)

    try:
        with open(my_pickle, 'rb') as f:
#            (b1p, b2p, b1a, b2a, b1n, b2n, prob_all) = pickle.load(f)
            return pickle.load(f)
    except Exception as e:
        pass

    sname_a = ['winter','spring','summer','fall']
    atype_string = ['diff','mono','wave','ions']
    aname = ['diff_','mono_','wave_','ions_']

    # When probability fits fail, this array is used as a backup
    prob_all = np.zeros((4,3,ndF,nmlt,nmlat))
    b1p = np.zeros((4,3,nmlt,nmlat))
    b2p = np.zeros((4,3,nmlt,nmlat))
    b1a = np.zeros((4,4,nmlt,nmlat))
    b2a = np.zeros((4,4,nmlt,nmlat))
    b1n = np.zeros((4,4,nmlt,nmlat))
    b2n = np.zeros((4,4,nmlt,nmlat)) # iseason,atype,nmlt,nmlat

    for atype in range(4): #loop on electron auroral types
        for iseason in range(4):
            sname = sname_a[iseason]
            afile = dmsp_path + sname + '_' +  atype_string[atype]
            nafile = afile + '_n'
            afile = afile + '.txt'
            nafile = nafile + '.txt'
            pfile = dmsp_path + sname + '_' + 'prob_b_'
            pfile = pfile + atype_string[atype] + '.txt'

            #read regression coefficients for auroral type probabilities
            #ions (atype=3) are not classified, so probability = 1.0 always
            if atype <= 2:
                b1 = 0.
                b2 = 0.
                yend = 1900
                dend = 1
                y0 = 1900
                d0 = 1
                files_done = 0
                sf0 = 0
                # reading what's in the file
                with open( pfile, 'r') as f:
                    # reading just the first line in the file
                    fline=f.readline().rstrip()
                    #this header gives model compilation statistics.  Model was originally
                    #complied from y0,d0 through yend,dend.  files_done is the number of
                    #DMSP satellite-day files = sum over (days x satellites).  sf0 is also
                    #book keeping, season flag refering to model development, not end usage
                    y0, d0, yend, dend, files_done, sf0 = np.array(fline.split())

                    # reading line by line
                    for i in range(nmlt):
                        for j in range(nmlat):
                            l = f.readline().rstrip()
                            b1, b2 = np.array(l.split())
                            b1p[iseason,atype,i,j] = b1
                            b2p[iseason,atype,i,j] = b2
                    #first resort is coefficients above.  Probability table below is just
                    #for back-up
                    for i in range(nmlt):
                        for j in range(nmlat):
                            l = f.readline().rstrip()
                            ## this absolutely does not work as intended ##
                            prob_all[iseason,atype,0:ndF-1,i,j] = float(l)

            # now read in regression coefficients for auroral flux

            with open( afile, 'r') as f:
                # reading just the first line in the file
                fline=f.readline().rstrip()
                y0, d0, yend, dend, files_done, sf0 = np.array(fline.split())
                i = 0
                j = 0
                b1 = 0.
                b2 = 0.
                rFa = np.zeros((nmlt,nmlat))
                while True:
                    l=f.readline().rstrip()
                    #Prob[0]print 'line', l, np.array(l.split())
                    if l == '':
                        # either end of file or just a blank line.....
                        # we'll assume EOF
                        break
                    i,j,b1,b2,rF = np.array(l.split())
                    # turning into int to use as indices
                    i, j = int(i), int(j)
                    rFa[i,j] = rF
                    b1a[iseason,atype,i,j] = b1
                    b2a[iseason,atype,i,j] = b2

            #number flux regress
            with open( nafile, 'r') as f:
                fline=f.readline().rstrip()
                y0, d0, yend, dend, files_done, sf0 = np.array(fline.split())
                i = 0
                j = 0
                b1 = 0.
                b2 = 0.
                rFan = np.zeros((nmlt,nmlat))
                while True:
                    l=f.readline().rstrip()
                    #Prob[0]print 'line', l, np.array(l.split())
                    if l == '':
                        # either end of file or just a blank line.....
                        # we'll assume EOF
                        break
                    i,j,b1,b2,rF = np.array(l.split())
                    # turning into int to use as indices
                    i, j = int(i), int(j)
                    rFan[i,j] = rF
                    b1n[iseason,atype,i,j] = b1
                    b2n[iseason,atype,i,j] = b2

            # remove salt and pepper noise and close postmidnight data gap
            b1a_temp = b1a[iseason,atype]
            b2a_temp = b2a[iseason,atype]
            b1n_temp = b1n[iseason,atype]
            b2n_temp = b2n[iseason,atype]

            if atype < 2:
                b1p_temp = b1p[iseason,atype]
                b2p_temp = b2p[iseason,atype]

            #salt and pepper noise removal
            #~ b1a_s,b2a_s,b1p_s,b2p_s,b1n_s,b2n_s = s_and_p(b1a_temp,b2a_temp,b1p_temp,b2p_temp,b1n_temp,b2n_temp)

            #use function in extrapolate_gap.py to close gap
            extrapolate_gap(b1a_temp)
            extrapolate_gap(b2a_temp)
            extrapolate_gap(b1n_temp)
            extrapolate_gap(b2n_temp)

            if atype < 2:
                extrapolate_gap(b1p_temp)
                extrapolate_gap(b2p_temp)

            # smoothing
            #~ b1a_smooth = simple_smooth_j(b1a_temp)
            #~ b2a_smooth = simple_smooth_j(b2a_temp)
            #~ b1n_smooth = simple_smooth_j(b1n_temp)
            #~ b2n_smooth = simple_smooth_j(b2n_temp)

            ## is it < 2 or <=2? only < 2 works to match the output of the original branch. ##
            if atype < 2:
                simple_smooth_j(b1p_temp)
                simple_smooth_j(b2p_temp)

    with open(my_pickle, 'wb') as f:
        pickle.dump((b1p, b2p, b1a, b2a, b1n, b2n, prob_all), f)

    return (b1p, b2p, b1a, b2a, b1n, b2n, prob_all)

def get_dmsp_smooth(b1p, b2p, b1a, b2a, b1n, b2n, prob_all, dFdt, day_of_year):

    je_a = np.zeros((4,nmlt,nmlat))     #returns four kinds of aurora
    jn_a = np.zeros((4,nmlt,nmlat))
    
    #atype=0 for diffuse,  atype=1=mono,  2=wave,  3=ions (in je_a and jn_a)
    #some routines, especially plotting, need to know whether a passed
    #array is energy flux or number, e- or ions.  jtype distinguishes:
    #jtype=1 for electron energy flux; 2=ion energy flux; 3=e- number flux
    #jtype=4=ion number flux;  5=e- average energy  6=ion average energy

    # loading finished

    for e_or_n in range(2):     #energy then number flux
        for atype in range(4):

            if (atype <= 2) and (e_or_n == 0):
                jtype = 1
            if (atype == 3) and (e_or_n == 0):
                jtype = 2
            if (atype <= 2) and (e_or_n == 1):
                jtype = 3
            if (atype == 3) and (e_or_n == 1):
                jtype = 4

            je_all = np.zeros((4,nmlt,nmlat))        #all refers to seasons

            # same code as 2010 version to calculate je from coefficients
            # running linear regressions in premodel folder to calculate aurora probabilities
            for iseason in range(4):
                #~ print('iseason',iseason)
                prob_curr = np.ones((nmlt, nmlat))

                if atype<=2:
                    for i in range(nmlt):
                        for j in range(nmlat):
                            b1t = b1p[ iseason, atype, i,j]
                            b2t = b2p[ iseason, atype, i,j]
                            # prob_estimate -> function in utilities_functions_season_epoch.py
                            prob_curr[i,j] = prob_estimate(b1t,b2t,dFdt,atype,i,j,prob_all[iseason])

                # energy fluxes at nmlt, nmlat
                je = np.zeros((nmlt,nmlat))

                # populating energy fluxes array

#               if jtype <=2: je = dFdt*b2a[iseason,atype,:,:]+b1a[iseason,atype,:,:]*prob_curr
#               if jtype >=3: je = dFdt*b2n[iseason,atype,:,:]+b1n[iseason,atype,:,:]*prob_curr

#               np.where(je <0.0, je, 0) # this does not do anything

#               je1 = spike_removal(je)
#               plt.imshow(je1)

                if jtype <= 2:
                    je = np.where(prob_curr > 0., (dFdt*b2a[iseason, atype] + b1a[iseason, atype]) * prob_curr, 0.)
                    je = np.where(je < 0., 0, je)

                    if atype<=2 and jtype==1:
                        je = np.where(je > 10., 0.5, je)
                        je = np.where(je > 5., 5., je)

                    if atype==3:
                        je = np.where(je > 6., 0.25, je)
                        je = np.where(je > 4., 4, je)

                else: # jtype == 3 or jtype == 4
                    je = np.where(prob_curr > 0., (dFdt*b2n[iseason, atype] + b1n[iseason, atype]) * prob_curr, 0.)

                    if atype<=2:
                        je = np.where(je < 0., 0, je)
                        je = np.where(je > 4.e10, 0, je)
                        je = np.where(je > 4.e9, 2.e9, je)

                    if atype==3:
                        je = np.where(je < 0., 0, je)
                        je = np.where(je > 5.e8, 0, je)
                        je = np.where(je > 1.e8, 1.e8, je)

                je_all[iseason,:,:] = je[:,:]
                # end of loop over iseason

            w0, w1, w2, w3 = season_weights(day_of_year)
#            for i in range(nmlt):
#                for j in range(nmlat):
            je =w0*je_all[0] + w1*je_all[1] + w2*je_all[2] + w3*je_all[3]

            if e_or_n == 0:
                je_a[atype,:,:] = je[:,:]
            else:
                jn_a[atype,:,:] = je[:,:]
            # end of aurora type loop
            # end of number or energy flux loop
#   plt.subplot(2,2,1)
#   plt.imshow(je_a[0,:,:])
#   plt.subplot(2,2,2)
#   plt.imshow(je_a[1,:,:])
#   plt.subplot(2,2,3)
#   plt.imshow(je_a[2,:,:])
#   plt.subplot(2,2,4)
#   plt.imshow(je_a[3,:,:])
#   plt.show()



    return je_a, jn_a

#get_guvi_dfdt.py:  Return an energy flux array, jea, for a
#specified dFdt -- with seasonal variations incorporated
#Uses GUVI data from the satellite TIMED
#Geolocation, dayglow removal, conversion to physical units by Kan Liou
#Model construction by Patrick Newell, November 2012
#The returned jea is calibrated to DMSP.
#Calibration to DMSP is based on comparisons at moderate activity,
#namely dFdt = 4300, especially premidnight (where both have data)
#to change calibration, change "convert_guvi" factor
#
#Based of IDL code from Patrick Newell, May 2013

def get_guvi_data( guvi_path ):
    my_pickle = '{}/guvi.pkl'.format(guvi_path)

    try:
        with open(my_pickle, 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        pass

    #just one set of regression coefficients, each nmlt x nmlat
    convert_guvi = 1.25
    dmsp_calibrate = True   #apply GUVI to DMSP conversions (MLT dependent)

    sname_a = ['winter','spring','summer','fall']

    b1a = np.zeros((4,nmlt,nmlat))           #4 for four seasons
    b2a = np.zeros((4,nmlt,nmlat))

    for iseason in range(4):
        sname = sname_a[iseason]
        afile = guvi_path + sname + '_totl.txt'

        # now read in regression coefficients for auroral flux
        with open( afile, 'r') as f:
        # reading just the first line in the file
            fline=f.readline().rstrip()
            y0_guvi, d0_guvi, yend_guvi, dend_guvi, files_guvi, sf0 = np.array(fline.split())
            i = 0
            j = 0
            b1 = 0.
            b2 = 0.
            rFa = np.zeros((nmlt,nmlat))
            while True:
                l=f.readline().rstrip()
                #Prob[0]print 'line', l, np.array(l.split())
                if l == '':
                    # either end of file or just a blank line.....
                    # we'll assume EOF
                    break
                i,j,b1,b2,rF = np.array(l.split())
                i,j = int(i), int(j)
                rFa[i,j] = rF
                b1a[iseason,i,j] = b1
                b2a[iseason,i,j] = b2


        # smoothing step here with simple_smooth_j
        smooth = 0
        if smooth==1:
            b1_temp = np.zeros((nmlt,nmlat))
            b2_temp = np.zeros((nmlt,nmlat))
            b1_temp[:,:] = b1a[iseason,:,:]
            b2_temp[:,:] = b2a[iseason,:,:]
            b1_smooth = simple_smooth_j(b1_temp)
            b2_smooth = simple_smooth_j(b2_temp)
            b1a[iseason,:,:] = b1_smooth[:,:]
            b2a[iseason,:,:] = b2_smooth[:,:]

    if dmsp_calibrate:
        b1a = b1a*convert_guvi
        b2a = b2a*convert_guvi

    with open(my_pickle, 'wb') as f:
        pickle.dump((b1a, b2a), f)

def get_guvi_dfdt( b1a, b2a, dFdt, day_of_year ):
    print ('Entering get_guvi_dfdt')
    
    # grid size of model (fixed)
    je_all = np.zeros((4,nmlt,nmlat))

    #~ smooth = 1            #apply smoothing before using coefficients
    smooth = 0          #smoothing can still be applied to final product
    #dmsp_calibrate = 0 #don't convert to DMSP values

    
    for iseason in range(4):
        #~ je=np.zeros((nmlt,nmlat))
        #~ for i in range(nmlt):
            #~ for j in range(nmlat):
                #~ je[i,j] = dFdt*b2a[iseason,i,j]+ b1a[iseason,i,j]
                #~ if je[i,j]<0.:
                    #~ je[i,j] = 0.
                #~ if je[i,j]>25.:
                    #~ je[i,j] = 25.
                #~ if je[i,j]>15.:
                    #~ je[i,j] = 15.
        #~ je_all[iseason,:,:] = je[:,:]
        je_all[iseason,:,:] = dFdt*b2a[iseason,:,:]+ b1a[iseason,:,:]
        
        je_all[iseason,:,0:hmlat] = je_all[iseason,:,hmlat:hmlat+hmlat]

    
    #~ plt.subplot(2,2,1)
    #~ plt.imshow(je_all[0,:,:])
    #~ plt.subplot(2,2,2)
    #~ plt.imshow(je_all[1,:,:])
    #~ plt.subplot(2,2,3)
    #~ plt.imshow(je_all[2,:,:])
    #~ plt.subplot(2,2,4)
    #~ plt.imshow(je_all[3,:,:])
    #~ plt.show()
    #~ quit()

    sday_of_year = day_of_year + 182
    if sday_of_year > 365: sday_of_year = sday_of_year - 365

    #~ Adjust the northern hemisphere for season    
    w0, w1, w2, w3 = season_weights(day_of_year)
    jea = w0*je_all[0,:,:] + w1*je_all[1,:,:] + w2*je_all[2,:,:] + w3*je_all[3,:,:] # nmlt x nmlat array
    
    #~ Adjust the southern hemisphere for season
    w0, w1, w2, w3 = season_weights(sday_of_year)
    jea[:,hmlat:2*hmlat] = w0*je_all[0,:,hmlat:2*hmlat] + w1*je_all[1,:,hmlat:2*hmlat] + w2*je_all[2,:,hmlat:2*hmlat] + w3*je_all[3,:,hmlat:2*hmlat]

    

    return jea


# combined_get.py:  Returns energy and number flux appropriate to a
# specified dF/dt with DMSP supplemented by GUVI
# seasonal variations incorporated
# Below threshold (dF/dt <= 12000.) DMSP data is used
# Above threshold, energy flux is from GUVI (previously normalized to
# DMSP) but average energy is calculated from DMSP (using the ratio of
# energy flux to number flux for each type of aurora).  Number flux is
# then taken from GUVI energy flux divided by the average energy for
# each auroral type.
# The partition of GUVI energy among auroral types is also based on the
# DMSP ratios
#
# Based on IDL code by Patrick Newell, January 2013 and Updated May 2013

def combined_get(dmsp_data, guvi_data, dFdt, day_of_year):
#   print ('Entering combined_get')

    #call supplies dFdt, day_of_year.  je_a, jn_a are returned
    #je_a = np.zeros(atype,nmlt,nmlat)
    #atype=0 is diffuse, 1=mono, 2=wave, 3=ions
    #for southern hemisphere, call with:
    #south_day = ((real_day + 365/2) mod 366) where
    #real_day = actual day of the year

    nmlt = 96                        #96 mlt bins
    nmlat = 160                    #80 mlat bins 50-90 x 2 hemispheres
    dFdt_change = 12000.              #change to GUVI above this
    je_a = np.zeros((4,nmlt,nmlat))   #energy flux with 4 types of aurora
    jn_a = np.zeros((4,nmlt,nmlat))   #number flux

    # get dmsp data
    jea_dmsp,jna_dmsp = get_dmsp_smooth( *dmsp_data, dFdt, day_of_year )

    # if below treshold, only use DMSP data
#   print ('Solar Wind Energy Transfer Function dFdt = ', int(dFdt))
    #dFdt = 11000.
    if dFdt < dFdt_change:
        je_a = jea_dmsp
        jn_a = jna_dmsp
    
        #~ print ('(1) in combined_get jn_a 0 ',jn_a[0,:,:])
    
        # above dFdt_change, need to obtain GUVI data
    else:
        # get guvi data
        jea_guvi = get_guvi_dfdt(  *guvi_data, dFdt, day_of_year )

        #~ print (' jea_guvi  ',jea_guvi[:,80])
        #~ quit()
        
        #jea_guvi is np.zeros(nmlt,nmlat) (No individual auroral types)
        #to partition jea_guvi into four auroral types and to get number flux,
        #we need still need DMSP
        for i in range(nmlt):
            for j in range(nmlat):
                je_sum=jea_dmsp[0,i,j] + jea_dmsp[1,i,j] + jea_dmsp[2,i,j] + jea_dmsp[3,i,j]
                if je_sum>0:
                    f0 = jea_dmsp[0,i,j]/je_sum     #fraction due to diffuse aurora
                    f1 = jea_dmsp[1,i,j]/je_sum     #monenergetic fraction
                    f2 = jea_dmsp[2,i,j]/je_sum     #broadband fraction
                    f3 = jea_dmsp[3,i,j]/je_sum     #ion fraction
                    jna0 = jna_dmsp[0,i,j]
                    eave0 = 10000.          #typical noise value
                    if jna0 > 0:
                        eave0 = (jea_dmsp[0,i,j]/jna0)/1.602e-12
                        # DMSP measures 32 to 30keV. Outside this range is a quirk
                        if eave0 < 32.:
                            eave0 = 32.
                        if eave0 > 30000.:
                            eave0 = 30000.
                        
                    jna1 = jna_dmsp[1,i,j]
                    eave1 = 10000.
                    if jna1 > 0:
                        eave1 = (jea_dmsp[1,i,j]/jna1)/1.602e-12
                        # DMSP measures 32 to 30keV. Outside this range is a quirk
                        if eave1 < 32.:
                            eave1 = 32.
                        if eave1 > 30000.:
                            eave1 = 30000.
                        
                    jna2 = jna_dmsp[2,i,j]
                    eave2 = 10000.
                    if jna2 > 0:
                        eave2 = (jea_dmsp[2,i,j]/jna2)/1.602e-12
                        #DMSP measures 32 to 30keV. Outside this range is a quirk
                        if eave2 < 32.:
                            eave2 = 32.
                        if eave2 > 30000.:
                            eave2 = 30000.

                    jna3 = jna_dmsp[3,i,j]
                    eave3 = 10000.
                    if jna3 > 0:
                        eave3 = (jea_dmsp[3,i,j]/jna3)/1.602e-12
                        #DMSP measures 32 to 30keV. Outside this range is a quirk
                        if eave3 < 32.:
                            eave3 = 32.
                        if eave3 > 30000.:
                            eave3 = 30000.

                    je_a[0,i,j] = f0*jea_guvi[i,j]
                    je_a[1,i,j] = f1*jea_guvi[i,j]
                    je_a[2,i,j] = f2*jea_guvi[i,j]
                    je_a[3,i,j] = f3*jea_guvi[i,j]
                    jn_a[0,i,j] = je_a[0,i,j]/eave0
                    jn_a[1,i,j] = je_a[1,i,j]/eave1
                    jn_a[2,i,j] = je_a[2,i,j]/eave2
                    jn_a[3,i,j] = je_a[3,i,j]/eave3


    #the conversion of ave energy to eV was for convenience since those
    #values are quite familar.  But energy flux is really ergs/cm2 s so we
    #have to remove that factor
    
        #~ print ('(2) in combined_get jn_a 0 ',jn_a[0,:,:])
    
        jn_a = jn_a/1.602e-12     #number flux in electrons/cm2-s
    
    #~ print ('(3) in combined_get jn_a 0 ',jn_a[0,:,:])
    #~ quit()
    
    return je_a, jn_a

# season_epoch function
# inputs:
#
#   aname : auroral name ('diff','mono','wave','ions')
#   atype : auroral type ( 0, 1, 2, 3 )
#   jtype : output types ( 1, 1, 1, 2 ); 1 = electron energy flux; 2 = ion energy flux
#   urlpath: path to realtime solar wind data
#   mode: 'FORECAST' or 'NOWCAST'
#
# outputs:
#
#   energy flux: je_all

def season_epoch_Kp(dmsp_data, guvi_data, sw_avg, mode, day_of_year,Kp_1):
    
    nmlt = 96  #number of mag local times in arrays
    nmlat = 160  #number of mag latitudes in arrays
    

    # check if real-time data was succesfully extracted


        # solar wind coupling function

    if mode == "FORECAST":      
    #    ********************    Convert from KP to Solar Wind Transfer Value  ********
        dFdt =  800.  -  224. * Kp_1    +  1154. * Kp_1**2  -   124. * Kp_1**3 + 20. *Kp_1**4
    else:

        if len(sw_avg) > 0:
            Bmag,Bx,By,Bz,v,ni,Ec=ap_inter_sol_realtime(sw_avg)  #calculating interpolated solar wind hourly averages and coupling founction
#           Bmag = 0

        if len(sw_avg) == 0 or Bmag == 0:
            print ('No real-time data available...  Using real-time Kp to run the model')
            time_now = datetime.datetime.utcnow() # Set Current Time
            
            time_Kp, Kp_1 , time_issue = swpc_get_Kp_NOWCAST_data_json()    
            sw_avg = { 'current_time' : time_now, 'time_latest_solar_wind' : time_issue,  'forecast_time': time_Kp, 'Bx' : 0, 'By' : 0, 'Bz' : 0, 'B_average' : 0, 'v' : 0, 'ni' : 0 }          
            dFdt =  800.  -  224. * Kp_1    +  1154. * Kp_1**2  -   124. * Kp_1**3 + 20. *Kp_1**4
        else:
            dFdt = Ec[19] # the 20th value of the coupling function is used

    # check for magnetic field value


    je_a, jn_a = combined_get(dmsp_data, guvi_data, dFdt,day_of_year)
    

    #~ print('je_a  ', je_a)
    

    #~ print('jn_a  ', jn_a)
    
    #~ quit()

    # returns auroral fluxes je_all
    #~ return Bx, By, Bz, v, ni, je_a
    return je_a


