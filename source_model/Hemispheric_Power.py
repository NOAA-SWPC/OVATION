# Hemispheric_Power.py
# Calculates Hemispheric Power

# Rodney Viereck (2018)

import numpy as np
#import matplotlib.pyplot as plt

def Hemispheric_Power(NS, mlt_bin, mlat_bin, je):
    
    #~ print ('1 ',je[0,10:20,100:110])
        
        
    nmlt = 96  #number of mag local times in arrays
    nmlat = 160  #number of mag latitudes in arrays
    ndF = 12

    mlt_array = []
    mlat_array = []
    je_array = []
    je_d = []
    je_m = []
    je_w = []
    je_i = []
    je_s = []

    power_hemi = 0
    
        # combining all 4 auroral types
    je_combined = je[0,:,:] + je[1,:,:] + je[2,:,:]  + je[3,:,:]  
    je_sum =  je[0,:,:] + je[1,:,:] + je[2,:,:] #   Electrons Only

    ############################
    # from draw_je.pro, line 325
    # combining north and south hemisphere
    # or extract just northern/southern hemisphere
    # to obtain average auroral flux

    j = 0
    
    for i in range(nmlt):
        mlt = 0.001 + 0.25*float(i)
        mlt_real = mlt
        mlt = mlt - 6.0 # rotate noon to top
        if mlt < 0.:
            mlt = mlt + 24.
        if mlt > 24.:
            print ('odd mlt (mlt,mlat)=',mlt,' i,j = ',i,j)
        for j in range(int(nmlat/2)):
            mlat = 50.001 + 0.5*float(j)
            mlat_real = mlat
            ibin = mlt_bin(mlt_real)
            jbin = mlat_bin(mlat_real)

            # flux array
            jeave = je_combined[ibin,jbin] # average je

            # combining north and south hemisphere to smooth noisier south
            jbin2 = mlat_bin(-mlat_real)
            je1 = jeave
            je2 = je_combined[ibin,jbin2]
            jeave = (je1 + je2)/2.
            if je1*je2 == 0.:
                jeave = je1 + je2

            # calculating hemispheric power
            mlat1 = mlat
            mlat2 = mlat + 0.5
            mlt1 = mlt_real
            mlt2 = mlt1 + 0.25
            # area described by solid angle on Earth returned in km2
            Re = 6371.      # Earth radius in km
            delta_phi = (2.*np.pi/24.)*np.absolute(mlt2-mlt1)
            theta1 = (90. - mlat1)*np.pi/180.# in radians
            theta2 = (90. - mlat2)*np.pi/180.
            A_earth = delta_phi*Re*Re*np.absolute(np.cos(theta1)-np.cos(theta2))

            power_hemi = power_hemi + 1.e-6*A_earth*jeave

            # combining coordinates and fluxes into 1D arrays: mlt, mlat, je
            #~ print(type(mlt_array), type(mlt_real))

            if NS!=0:
                mlat_real = -mlat_real
            mlt_array.append(mlt_real)
            mlat_array.append(mlat_real)
            je_array.append(jeave)
            je_d.append(je[0,ibin,jbin2])
            je_m.append(je[1,ibin,jbin2])
            je_w.append(je[2,ibin,jbin2])
            je_i.append(je[3,ibin,jbin2])
#           je_s.append(je_sum[ibin,jbin2])
            #~ print (ibin,jbin2,je[0,ibin,jbin2])
            #~ print (je_d[-1])
        
    #~ plt.imshow(je[0,:,:])
    #~ plt.plot(je_d)
    #~ plt.show()
    #~ quit()

    # turning data into numpy arrays
    mlt_array = np.asarray( mlt_array )
    mlat_array = np.asarray( mlat_array )
    je_array = np.asarray( je_array )
    je_d = np.asarray(je_d)
    je_m = np.asarray(je_m)
    je_w = np.asarray(je_w)
    je_i = np.asarray(je_i)
#   je_s = np.asarray(je_s)
    

    return (mlt_array, mlat_array, je_array, je_d, je_m, je_w, je_i, power_hemi)



