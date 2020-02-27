# solar_wind_coupling.py calculates the solar wind coupling function
# using magnetic field, speed and density as inputs
#
# input parameters:
#   Bx: in, required, type='numpy.ndarray'
#   By: in, required, type='numpy.ndarray'
#   Bz: in, required, type='numpy.ndarray'
#   v: in, required, type='numpy.ndarray'
#   ni: in, required, type='numpy.ndarray'
#
# returns:
#   Ec: coupling function
#
# Author: Diana Morosan
# morosand@tcd.ie
#
# Based on IDL code by Patrick Newell, 2009 and
# Ovation Model Operational IDL code


import numpy as np

# boro.py function to evaluate Borovsky's function
# based on boro.pro by Patrick Newell, March 2008 from
# Ovation Model Operational IDL code
def boro(sintc,ni,v,B):
    if ni<0:
        Eb = 0
        return Eb
    
    mp = 1.
    rhom = 0.
    rho = ni*mp
    VA = (2.18e6)*B/np.sqrt(ni) # in cm/s
    VA = VA/(1.0e5) # in km/s
    MA = v/VA
    MMS = MA
    C = ((0.25)**6. + (1./(1. + 1.38*np.log(MA)))**6.)**(-0.166667)
    betaS = 0.032*(MA**1.92)
    fact1 = 1.6*sintc
    fact2 = rho*v*v
    fact3 = (1 + 0.5/(MMS*MMS))
    fact4 = 1./np.sqrt(1. + betaS)
    fact5 = (C*rho + (1./np.sqrt(1.+betaS))*rhom)**(-0.5)
    fact6 = (np.sqrt(1.+betaS)+1.)**(-0.5)
    R = fact1*fact2*fact3*fact4*fact5*fact6
    Eb = R
    
    return Eb


# solar wind coupling function
def sol_coup(Bx, By, Bz, v, ni):
    ncoup = 33
    Ec = np.zeros(ncoup)
    B = np.sqrt(Bx**2 + By**2 + Bz**2)
    
    #if ni*B*v==0:
    #   return
    
    BT = np.sqrt(By**2 + Bz**2)
    Bztemp = Bz
    if Bz == 0:
        Bztemp = 0.001
    
    tc = np.arctan2(By,Bztemp)        #tc = np.arctan(By,Bztemp) -- sol_coup.pro       
    if BT*np.cos(tc)*Bz < 0:
        tc = tc + 3.14159265
    
    sintc = np.absolute(np.sin(tc/2.))
    Bs = Bz
    if Bs > 0:
        Bs = 0.
    EKL = v*BT*sintc*sintc
    Ewav = EKL*sintc*sintc

    Ec[0] = Bz
    Ec[1] = EKL
    eps = v*B*B*sintc*sintc*sintc*sintc
    #Ec(2) = eps
    Ec[2] = v*np.sqrt(BT)*sintc*sintc       #Lyatsky
    Ec[3] = v*Bs
    p = ni*v*v*(1.67e-6)
    Ec[4] = p                         #dynamic pressure in nPa
    Ec[5] =sintc*sintc*EKL*(v**0.3333)*(p**0.16667)       #Vasylinas
    Ec[6] =sintc*sintc*BT*(p**0.33333)        #Another Vasyliunas formula
    Ec[7] = Ewav
    #Ec[8] = ni
    Ec[8] = np.sqrt(ni)*v*v          #8/2007 change for twovar paper plot
    Ec[9] = v
    Ec[10] = Bs
    Ec[11] = v*BT
    Ec[12] = v*BT*BT*sintc*sintc*sintc*sintc
    Ec[13] = v*B*sintc*sintc*sintc*sintc
    Ec[14] = Ewav*Ewav
    Ec[15] = np.sqrt(Ewav)
    Ec[16] = Ewav**0.6667
    Ec[17] = np.sqrt(EKL)
    Ec[18] = EKL*(p**0.166667)*(v**0.3333)
    Ec[19] = (v**1.33333)*(sintc**2.66667)*(BT**0.66667)
    #Ec[19] = (v^1.33333)*exptc*(BT^0.66667)
    Ec[20]  = Ec[19]/(v**0.33333)
    Ec[21] = (v**0.166667)*Ec[19]
    Ec[22] = (sintc**0.333333)*Ec[19]
    Ec[23] = 0.
    if sintc!=0.:
        Ec[23] = Ec[19]/(sintc**0.66667)
    Ec[24] = (BT**0.3333)*Ec[19]
    Ec[25] = 0.
    if BT > 0.:
        Ec[25] = Ec[19]/(BT**0.166667)
    Ec[26] = Ewav*np.sqrt(p)
    Ec[27] = v*(BT**0.67)*sintc*sintc      #Weimer
    Ec[28] = np.sqrt(ni)*v*v*BT*sintc*sintc*sintc*sintc*sintc*sintc
    Ec[29] = Ec[19]*np.sqrt(p)
    Ec[30] = Ec[19]*(ni**0.16667)
    #Ec[31] = sintc*(n**(-0.22))*(v**0.56)*(B**1.44)       #Borovosky
    Eb = boro(sintc,ni,v,B) # Borovosky function
    Ec[31] = Eb
    Ec[32] = Ec[31]*(sintc**1.667)
    
    return Ec

