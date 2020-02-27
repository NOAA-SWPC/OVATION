# Useful functions used throughout season_epoch.py
# Based on OVATION Model Operational
#
# function: dF_bin(dF_call)
#
# function: prob_estimate(b1,b2,dF,atype,i,j,Prob)
#
# function: season_weights(day_of_year)


##################################################
# used in prob_estimate
def dF_bin(dF_call):
    ndF = 12
    dFave = 4421.
    dFstep = dFave/8.
    dF = dF_call
    bin = int(dF/dFstep) # converts to integer
    if bin<0:
        bin = 0
    if bin>(ndF-1):
        bin = ndF-1
    return bin

##################################################
# probability estimate function
def prob_estimate(b1,b2,dF,atype,i,j,Prob):
    ndF = 12
    
    if atype == 2:
        j_n50 = mlat_bin(-50.)
        j_n60 = mlat_bin(-60.)
        j_50 = mlat_bin(50.)
        j_60 = mlat_bin(60.)

        #Zero out wave aurora (atype = 2) below 60 MLAT since statistics are so poor
        if j > j_n50 and j < j_n60:
            p = 0.
            return p
        if j > j_50 and j < j_60:
            p = 0.
            return p

    p = b1 + b2*dF
    if p<0.:
        p = 0.
    if p>1.:
        p = 1.
    
    if b1!=0. or b2!=0.:
        return p
    
    dFbin = dF_bin(dF)
    p0 = Prob[atype,dFbin,i,j]

    if p0!=0.:
        return p0
    
    
    dFbin_m1 = dFbin -1
    dFbin_p1 = dFbin + 1
    if dFbin_m1<0:
        dFbin_m1 = dFbin + 2

    if dFbin_p1>(ndF-1):
        dFbin_p1 = dFbin - 2
    
    p1 = (Prob[atype,dFbin_m1,i,j] + Prob[atype,dFbin_p1,i,j])/2.
    return p1


##################################################
# season weights
def season_weights(day_of_year):
    winter_w = 0.
    spring_w = 0.
    summer_w = 0.
    fall_w = 0.

    if day_of_year>=79 and day_of_year<171:
        summer_w = 1. - float(171-day_of_year)/92.
        spring_w = 1. - summer_w

    elif day_of_year>=171 and day_of_year<263:
        fall_w = 1. - float(263-day_of_year)/92.
        summer_w = 1. - fall_w

    elif  day_of_year>=263 and day_of_year<354:
        winter_w = 1. - float(354-day_of_year)/91.
        fall_w = 1. - winter_w
    else:
        # must be in the range 354 to 78. Convert 354-365 to negative numbers
        day_of_year0 = day_of_year
        if day_of_year>=354:
            day_of_year0 = day_of_year - 365

        spring_w = 1. - float(79-day_of_year0)/90.
        winter_w = 1. - spring_w

    return winter_w, spring_w, summer_w, fall_w

####################################################
# fuctions below are used to calculate the bin sizes
# in magnetic local time and latitude:

###################################################
def mlt_bin(mlt_call):
    mlt = mlt_call
    if mlt < 0.:
        mlt = mlt + 24.
    if mlt > 24.:
        mlt = mlt - 24.
    if mlt < 0. or mlt > 24.:
        return -1
 
    #bin = fix(mlt*2)
    bin = int(mlt*4)
    #if(bin gt 47) then bin =47
    if bin > 95:
        bin = 95
    if bin < 0:
        bin = 0

    return bin

###################################################
def mlat_bin(mlat_call):
    mlat=mlat_call
    
    if mlat < 0.:
        mlat = -mlat
        if mlat < 50. or mlat > 90.:
            return -1
        bin = int( (mlat-50.)*2 )
        if bin < 0.:
            bin = 0
        if bin > 79:
            bin = 79
    elif mlat > 0.:
        if mlat < 50. or mlat > 90.:
            return -1
        bin = int( (mlat-50.)*2 )
        if bin < 0.:
            bin = 0
        if bin > 79:
            bin = 79
        bin = bin + 80

    return bin


################################################
def mlt_plus(mlt_call, delta_t):
    nmlt = 96
    mlt_current = mlt_call

    if mlt_current >= nmlt:
        mlt_current = mlt_current - nmlt
    if mlt_current < 0:
        mlt_current = mlt_current + nmlt

    mlt_new = mlt_current + delta_t

    if mlt_new >= nmlt:
        mlt_new = mlt_new - nmlt
    if mlt_new < 0:
        mlt_new = mlt_new + nmlt

    return mlt_new

################################################
def mlat_plus(mlat_call, delta_lat):
    nmlat = 160
    ns_divide = nmlat/2
    old_lat = mlat_call
    if old_lat < ns_divide:
        new_lat = old_lat + delta_lat
        if new_lat < 0:
            new_lat=0
        if new_lat >= ns_divide:
            new_lat = ns_divide - 1
    else:
        new_lat = old_lat + delta_lat
        if new_lat<ns_divide:
            new_lat=ns_divide
        if new_lat>=nmlat:
            new_lat=nmlat-1

    return new_lat





