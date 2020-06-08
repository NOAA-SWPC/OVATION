################################################
# Converts magnetic local time and magnetic
# latitude to geographic longitude and latitude
# using Dipole Approximation.
# Based on ov_cgm2geo.pro (IDL OVATION 2010)
#
#Author: Diana Morosan
#   morosand@tcd.ie
#

import numpy as np

def mag2geo(time_sw,mlt_array,mlat_array):

    ut_hour = time_sw.hour + time_sw.minute/60. + time_sw.second/3600

    # dipole approximation
    prime_magnetic_meridian = 71.
    # converting magnetic local time to magnetic longitude
    mag_lon = ( mlt_array - ut_hour )*15. + prime_magnetic_meridian

    idx_lt_neg180 = np.where( mag_lon < -180 )
    idx_ge_pos180 = np.where( mag_lon >  180 )

    if len(idx_lt_neg180) > 0:
        mag_lon[ idx_lt_neg180 ] = mag_lon[ idx_lt_neg180 ] + 360.
    if len(idx_ge_pos180) > 0:
        mag_lon[ idx_ge_pos180 ] = mag_lon[ idx_ge_pos180 ] - 360.

    #
    # transforming cgm to geo using the dipole assumption
    # based on mag2geo.pro
    #
    #
    # Some 'constants':
    Dlong=288.59   # longitude (in degrees) of Earth's magnetic south pole
    # (which is near the geographic north pole!) (1995)
    Dlat=79.30     # latitude (in degrees) of same (1995)
    R = 1.         # distance from planet center (value unimportant --
    # just need a length for conversion to rectangular coordinates)

    #convert first to radians
    Dlong=Dlong*np.pi/180.
    Dlat=Dlat*np.pi/180.

    mlat=mlat_array*np.pi/180.
    mlon=mag_lon*np.pi/180.
    malt=mlat * 0. + R

    coord=[mlat,mlon,malt]

    #convert to rectangular coordinates
    #     X-axis: defined by the vector going from Earth's center towards
    #             the intersection of the equator and Greenwich's meridian.
    #     Z-axis: axis of the geographic poles
    #     Y-axis: defined by Y=Z^X
    x=coord[2]*np.cos(coord[0])*np.cos(coord[1])
    y=coord[2]*np.cos(coord[0])*np.sin(coord[1])
    z=coord[2]*np.sin(coord[0])

    #First rotation matrix: in the plane of the current meridian from
    #magnetic pole to geographic pole.
    togeolat=np.zeros((3,3))
    togeolat[0,0]=np.cos(np.pi/2-Dlat)
    togeolat[0,2]=np.sin(np.pi/2-Dlat)
    togeolat[2,0]=-np.sin(np.pi/2-Dlat)
    togeolat[2,2]=np.cos(np.pi/2-Dlat)
    togeolat[1,1]=1.
    # Computes array elements by multiplying the columns
    # of the first array by the rows of the second array
    # equivalent of # in IDL
    out = np.dot( togeolat, np.asarray([x,y,z]) ).T


    #Second rotation matrix : rotation around plane of the equator, from
    #the meridian containing the magnetic poles to the Greenwich meridian.
    maglong2geolong=np.zeros((3,3))
    maglong2geolong[0,0]=np.cos(Dlong)
    maglong2geolong[0,1]=-np.sin(Dlong)
    maglong2geolong[1,0]=np.sin(Dlong)
    maglong2geolong[1,1]=np.cos(Dlong)
    maglong2geolong[2,2]=1.
    out=np.dot(maglong2geolong , out.T).T

    #convert back to latitude, longitude and altitude
    glat=np.arctan2(out[:,2],np.sqrt(out[:,0]**2+out[:,1]**2))
    glat=glat*180./np.pi # geographic latitude
    glon=np.arctan2(out[:,1],out[:,0])
    glon=glon*180./np.pi # geographic longitude

    return glon, glat
