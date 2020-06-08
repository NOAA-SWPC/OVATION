#  Getting Data from the swdsst database
#Collects Data from the Mag table and the Plasma table to create an output of
# Date     Time      Bx      By     Bz     Bt     Speed     Dens      SatPos_x     SatPos_y

import configparser
import pyodbc 
import datetime as dt
import numpy as np
from datetime import datetime, date, time



#**********    Read the Config File Input  *************************

config_file = 'C:\Docs\Python\GOES_Stand_Bands\Py_Source\DB_Config\Data_Base_Config.cfg'

config = configparser.ConfigParser()
config.read(config_file)

databaseServer = config.get('Database', 'database-server')
database = config.get('Database', 'database')

opath = "C:/Docs/Python/Ovation_New_Realtime_2019/SW_Data/"
output_file = "sw_data_2019_Sept.dat"

ofile = open(opath + output_file, "w")

ofile.write('Date     Time      Bx      By     Bz     Bt     Speed     Dens      SatPos_x     SatPos_y \n')



#  ******************  Set up database  ****************

conn_str = (
    r'DRIVER={SQL Server};'
    r'SERVER='+databaseServer+';'
    r'DATABASE='+database+';'
    r'Trusted_Connection=yes;'
)

cnxn = pyodbc.connect(conn_str)
cursor = cnxn.cursor()

#~ # ************    CreateDatabase Call   ***********************

Ys = 2019
Ms = 9
Ds = 1
Hs = 0
Mins = 0

Ye = 2019
Me = 10
De = 1
He = 0
Mine = 0

SD = dt.datetime(Ys, Ms, Ds, Hs, Mins)
ED = dt.datetime(Ye, Me, De, He, Mine)

start_time = ("'"+SD.strftime("%Y-%m-%d %H:%M:%S")+"'")
end_time = ("'"+ED.strftime("%Y-%m-%d %H:%M:%S")+"'")


#  ************   Transition Time   *********************

TD = dt.datetime(2016,3,9,15,27)

trans_time =  ("'"+TD.strftime("%Y-%m-%d %H:%M:%S")+"'")



meq1 = ("SELECT [time_tag],[gse_bx],[gse_by],[gse_bz],[bt],[gsm_lat],[gsm_lon] FROM [dbo].[tb_ace_mag_1m] ")
meq2 = (" WHERE [time_tag] BETWEEN ")
meq3 = ("%s AND %s ORDER BY [time_tag]" % (start_time , end_time))

mlq1 = ("SELECT [time_tag],[bx_gse],[by_gsm],[bz_gsm],[bt] FROM [rtsw].[all_magnetic_field] ")
mlq2 = (" WHERE [active]=1 AND [cadence]=60 AND [time_tag] BETWEEN ")
mlq3 = ("%s AND %s ORDER BY [time_tag]" % (start_time , end_time))

peq1 = ("SELECT [time_tag],[speed], [dens] FROM [dbo].[tb_ace_sw_1m] ")
peq2 = (" WHERE [time_tag] BETWEEN ")
peq3 = ("%s AND %s ORDER BY [time_tag]" % (start_time , end_time))

plq1 = ("SELECT [time_tag],[proton_speed], [proton_density] FROM [RA].[rtsw].[all_solar_wind] ")
plq2 = (" WHERE [active]=1 AND [time_tag] BETWEEN ")
plq3 = ("%s AND %s ORDER BY [time_tag]" % (start_time , end_time))


if SD > TD:
	querym = mlq1+mlq2+mlq3
	queryp = plq1+plq2+plq3	
	
if ED < TD:
	querym = meq1+meq2+meq3
	queryp = peq1+peq2+peq3
	
	
#  ****************   Get Mag Data   *************************

print (querym)

cursor.execute(querym)

row = cursor.fetchone()

mdate = []
Bx = []
By = []
Bz = []
Btot = []
Glat = []
Glon = []

while row:
	mdate.append(row[0])
	Bx = np.append(Bx,row[1])
	By = np.append(By,row[2])
	Bz = np.append(Bz,row[3])
	Btot = np.append(Btot,row[4])
	#~ Glat = np.append(Glat,row[5])
	#~ Glon = np.append(Glon,row[6])

	row = cursor.fetchone()

#  ****************   Get Solar Wind Data   *************************

print (queryp)

cursor.execute(queryp)

row = cursor.fetchone()

pdate = []
V = []
N = []


while row:
	if row[1]:
		pdate.append(row[0])
		V = np.append(V,row[1])
		N = np.append(N,row[2])

	row = cursor.fetchone()


#*********************************************************************

print (len(mdate), len(pdate))

im= 0
ip = 0

for i1 in range(0,len(mdate)-1):

	if mdate[im] == pdate[ip]:
		mtime = mdate[im].strftime("%Y-%m-%d %H:%M:%S")
		ptime = pdate[ip].strftime("%Y-%m-%d %H:%M:%S") 
		#~ print (' %20s %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e   %8.4e \n' %(mtime, Bx[im], By[im], Bz[im], Btot[im], V[ip], N[ip], Glat[im], Glon[im]))
		#~ ofile.write(' %20s %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e   %8.4e \n' %(mtime, Bx[im], By[im], Bz[im], Btot[im], V[ip], N[ip], Glat[im], Glon[im]))
		if N[ip]: ofile.write(' %20s %8.4e  %8.4e  %8.4e  %8.4e  %8.4e   %8.4e \n' %(mtime, Bx[im], By[im], Bz[im], Btot[im], V[ip], N[ip]))
		im = im+1
		ip = ip+1
	else:
		mtime = mdate[im].strftime("%Y-%m-%d %H:%M:%S")
		ptime = pdate[ip].strftime("%Y-%m-%d %H:%M:%S") 
		if mdate[im] < pdate[ip]:
			im=im+1
		elif mdate[im]> pdate[ip]:
			ip = ip+1




ofile.close()


exit()
