### This code extracts the trajectry of Argo data
### and divides position-interpolated data for Argo netcdf

# coding: utf-8
import numpy as np
import glob
import netCDF4
import re
import pandas as pd
import datetime

# path = 'B:/Data/Argo/Argo_60S30-160E'
names = glob.glob(path+'/ALL/*/*_prof.nc')
names.sort()
n_float = len(names)

dfs = []
for Float in range(n_float):
    # if Float%20 == 0:
    #     print(Float, 'th data of ', n_float)
    Prof = netCDF4.Dataset(names[Float])
    idnum = str(Prof.variables['PLATFORM_NUMBER'][0])
    idnum = re.sub(r'\D', '', idnum) # '\D' means strings except number
    # print(idnum)
    cast = Prof.variables['CYCLE_NUMBER'][:]
    # if len(cast) < 2:
    #     return
    juld = Prof.variables['JULD'][:]
    # juld.units = 'days since 1950-01-01 00:00:00 UTC'
    # juld = netCDF4.num2date(juld[:],juld.units)
    lons = Prof.variables['LONGITUDE'][:]
    lats = Prof.variables['LATITUDE'][:]
    qpos = str(Prof.variables['POSITION_QC'][:]) # positioning quality flag
    qpos = re.sub(r'\D', '', qpos) # '\D' means strings except number
    ids = np.array([idnum]*len(cast))
    df = pd.DataFrame([ids,cast,qpos,lons,lats,juld]).T
    df.columns = ['id','cast','qpos','lons','lats','juld']
    dfs.append(df)

df_ = pd.concat(dfs, ignore_index=True)
df_.to_csv(path+'/TrajectoryCSV/trajAll.csv')
df_.describe()
