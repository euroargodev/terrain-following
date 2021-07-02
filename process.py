### selection for good data(quality flag <= 2)
### if there is discontinuity, rename id

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path = 'B:/Data/Argo/TrajectoryCSV'
df = pd.read_csv(path + '/trajAll.csv').values[:,1:]
df[np.where(df=='--')] = 999; df = np.float64(df)
df[np.isnan(df)] = 999
df = pd.DataFrame(df)
df.columns = ['id','cast','qpos','lons','lats','juld']
### make lons [-180,180] --> [0,360]
df.lons[np.where(df.lons.values<0)[0]] += 360

ids = np.unique(df.id)
dfs = []
for idi in ids:
    dfi = df[df.id==idi].sort_values(by='cast')
    dfi = dfi.reset_index(drop=True)

    # selection for good data
    discont = []
    indx = np.where(dfi.qpos<=2)[0]
    for i in range(len(indx)-1):
        if indx[i+1] - indx[i] > 1: # search discontinuity for good data
            discont.append(i+1) # corresponds to i_bgn

    idiscont = 0
    for i in range(len(indx)-1):
        indxi = indx[i]
        if len(discont) == 0:
            ### if there is discontinuity, change id
            if (np.abs(dfi.lons[indxi+1]-dfi.lons[indxi]) > 10)\
              or (np.abs(dfi.lats[indxi+1]-dfi.lats[indxi]) > 5)\
              or (dfi.cast[indxi+1] - dfi.cast[indxi] > 1):
                dfi.id[indxi+1:] += 1e+7
        elif indxi == indx[discont[idiscont]-1]: # at discontinuity
            if idiscont < len(discont)-1:
                idiscont += 1
            continue
        else:
            ### if there is discontinuity, change id
            if (np.abs(dfi.lons[indxi+1]-dfi.lons[indxi]) > 10)\
              or (np.abs(dfi.lats[indxi+1]-dfi.lats[indxi]) > 5)\
              or (dfi.cast[indxi+1] - dfi.cast[indxi] > 1):
                dfi.id[indxi+1:] += 1e+7
    dfs.append(dfi)

df = pd.concat(dfs)
df = df.sort_values(by='id')
df = df.reset_index(drop=True)
df.to_csv(path + '/TrajectoryCSV/traj_processed.csv')
