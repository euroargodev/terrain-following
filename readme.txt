This is a readme file for the terrain-following interpolation algorithm 
documented in Yamazaki+(2020) https://doi.org/10.1029/2019JC015406.

The codes are originally designed to be applied for Argo data in *_prof.nc format fetched from Argo Global Data Assembly Center.

interpolation_v2.py
    ...the interpolation algorithm with some plotting function. This is slightly modified from Yamazaki+(2020) for clarity of the scheme and customizability.

Test trajectory data (testData.csv) is included in this folder, which can be directly loaded to interpolation_v2.py.
The result of interpolation for test data is shown in two .png files.

Following codes would help pre-processing of the trajectory data:

nc2csv.py
    ...extract trajectory data and position quality flag (qpos) from netcdf file in the Argo profile format (*_prof.nc).

process.py
    ...change id for data seeming sequential as trajectory by id-cast but problematic in lon-lat coords.	

Correspondence: Kaihe Yamazaki (kaiheyamazaki@gmail.com)
written in 2019 September, modified in 2021 July.