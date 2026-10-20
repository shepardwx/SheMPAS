import uxarray as ux
import glob
import pandas as pd
import numpy as np
import xarray as xr

" returns a uxds that originally had a 'rain_total' variable (1 hourly precipitation since init), with a new 'pcp{interval}' variable "
def rain_total_to_interval(uxds,interval):
    rain=uxds['rain_total']

    time_dim = rain.dims[0]
    later = rain.isel({time_dim: slice(interval, len(rain[time_dim]), interval)})
    earlier = rain.isel({time_dim: slice(0, len(rain[time_dim])-interval, interval)})
    pcp_interval_data = later.data - earlier.data
    other_dims = rain.dims[1:]
    pcp_interval_coords = {}

    for dim in other_dims:
        if dim in rain.coords and rain.coords[dim].dims == (dim,):
            pcp_interval_coords[dim] = rain.coords[dim]

    pcp_interval_coords[f"time{interval}"] = (
        f"time{interval}",
        uxds[time_dim].isel({time_dim: slice(interval, len(rain[time_dim]), interval)}).data
    )

    pcp = xr.DataArray(
        pcp_interval_data,
        dims=(f"time{interval}",) + other_dims,
        coords=pcp_interval_coords,
        name=f"pcp{interval}",
        attrs={
            "long_name": f"precipitation during the previous {interval} hours",
            "units": "inches",
        },
    )

    uxds[f"pcp{interval}"] = pcp

    return uxds


def hist_to_uxds(date,note,model,scope='regional',resolution='4',intervals=[6,24]):

    # FILE PATHS
    grid_path=f'/home/disk/hot/shepard0/mpas/meshes/{scope}/{resolution}/{note}/{note}.grid.nc'
    data_path_pfx=f'/home/disk/hot/shepard0/mpas/{scope}/{note}/{resolution}/{date}/{model}'

    # LIST OF HIST FILES
    hist_files=[]
    for file in glob.glob(f'{data_path_pfx}/'+'history*'):
        hist_files.append(file)

    hist_files=sorted(hist_files)

    # print(hist_files)

    # OPEN DATASET
    uxds=ux.open_mfdataset(
        grid_path,
        hist_files,
        combine='nested',
        concat_dim='time'
    )

    # FIX TIME
    times=[]
    for hist_file in hist_files:
        times.append(pd.to_datetime(hist_file[-22:-3], format='%Y-%m-%d_%H.%M.%S'))
    uxds=uxds.assign_coords(time=('time',times)).drop_vars('xtime')
    uxds=uxds.isel(Time=0)


    # GET rain_total
    if 'rainc' in uxds.data_vars:
        uxds['rain_total']=uxds['rainc']+uxds['rainnc']
    else:
        uxds['rain_total']=uxds['rainnc']

    # GET CREF
    if 'refl10cm' in uxds.data_vars:
        uxds['cref']=np.max(uxds['refl10cm'],axis=2)

    # INIT TIME, MODEL, RESOLUTION, NOTE, SCOPE
    uxds['init']=uxds['time'][0]
    uxds['model']=model
    uxds['resolution']=resolution
    uxds['note']=note
    uxds['scope']=scope

    # GET RAIN INTERVALS
    for interval in intervals:
        uxds=rain_total_to_interval(uxds,interval)

    # RETURN UXDS    
    return uxds


# FUNCTION TO GRAB ONLY LIST vars FROM A UXDS
def select_vars(uxds,vars=[]):
    return uxds[vars]

if __name__=='__main__':

    plotting_vars=['t2m',
                   'mslp',
                   'rain_total',
                   'cref',
                   'init',
                   'model',
                   'resolution',
                   'note',
                   'scope']

    print(
        select_vars(hist_to_uxds(
            date='2026081512',
            note='hawaii_shifted',
            model='RCON'
        ), vars=plotting_vars)
        )