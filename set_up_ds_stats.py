import hawaii_mesonet_downloader
from datetime import datetime
import uxarray as ux
import numpy as np
import xarray as xr
import hist_to_uxds
import latlon_verify

def set_up_ds(init_dates):
    # dates of MPAS initializations
    # init_dates=['2026031000','2026031012','2026031100','2026031112','2026031200','2026031212','2026031300','2026031312','2026031400','2026031412','2026031500','2026031512','2026031600']#[:11]
    init_ts_dim=[]
    for date in init_dates:
        init_ts_dim.append(np.datetime64(datetime.strptime(date, "%Y%m%d%H")))
    init_ts_dim=np.array(init_ts_dim)
    model=np.array(['RCON','Thompson'])
    ids=hawaii_mesonet_downloader.which_ids('all')
    station_id=[]
    for id in ids:
        station_id.append(f'0{id}')
    station_id=np.array(station_id)

    coords = {
        "init_ts": init_ts_dim,
        "station_id": station_id,
        "model": model,
    }

    dims = ("init_ts", "station_id", "model")

    variable_names = [
        # "bias_18_to_24",
        # "bias_00_to_24",
        'model_pcp_18_to_24',
        'model_pcp_00_to_24',
        'obs_pcp_18_to_24',
        'obs_pcp_00_to_24',
        'model_elevation',
        'obs_elevation',
        'lat',
        'lon'
    ]

    # Create a NaN-filled array with the full dimensional shape
    shape = tuple(len(coords[dim]) for dim in dims)
    nan_values = np.full(shape, np.nan, dtype=np.float64)

    ds_stats = xr.Dataset(
        data_vars={
            name: (dims, nan_values.copy())
            for name in variable_names
        },
        coords=coords,
    )

    ds_stats['model_elevation']=ds_stats['model_elevation'].isel(model=0,init_ts=0,drop=True)
    ds_stats['obs_elevation']=ds_stats['obs_elevation'].isel(model=0,init_ts=0,drop=True)
    ds_stats['lat']=ds_stats['lat'].isel(model=0,init_ts=0,drop=True)
    ds_stats['lon']=ds_stats['lon'].isel(model=0,init_ts=0,drop=True)

    return ds_stats




def fill_ds_stats(init_dates):

    ds_stats=set_up_ds(init_dates)

    grid_path='/home/disk/hot/shepard0/mpas/meshes/regional/4/hawaii_shifted/hawaii_shifted.grid.nc'
    init_path='/home/disk/hot/shepard0/mpas/regional/hawaii_shifted/4/2026031412/hawaii_shifted.init.nc'

    uxds_ter=ux.open_dataset(grid_path,init_path)

    for model in ds_stats['model'].values:
        for date_i,date in enumerate(ds_stats['init_ts'].values):
            #DATE IN STR FORM YYYYMMDDHH
            date_str=str(date).replace('-','').replace('T','')[:10]

            # MODEL DATA FOR THIS DATE
            # ds_model=latlon_verify.fix_time_return_ds(date=date_str,model=f'{model}_gf')

            uxds_model=hist_to_uxds.hist_to_uxds(date=date_str,
                        note='hawaii_shifted',
                        model=f'{model}_gf',
                        scope='regional',
                        resolution='4',
                        intervals=[6,24])

            # # OBS FOR THIS DATE
            end_date_str=str(uxds_model['time'].values[-1]).replace('-','').replace('T','')[:10]
            ds_5min=hawaii_mesonet_downloader.ds_5min(start_date=date_str,end_date=end_date_str)
            ds_obs=hawaii_mesonet_downloader.ds_interval(ds_5min)

            uxds_model_closest=latlon_verify.return_closest_cells(uxds_model,ds_obs)
            uxds_ter_closest=latlon_verify.return_closest_cells(uxds_ter,ds_obs)

            # 6 HOUR PRECIP STATS
            for id in ds_obs['station_id'].values:

                ds_stats['obs_elevation'].loc[
                    dict(
                        station_id=id
                    )
                ] = ds_obs['elevation'].sel(station_id=id)[0].values/12

                ds_stats['model_elevation'].loc[
                    dict(
                        station_id=id
                    )
                ] = uxds_ter_closest['ter'].sel(station_id=id).values


                ds_stats['lat'].loc[
                                    dict(
                                        station_id=id
                                    )
                                ] = ds_obs['lat'].sel(station_id=id)[0].values/12

                ds_stats['lon'].loc[
                                    dict(
                                        station_id=id
                                    )
                                ] = ds_obs['lng'].sel(station_id=id)[0].values/12


                ds_stats[f'obs_pcp_18_to_24'].loc[
                    dict(
                        init_ts=date,
                        model=model,
                        station_id=id
                    )
                ] = ds_obs['value'].sel(station_id=id)[18:24].sum().values


                ds_stats[f'model_pcp_18_to_24'].loc[
                    dict(
                        init_ts=date,
                        model=model,
                        station_id=id
                    )
                ] = uxds_model_closest['pcp6'].sel(station_id=id).values[3]    

                ds_stats[f'obs_pcp_00_to_24'].loc[
                    dict(
                        init_ts=date,
                        model=model,
                        station_id=id
                    )
                ] = ds_obs['value'].sel(station_id=id)[:24].sum().values


                ds_stats[f'model_pcp_00_to_24'].loc[
                    dict(
                        init_ts=date,
                        model=model,
                        station_id=id
                    )
                ] = uxds_model_closest['pcp24'].sel(station_id=id).values[0]


    ds_stats['obs_pcp_00_to_24']=ds_stats['obs_pcp_00_to_24'][:,:,0].drop_vars('model')
    ds_stats['obs_pcp_18_to_24']=ds_stats['obs_pcp_18_to_24'][:,:,0].drop_vars('model')

    ds_stats['bias_pcp_18_to_24']=ds_stats['model_pcp_18_to_24']-ds_stats['obs_pcp_18_to_24']
    ds_stats['bias_pcp_00_to_24']=ds_stats['model_pcp_00_to_24']-ds_stats['obs_pcp_00_to_24']

    ds_stats['bias_elevation']=ds_stats['model_elevation']-ds_stats['obs_elevation']

    return ds_stats