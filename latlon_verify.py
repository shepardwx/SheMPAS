import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import glob
import uxarray as ux
from scipy.spatial import cKDTree
import hawaii_mesonet_downloader



# REPEATEDLY USEFUL DSs

grid_path='/home/disk/hot/shepard0/mpas/meshes/regional/4/hawaii_shifted/hawaii_shifted.grid.nc'
init_path='/home/disk/hot/shepard0/mpas/regional/hawaii_shifted/4/2026031412/hawaii_shifted.init.nc'

ds = xr.open_dataset(grid_path)
ds_init = xr.open_dataset(init_path)
uxds_init = ux.open_dataset(grid_path,init_path)

#   FUNCTION THAT WILL TAKE A DATE (yyyymmddhh) AND MODEL (RCON_gf, Thompson_gf, etc.) and return a ds from that model's latlon.nc, with time fixed

def fix_time_return_ds(date,model):
    latlon_dir=f'/home/disk/hot/shepard0/mpas/regional/hawaii_shifted/4/{date}/{model}'
    ds=xr.open_dataset(f'{latlon_dir}/latlon.nc')

    hist_pat=latlon_dir+'/hist*'
    hist_files = [(file, index) for index, file in enumerate(glob.glob(hist_pat))]

    xtimes=[]
    for tuple in hist_files:
        temp_ds=xr.open_dataset(tuple[0])
        xtimes.append(temp_ds['xtime'].values[0][:19])

    times=[]
    for xtime in xtimes:
        # time=pd.to_datetime(xtime,format='mixed')
        time=pd.to_datetime(xtime.decode(),format='%Y-%m-%d_%H:%M:%S')
        times.append(time)

    ds=ds.assign_coords(time=('Time',times))
    # os.system(f'ls {latlon_dir}')
    # print(hist_files)

    ds['rain_total']=ds['rainnc']
    return ds





def return_closest_cells(uxds,ds_obs):

    R = 6371229.0 

    xyz_cell = np.column_stack([ds.xCell, ds.yCell, ds.zCell])
    tree = cKDTree(xyz_cell)

    lat, lon = np.radians(ds_obs.lat.isel(timestamp=0).values/12),np.radians(ds_obs.lng.isel(timestamp=0).values/12)
    xyz_stn = np.column_stack([np.cos(lat)*np.cos(lon),
                            np.cos(lat)*np.sin(lon),
                            np.sin(lat)])

    dist, idx = tree.query(xyz_stn)
    dist_m = dist * R 

    nbrs = tree.query_ball_point(xyz_stn, r=15e3 / R)


    uxds_cells=uxds.isel(n_face=idx)
    uxds_cells=uxds_cells.assign_coords(station_id=('n_face',ds_obs['station_id'].values))

    return uxds_cells 



def bias_18_to_24_closest(uxds_model_closest,ds_obs):

    uxda_bias=uxds_model_closest['pcp6'].isel(time6=3)-ds_obs['value'][18:24].sum(dim='timestamp').values
    uxda_bias=uxda_bias.assign_attrs(long_name='6-hr Precipitation Bias, valid at FHR 24',units='mm')

    # bias=uxds_model_closest['pcp6'].isel(time6=3).sel(station_id=id)-ds_obs['value'][18:24].sum(dim='timestamp').sel(station_id=id).values
    # uxda_bias=uxda_bias.assign_attrs(long_name='6-hr Precipitation Bias, valid at FHR 24',units='mm')



    return uxda_bias











#   FUNCTION THAT TAKES A DATE (yyyymmddhh), MODEL (i.e. RCON_gf), id ('0115'), ds_obs and ds_mpas, and returns the index locations of the 
#   closest, and second closest lat and lons in the MPAS data to that id's lat/lon location

# def find_closest_lat_lon_i(id,ds_obs,ds_mpas):

#     # obs ds
#     ds_obs=ds_obs.copy()

#     #mpas lat/lon
#     mpas_lats=ds_mpas.latitude.values
#     mpas_lons=ds_mpas.longitude.values

#     #station lat/lon
#     station_lat=float(ds_obs['lat'].isel(timestamp=0).sel(station_id=id))
#     station_lon=float(ds_obs['lng'].isel(timestamp=0).sel(station_id=id))


#     # #differences
#     lat_difs=mpas_lats-station_lat
#     lon_difs=mpas_lons-station_lon


#     # #ABSOLUTE minimum args
#     lat_arg_min=np.argmin(np.abs(lat_difs))
#     lon_arg_min=np.argmin(np.abs(lon_difs))

#     # #NEXT smallest args
#     if mpas_lats[lat_arg_min] < station_lat:
#         lat_arg_min2=lat_arg_min+1
#     elif mpas_lats[lat_arg_min] > station_lat:
#         lat_arg_min2=lat_arg_min-1

#     if mpas_lons[lon_arg_min] < station_lon:
#         lon_arg_min2=lon_arg_min+1
#     elif mpas_lons[lon_arg_min] > station_lon:
#         lon_arg_min2=lon_arg_min-1

#     # print(lat_arg_min,lon_arg_min)
#     # print(f'STN                 : {station_lat}, {station_lon}')
#     # print(f'MPAS                : {mpas_lats[lat_arg_min]}, {mpas_lons[lon_arg_min]} [{mpas_lats[lat_arg_min]-station_lat}, {mpas_lons[lon_arg_min]-station_lon}]')
#     # print(f'MPAS, next closest  : {mpas_lats[lat_arg_min2]}, {mpas_lons[lon_arg_min2]} [{mpas_lats[lat_arg_min2]-station_lat}, {mpas_lons[lon_arg_min2]-station_lon}]')

#     return lat_arg_min,lon_arg_min,lat_arg_min2,lon_arg_min2


# #   FUNCTION that returns an array of MPAS values (dim time) at the location of the station, by bilinear interpolation

# def my_bilinear_interp(id,ds_obs,ds_mpas,var='rain_total'):

#     # ds_obs
#     ds_obs = ds_obs.copy()

#     #mpas lat/lon
#     mpas_lats=ds_mpas.latitude.values
#     mpas_lons=ds_mpas.longitude.values

#     #station lat/lon
#     station_lat=float(ds_obs['lat'].isel(timestamp=0).sel(station_id=id))
#     station_lon=float(ds_obs['lng'].isel(timestamp=0).sel(station_id=id))

#     lat_arg_min,lon_arg_min,lat_arg_min2,lon_arg_min2=find_closest_lat_lon_i(id,ds_obs,ds_mpas)

#     a_coords = [mpas_lats[lat_arg_min],mpas_lons[lon_arg_min]]
#     b_coords = [mpas_lats[lat_arg_min],mpas_lons[lon_arg_min2]]
#     c_coords = [mpas_lats[lat_arg_min2],mpas_lons[lon_arg_min2]]
#     d_coords = [mpas_lats[lat_arg_min2],mpas_lons[lon_arg_min]]

#     # a_stn_lat_dif=np.abs(station_lat-a_coords[0])
#     # b_stn_lat_dif=np.abs(station_lat-b_coords[0])
#     # c_stn_lat_dif=np.abs(station_lat-c_coords[0])
#     # d_stn_lat_dif=np.abs(station_lat-d_coords[0])

#     # a_stn_lon_dif=np.abs(station_lon-a_coords[1])
#     # b_stn_lon_dif=np.abs(station_lon-b_coords[1])
#     # c_stn_lon_dif=np.abs(station_lon-c_coords[1])
#     # d_stn_lon_dif=np.abs(station_lon-d_coords[1])


#     # find middle values at each lat
#     lat1_val=(ds_mpas[var].isel(latitude=lat_arg_min,longitude=lon_arg_min) *
#                     np.abs((station_lon-mpas_lons[lon_arg_min2])/(mpas_lons[lon_arg_min]-mpas_lons[lon_arg_min2])) ) + (ds_mpas[var].isel(latitude=lat_arg_min,longitude=lon_arg_min2) *
#                     np.abs((station_lon-mpas_lons[lon_arg_min])/(mpas_lons[lon_arg_min]-mpas_lons[lon_arg_min2])) )

#     lat2_val=(ds_mpas[var].isel(latitude=lat_arg_min2,longitude=lon_arg_min) *
#                     np.abs((station_lon-mpas_lons[lon_arg_min2])/(mpas_lons[lon_arg_min]-mpas_lons[lon_arg_min2])) ) + (ds_mpas[var].isel(latitude=lat_arg_min2,longitude=lon_arg_min2) *
#                     np.abs((station_lon-mpas_lons[lon_arg_min])/(mpas_lons[lon_arg_min]-mpas_lons[lon_arg_min2])) )

#     # THEN find the middle value at the right lon

#     lon_val=(lat1_val*
#                 np.abs((station_lat-mpas_lats[lat_arg_min2])/(mpas_lats[lat_arg_min]-mpas_lats[lat_arg_min2])))+(lat2_val*
#                 np.abs((station_lat-mpas_lats[lat_arg_min])/(mpas_lats[lat_arg_min]-mpas_lats[lat_arg_min2])))

#     final_val=lon_val

#     return(final_val)


def bias_18_to_24(uxds_model,ds_obs):

    # start_date=ds_model['time'].values[0].astype("datetime64[ms]").tolist().strftime("%Y-%m-%d")

    ds_both=ds_obs.copy()

    obs_sums=hawaii_mesonet_downloader.interval_18_to_24(ds_obs)

    interped_data=[]
    obs_list=[]
    biases=[]
    # for id in ds_obs['station_id']:

    #     this_interp=my_bilinear_interp(id=id,ds_obs=ds_obs,ds_mpas=ds_model)
    #     interp_18_to_24=this_interp.sel(time=this_interp['time'].values[0]+np.timedelta64(24, "h"))-this_interp.sel(time=this_interp['time'].values[0]+np.timedelta64(18, "h"))

    #     obs_sum=np.round(obs_sums.sel(station_id=id)['value'].values,2)

    #     obs_list.append(obs_sum)
    #     interped_data.append(interp_18_to_24.values)
    #     biases.append(interp_18_to_24.values-obs_sum)

        

    da_model = xr.DataArray(
    interped_data,
    dims=["station_id"],
    coords={"station_id": ds_both["station_id"]},
    name="model_pcp_18_to_24",
    )

    da_obs = xr.DataArray(
    obs_list,
    dims=["station_id"],
    coords={"station_id": ds_both["station_id"]},
    name="obs_pcp_18_to_24",
    )

    da_bias = xr.DataArray(
    biases,
    dims=["station_id"],
    coords={"station_id": ds_both["station_id"]},
    name="bias_pcp_18_to_24",
    )

    ds_both['model_pcp_18_to_24']=da_model
    ds_both['obs_pcp_18_to_24']=da_obs
    ds_both['bias_pcp_18_to_24']=da_bias

    return(ds_both)

# def bias_00_to_24(ds_model,ds_obs):

#     # start_date=ds_model['time'].values[0].astype("datetime64[ms]").tolist().strftime("%Y-%m-%d")

#     ds_both=ds_obs.copy()

#     obs_sums=hawaii_mesonet_downloader.interval_0_to_24(ds_obs)

#     interped_data=[]
#     obs_list=[]
#     biases=[]
#     for id in ds_obs['station_id']:

#         this_interp=my_bilinear_interp(id=id,ds_obs=ds_obs,ds_mpas=ds_model)
#         interp_00_to_24=this_interp.sel(time=this_interp['time'].values[0]+np.timedelta64(24, "h"))-this_interp.sel(time=this_interp['time'].values[0]+np.timedelta64(0, "h"))

#         obs_sum=np.round(obs_sums.sel(station_id=id)['value'].values,2)

#         obs_list.append(obs_sum)
#         interped_data.append(interp_00_to_24.values)
#         biases.append(interp_00_to_24.values-obs_sum)

        

#     da_model = xr.DataArray(
#     interped_data,
#     dims=["station_id"],
#     coords={"station_id": ds_both["station_id"]},
#     name="model_pcp_00_to_24",
#     )

#     da_obs = xr.DataArray(
#     obs_list,
#     dims=["station_id"],
#     coords={"station_id": ds_both["station_id"]},
#     name="obs_pcp_00_to_24",
#     )

#     da_bias = xr.DataArray(
#     biases,
#     dims=["station_id"],
#     coords={"station_id": ds_both["station_id"]},
#     name="bias_pcp_00_to_24",
#     )

#     ds_both['model_pcp_00_to_24']=da_model
#     ds_both['obs_pcp_00_to_24']=da_obs
#     ds_both['bias_pcp_00_to_24']=da_bias

#     return(ds_both)
    


# if __name__ == "__main__":

#     date='2026031512'
#     model1='RCON_gf'
#     model2='Thompson_gf'
#     ds_rcon=fix_time_return_ds(date=date,model=model1)
#     ds_thompson=fix_time_return_ds(date=date,model=model2)
#     ds_obs=hawaii_mesonet_downloader.ds_5min(start_date=date,end_date='2026031712')


#     rcon_verification_18_to_24=bias_18_to_24(ds_model=ds_rcon,ds_obs=ds_obs)
#     rcon_verification_00_to_24=bias_00_to_24(ds_model=ds_rcon,ds_obs=ds_obs)

#     thompson_verification_18_to_24=bias_18_to_24(ds_model=ds_thompson,ds_obs=ds_obs)
#     thompson_verification_00_to_24=bias_00_to_24(ds_model=ds_thompson,ds_obs=ds_obs)

#     plt.plot(rcon_verification_18_to_24['station_id'],rcon_verification_18_to_24['bias_pcp_18_to_24'],label='RCON')
#     plt.plot(rcon_verification_18_to_24['station_id'],thompson_verification_18_to_24['bias_pcp_18_to_24'],label='Thompson')
#     plt.title(f'Init {date}. 6-hr precipitation bias FHR 24')
#     plt.legend()
#     plt.show()


#     plt.plot(rcon_verification_18_to_24['station_id'],rcon_verification_00_to_24['bias_pcp_00_to_24'],label='RCON')
#     plt.plot(rcon_verification_18_to_24['station_id'],thompson_verification_00_to_24['bias_pcp_00_to_24'],label='Thompson')
#     plt.title(f'Init {date}. 24-hr precipitation bias FHR 24')
#     plt.legend()
#     plt.show()

#     # for id in ds_obs['station_id'][:]:
#     #     print('\nID:',str(id.values))
#     #     # print(find_closest_lat_lon_i(date=date,model=model,id=id,ds_obs=ds_obs,ds_mpas=ds_mpas))
#     #     interped_data_rcon=my_bilinear_interp(date=date,model=model1,id=id,ds_obs=ds_obs,ds_mpas=ds_rcon)
#     #     interped_data_thompson=my_bilinear_interp(date=date,model=model2,id=id,ds_obs=ds_obs,ds_mpas=ds_thompson)
#     #     # obs_data=hawaii_mesonet_downloader.ds_interval(ds_obs.sel(station_id=id),interval=1)['value'].values*25.4
#     #     obs_data=hawaii_mesonet_downloader.interval_18_to_24(ds=ds_obs).sel(station_id=id)


#     #     print('RCON:',interped_data_rcon[24]-interped_data_rcon[18])
#     #     print('Thompson:',interped_data_thompson[24]-interped_data_thompson[18])
#     #     print('Obs:',obs_data['value'].values)