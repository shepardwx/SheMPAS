import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import json
import os
import argparse


# THIS SCRIPT DOWNLOADS OBS FROM THE HAWAII MESONET WEBSITE, USING MY API KEY (hidden)
with open('/home/disk/hot/shepard0/keys/hawaii_mesonet.txt') as f:
    lines = f.readlines()
key=lines[0]



# location='all'


ids_maui=[115,
          116,
          118,
          119,
          121,
          122,
          131,
          132,
          133,
          134,
          141,
          143,
          144,
          145,
          151,
          152,
          153,
          154,
          155,
          156,
          157,
          160,
          161,
          165,
          166]


ids_hawaii=[201,
     202,
     204,
     211,
     212,
     213,
     214,
     221,
     222,
     231,
     235,
     241,
     242,
     243,
     244,
     246,
     251,
     252,
     253,
     254,
     255,
     256,
     257,
     258,
     281,
     282,
     283,
     286,
     287,
     288]



ids_unknown=[301,
             302]

ids_molokai=[411,
             412,
             421,
             422,
             431]

ids_oahu=[501,
          502,
          503,
          504,
          505,
          506,
          521,
          531,
          533,
          541,
          551,
          552]

ids_kauai=[601,
           602,
           603,
           611,
           621,
           641]




def which_ids(location):
    if location=='all':
        ids=ids_maui+ids_hawaii+ids_unknown+ids_molokai+ids_oahu+ids_kauai
    elif location=='maui':
        ids=ids_maui
    elif location=='hawaii':
        ids=ids_hawaii
    elif location=='unknown':
        ids=ids_unknown
    elif location=='molokai':
        ids=ids_molokai
    elif location=='oahu':
        ids=ids_oahu
    elif location=='kauai':
        ids=ids_kauai

    return ids



def ds_5min(start_date,end_date,location='all'):
    var_id='RF_1_Tot300s'
    fname=f'response_{start_date}_{end_date}.json'

    station_ids=which_ids(location)

    stns_str=f'0{station_ids[0]}'
    for id in station_ids[1:]:
        stns_str=stns_str+f'%2C0{id}'


    curl_fstr=f'curl -X \'GET\' \'https://api.hcdp.ikewai.org/mesonet/db/measurements?station_ids={stns_str}&start_date={start_date[:4]}-{start_date[4:6]}-{start_date[6:8]}T{start_date[8:10]}%3A00%3A00&end_date={end_date[:4]}-{end_date[4:6]}-{end_date[6:8]}T{end_date[8:10]}%3A00%3A00&var_ids={var_id}&location=hawaii&limit=1000000&reverse=reverse&join_metadata=true&row_mode=array\' -H \'accept: */*\' -H \'Authorization: Bearer {key}\' -o {fname}'

    os.system(curl_fstr)

    # Load the JSON file
    with open(fname, "r", encoding="utf-8") as f:
        raw = json.load(f)

    # Column names matching the values in each data row
    columns = [
        "timestamp",
        "station_id",
        "variable",
        "value",
        "flag",
        "units",
        "units_plain",
        "units_expanded",
        "variable_display_name",
        "interval_seconds",
        "station_name",
        "lat",
        "lng",
        "elevation",
        "extra_1",
        "extra_2",
        "extra_3",
        "extra_4",
    ]

    # Convert rows to a DataFrame
    df = pd.DataFrame(raw["data"], columns=columns)

    # Convert columns to appropriate types
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["interval_seconds"] = pd.to_numeric(df["interval_seconds"], errors="coerce")
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lng"] = pd.to_numeric(df["lng"], errors="coerce")
    df["elevation"] = pd.to_numeric(df["elevation"], errors="coerce")

    # Sort chronologically
    df = df.sort_values(["station_id", "timestamp"])

    # Convert to xarray
    ds = df.set_index(["station_id", "timestamp"]).to_xarray()

    os.system(f'rm -rf {fname}')

    return ds

def ds_interval(ds,interval=1):


    my_ds=ds.copy()
    my_ds['value']=my_ds['value']
    my_ds['value']=np.round(my_ds['value'],2)

    # shifted = my_ds.assign_coords(timestamp=my_ds.timestamp - pd.Timedelta(minutes=5))

    hourly = my_ds.resample(
        timestamp=f"{interval}h",
        offset="5min",
        closed="left",
        label="right"
    ).sum(
        dim="timestamp",
        min_count=1
    )

    index_array = xr.DataArray(np.arange(len(hourly['timestamp']))[1:], dims="timestamp")
    indexed_data = hourly.isel(timestamp=index_array)
    final_ds = indexed_data.assign_coords(
        timestamp=indexed_data.timestamp - pd.Timedelta(minutes=5)
    )

    return final_ds


def interval_18_to_24(ds):
    ds_1hr=ds_interval(ds)

    result = ds_1hr.isel(timestamp=slice(18, 24)).sum(dim="timestamp")

    return result


def interval_0_to_24(ds):
    ds_1hr=ds_interval(ds)

    result = ds_1hr.isel(timestamp=slice(0, 24)).sum(dim="timestamp")

    return result



def print_sums(ds,location='all'):

    station_ids=which_ids(location)

    for id in station_ids:
        if f'0{id}' in ds.station_id:
            print(ds['station_name'].sel(station_id=f'0{id}')[0].values)
            print(f'{ds['lat'].sel(station_id=f'0{id}')[0].values}, {ds['lng'].sel(station_id=f'0{id}')[0].values}')
            print(f'Rain Sum: {str(np.nansum(ds['value'].sel(station_id=f'0{id}')[1:])/25.4)} in'+'\n')
        else:
            print(f'{id} NOT FOUND')


if __name__ == "__main__":

    # COMMAND LINE ARGUMENTS
    parser=argparse.ArgumentParser()
    parser.add_argument('-s','--start_date',type=str)
    parser.add_argument('-e','--end_date',type=str)
    parser.add_argument('-i','--interval',type=int,default=24)
    parser.add_argument('-l', '--location',type=str,default='all')

    args = parser.parse_args()

    start_date = args.start_date
    end_date=args.end_date
    interval=args.interval
    location=args.location


    ds=ds_5min(start_date=start_date,end_date=end_date,location=location)
    print_sums(ds,location=location)
    # print(ds['value'].sel(station_id=f'0{552}'))