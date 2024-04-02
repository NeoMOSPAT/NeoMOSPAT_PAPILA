# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr
from os import makedirs
import pandas as pd
import netCDF4
from datetime import datetime




##############################################
## Saves netcdf files
def write_ncfile(new_model, variables, mod_name, date, mesh, direc="", freq_groupby="", freq_data=""):
    ############
    new_model.sort_values(by=["i_lat", "i_lon"], inplace=True)
    mesh.sort_values(by=["i_lat", "i_lon"], inplace=True)
    units=new_model.meta.units

    ## DateTime index
    time=new_model.index.drop_duplicates(keep="first")
    ## DateTime objectsList
    time=list(map(lambda x: datetime.strptime(x, "%Y-%m-%d %H:%M:%S"), time.strftime("%Y-%m-%d %H:%M:%S").tolist())) 

    ## In numbers for netcdf
    time_unit_out="hours since %s"%(time[0].strftime("%Y-%m-%d %H:%M:%S"))
    time=netCDF4.date2num(time, time_unit_out)

    ##
    lats=mesh.drop_duplicates(["i_lat"], keep="first")["i_lat"]
    lons=mesh.drop_duplicates(["i_lon"], keep="first")["i_lon"]


    ## Dimensions and coordinates
    coords = dict(time=time, x=lats, y=lons)
    coords = {'x': ('x', lats), 'y': ('y', lons), 'time': ('time', time, {'units':time_unit_out, 'long_name':'Time', "calendar":"standard" }) }
    dims = ("time", 'x', 'y')
  

    

    ## Temporal store for DataArrays
    temporal_dataset = {}
    ds = xr.DataArray(np.reshape(mesh["lat"].values, (-1, len(lons))), dims=('x', 'y'), coords=dict(x=lats, y=lons))
    ds.attrs['units'] = 'degrees_north'
    ds.attrs['long_name'] = "Latitude"
    temporal_dataset["lat"]=ds

    ds = xr.DataArray(np.reshape(mesh["lon"].values, (-1, len(lons))), dims=('x', 'y'), coords=dict(x=lats, y=lons))
    ds.attrs['units'] = 'degrees_east'
    ds.attrs['long_name'] = "Longitude" 
    temporal_dataset["lon"] = ds


    ## Makes it so that the reshaping is correct, orders by dimension
    new_model=new_model.reset_index().sort_values(by=["time", "i_lat", "i_lon"])
    for var in variables:
        ## Reshapes data to have sahpe of lats and lons
        bcvar=new_model[var].replace(0, -999.9).values
        bcvar=np.reshape(bcvar, (len(time), -1, len(lons)))
        ds = xr.DataArray(bcvar, dims=dims, coords=coords)
        ds.attrs['_FillValue'] = -999.9
        ds.attrs['units'] = units[var]
        temporal_dataset[var] = ds

    ## Create DataSet
    DT=xr.Dataset(temporal_dataset, attrs = {"description":"Bias Corrected"})

    ## Save to file
    filename="%s/%s/%s%s_%s_%s.nc"%(direc, mod_name, freq_groupby, freq_groupby, mod_name, date)
    DT.to_netcdf(filename, format="NETCDF4_CLASSIC")
    print("   ...Saved to:  %s"%(filename), flush=True)


    return None


def main(models_data, models_info, dics_to_save, freq, direc=""):
    ############
    freq_data, freq_groupby=freq.split(",")
    for mod, variables in dics_to_save.items():
        mod_name=mod.replace(" ","-")
        mesh=models_info[mod]
        print(f"...Writing netCDF files of model  {mod}", flush=True)
        makedirs("%s/%s"%(direc, mod_name), exist_ok = True)

        units=models_data[freq_data][mod].meta.units 
        data=models_data[freq_data][mod].reset_index().set_index("time")
        groupby=data.groupby(pd.Grouper(freq=freq_groupby))
        for  group, df_group in groupby:
            df_group.meta.units=units
            date=group.strftime("%Y%m%d00")
            if df_group.empty:
                continue
            write_ncfile(df_group, variables, mod_name, date,mesh, direc=direc, freq_groupby=freq_groupby, freq_data=freq_data)
