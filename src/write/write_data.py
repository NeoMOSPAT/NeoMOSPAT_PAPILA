# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr
from os import makedirs
from os.path import join
import pandas as pd
import netCDF4
from datetime import datetime
import src.IncludeFile as IncF
import src.common as aux
from src.write import write_stations 


IncF.i_Write = [int(c) for c in str(IncF.i_Write)]


##############################################
def sjoin(x): 
    return  x[0]   

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


##############################################
def main(t_ObsStationData, t_ModelStationData, t_ModelData, t_Stations_Info, t_ModelStation_Info, t_Models_Info, t_Stations_Filters, t_Model_Filters):
 
    if 0 in IncF.i_Write:
        print("...Nothing was asked to be saved")
        return None 

    if 2 in IncF.i_Write:
        print("...Saving observation from stations")
        write_stations.main(t_ObsStationData, t_Stations_Info, t_Stations_Filters)

    if 3 in IncF.i_Write:
        print("...Saving model at observation points")
        write_stations.main(t_ModelStationData, t_ModelStation_Info, t_Stations_Filters, model=True)

    return None

    timeframes=IncF.c_TimeMeanRes
    networks=list(set(["_".join(net.split("_")[:-1]) for net in list(t_ObsStationData[timeframes[0]].keys())]))

    categories=aux.read_yaml("vars")
    for network in networks:
        for timeframe in timeframes:
            data=t_ObsStationData[timeframe]
            keys=[key for key in data.keys() if network=="_".join(key.split("_")[:-1])]
            units=data[keys[0]].meta.units 

            merged_data=pd.concat([data[k] for k in keys], axis=1).groupby(level=[0,1], axis=1).apply(lambda x: x.apply(sjoin, axis=1))
            

            #########
            stations=set(merged_data.columns.get_level_values(1))
            save_dir=join(IncF.c_FigureDir, "Data", network, timeframe)
            makedirs(save_dir, exist_ok=True)
            t_Stations_Info.loc[keys].drop_duplicates(keep="first").to_csv(join(save_dir, "Chile-Stations.csv"))

            for ID in stations:
                sta_data=merged_data.xs(ID, axis=1, level=1).dropna(how="all")
                
                cal_columns=[c+"_"+units[c]  for c in sta_data.columns if c in categories["Cal"].split(" ")]
                met_columns=[c+"_"+units[c]  for c in sta_data.columns if c in categories["Met"].split(" ")]
                sta_data.columns=[c+"_"+units[c] for c in sta_data.columns]
                cal_data=sta_data[cal_columns].dropna(how="all")
                met_data=sta_data[met_columns].dropna(how="all")
                if not cal_data.empty:
                    cal_data.to_csv(join(save_dir, f"ID-{ID}--Cal__{timeframe}.csv"))
                    print("   ...Writing data to ",join(save_dir, f"ID-{ID}--Cal__{timeframe}.csv"))
                if not met_data.empty:
                    met_data.to_csv(join(save_dir, f"ID-{ID}--Met__{timeframe}.csv"))
                    print("   ...Writing data to ",join(save_dir, f"ID-{ID}--Met__{timeframe}.csv"))
                

    """
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
    """
