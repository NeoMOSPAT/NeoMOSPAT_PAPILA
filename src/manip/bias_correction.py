# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.ops import nearest_points, transform
from shapely.geometry import Point, Polygon
import netCDF4
#
from functools import partial
import pyproj
import geopy.distance
import src.common as aux
from src.write  import write_netcdf 
from src import common as aux
from copy import deepcopy
from src.plot import Surface_Maps as maps
#
import os
os.environ["PYTHONBREAKPOINT"]="pudb.set_trace"
from datetime import datetime, timedelta

mean_type="end"



######################################################
## Creates the bias factor between grid point and station for a daily period and a window of 10 days of rolling mean.
def bias_station(station_data, model_data, var, mean_type="end", window=10):
    ## Place in one dataframe
    both=pd.merge(model_data, station_data, left_index=True, right_index=True)
    #.dropna()
    both=both[[col for col in both.columns if var in col]]

    ## Mean at end of periodfactors_stations.dropna(axis=1, how="all")
    both=both.fillna(method='ffill')
    ## Delete
    roll_mean=both.rolling(window=window, min_periods=1).mean()#min_periods=1
    
    if mean_type=="center":
        roll_mean.index=roll_mean.index.shift(-int(window/2))
    elif mean_type=="start":
        roll_mean.index=roll_mean.index.shift(-window)

    model_col=[col for col in roll_mean.columns if len(col.split())==2][0]
    obs_col=[col for col in roll_mean.columns if col==var][0]

    roll_mean["factor"]=roll_mean[obs_col]/roll_mean[model_col]
    return roll_mean[["factor"]]


######################################################
## Obtain bias correction for all stations in mini with closest grid point 
def bias_all_stations(mini,  model, model_name, data, var, mean_type="end", window=10):

    ####################################################
    factors=[]
    ### Get data nearest to stations
    for ID in mini["ID"]:
        row=mini.loc[mini["ID"]==ID]
        ## Find indexes to extract model data
        nearest=row[["i_lon", "i_lat"]].values[0]
        ## Extract Model data
        model_data=model.loc[:, nearest[1], nearest[0]]
        model_data.columns=['%s  %s'%(model_name, x) for x in model_data.columns]

        ## Extract Station data
        station_data=aux.get_station_data(data, ID)

        ## Obtain bias factor correction
        factor_station=bias_station(station_data, model_data, var, mean_type=mean_type, window=window)
        factor_station.columns=[ID]

        ## Add factor to dictionary
        factors.append(factor_station)

    factors=pd.concat(factors, axis=1)

    return factors.dropna(axis=1, how="all")


######################################################
## Find grid points that are within radius of all stations
proj_4326 = pyproj.Proj('+proj=longlat +datum=WGS84')
def grids_within_radius(mini, mesh, radius):
    ######################################################
    ## Return Polygon of circle centered at lat, lon with radius km
    def get_circle(lat, lon, km):
        ## Gets points in circle around lat, lon with radius of km
        def geodesic_point_buffer(lat, lon, km):
            aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
            project = partial(
                pyproj.transform,
                pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)), proj_4326)
            #'<authority>:<code>'
            buf = Point(0, 0).buffer(km * 1000)  # distance in metres
            return transform(project, buf).exterior.coords[:]

        points_list=geodesic_point_buffer(lat, lon, km)
        return Polygon(points_list)

    ####################################################
    ##########
    ## Add circular polygon rodenado station as geometry
    geometry=mini.apply(lambda a: get_circle(a.Latitud, a.Longitud, radius), axis=1)
    mini=gpd.GeoDataFrame(mini, crs='epsg:4326', geometry=geometry)
    mini["Station_coords"]=[Point(xy) for xy in zip(mini.Longitud, mini.Latitud)]

    ##########
    ## Find intersection points
    def get_points_inside(data):
        geometry=data.geometry
        df= mesh.loc[mesh["geometry"].within(geometry)]

        df["ID"]=data.ID
        df["Station_coords"]=data.Station_coords
        return df

    intersection=[get_points_inside(row) for index, row in mini.iterrows()]
    intersection=pd.concat(intersection).sort_values(by=["i_lat", "i_lon"]).reset_index(drop=True)

    if intersection.empty:
        intersection["distance_geodesic"]=None
        return intersection

    ##########
    ## Get distance between point and station
    intersection["distance_geodesic"]=intersection.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.Station_coords.y, a.Station_coords.x)), axis=1 )
    intersection=intersection.sort_values(by="distance_geodesic")

    return intersection


######################################################
## Generate grid with all their respective biases. Considers distance, and multiple biases factors from multiple stations. 
def get_factor_grid(intersection, factors, model, var, radius=10):
    factor_info=intersection[["i_lat", "i_lon", "lat", "lon", "distance_geodesic", "ID"]]
    factor_info=factor_info.reset_index(drop=True)


    ##########
    ## Find the perdentage with which it will go down. From 1 at the top for full strength to 1/factor when there is no strength
    ## Retuturn multiplier given the distance
    def line_eq(distance, rad=radius*1000):
        distance=distance.meters
        ## Equation_ multiplier=(1/f-1)/10*distance +1
        ## Equation total_factor=f(1-distance/10)+distance/10
        if distance>=rad:
            ## First number accompanies factor, second is number to sum
            return [0, 1]

        else:
            ## First number accompanies factor, second is number to sum
            return [1-distance/rad, distance/rad]


    ##########
    ## Find function: multiplier*f +sum of how the factor f changes with distance
    column=factor_info.apply(lambda a: line_eq(a.distance_geodesic), axis=1)
    factor_info[["multiplier", "sum"]]=pd.DataFrame(data=column.values.tolist(), columns=["multiplier", "sum"])
    factor_info=factor_info.sort_values(by="distance_geodesic")[["i_lat", "i_lon", "lat", "lon", "multiplier",  "sum", "ID", "distance_geodesic"]]



    ####################################################  
    ## Model at points within radious      
    if factor_info.empty:
        #corrected=model.reset_index(level=0).loc[intersection.indexes].reset_index().set_index(["time", "i_lat", "i_lon"])[[var]] 
        #corrected[var]=1
        corrected=pd.DataFrame(columns=[var, "time", "i_lat", "i_lon"]).set_index(["time", "i_lat", "i_lon"]) 
        corrected[var]=pd.to_numeric(corrected[var], errors="coerce")
        return corrected.sort_index()


    ####################################################  
    ## Model at points within radious  
    dfs=[]
    for name, group in factor_info.groupby("ID"):
        factor_id=factors[[name]]
        factor_id.index.name="time"

        group=group.reset_index(drop=True)

        ## change this lat lon for i_lat i_lon
        group=group.set_index(['i_lat', 'i_lon'])
        corrected=pd.DataFrame((factor_id.values * group["multiplier"].values + group["sum"].values)*1, index=factor_id.index, columns=group.index).unstack().to_frame(name=var)

        dfs.append(corrected)

    corrected=pd.concat(dfs, axis=0).sort_index()
    corrected=corrected.reorder_levels(["time", "i_lat", "i_lon"])
    corrected[var]=pd.to_numeric(corrected[var], errors="coerce")
    corrected=corrected.groupby(level=[0,1,2]).mean()

    ##
    return corrected.fillna(1)



######################################################
## Find new corrected var values
def corrected_model(mini, mini_mod, data, model, mesh, var, model_name, mean_type="end", window=1):
 
    ## Obtain bias correction factor of all stations in mini
    print("   ...Finding bias factor of stations", flush=True)
    
    factors_stations=bias_all_stations(mini_mod, model, model_name, data, var, mean_type=mean_type, window=window)
    mini_mod=mini_mod[mini_mod["ID"].isin(factors_stations.columns)]

    ## Obtain points per station that are within radius of window
    corrected=model[[var]]
    for radius in IncF.i_BC_Rads:
        print("   ...Finding grid points inside radius of ", radius, flush=True)

        intersection=grids_within_radius(mini_mod, mesh, radius)


        #############
        ## Obtain grid factor for all grid points inside at least one circle of a station. If point has multiple stations in range, then the average of the factor of each one (counting distance) is computed.
        print("      ...Obtaining bias correction factor for grid points", flush=True)
        grid_factors=get_factor_grid(intersection, factors_stations, model, var, radius=radius)
        final_factor=grid_factors.copy()

        """
        bounds=mesh.total_bounds
        colorbar=[np.nanpercentile(final_factor, 2), np.nanpercentile(final_factor, 98)]

        data_for_plot=pd.merge(final_factor.reset_index(), mesh, on=["i_lat", "i_lon"], how="inner").set_index(["time", "i_lat", "i_lon"]).sort_index().loc["2016-07-01 10:00:00":"2016-07-01 14:00:00", :, :].reset_index(level=[1,2])
        maps.plot_maps(data_for_plot, ["PM25"], {"PM25":"ug/m3"}, "R%s"%radius, bounds, mesh, data, mini, "H", colorbar, direc="Trial_Factors/", filter="")
        data_for_plot=pd.merge(final_factor.reset_index(), mesh, on=["i_lat", "i_lon"], how="inner").set_index(["time", "i_lat", "i_lon"]).sort_index().loc["2016-07-01 10:00:00":"2016-07-01 14:00:00", :, :].reset_index(level=[1,2]).replace({"PM25":{0:np.nan}})
        maps.plot_maps(data_for_plot, ["PM25"], {"PM25":"ug/m3"}, "R%s_0==nan"%radius, bounds, mesh, data, mini, "H", colorbar, direc="Trial_Factors/", filter="")
        #"""

        #############
        ## Obtain corrected Model data
        print("      ...Obtaining new ", var, " values  ", flush=True)

        corrected=grid_factors.mul(corrected, axis=0, fill_value=1)

    ## Removes the days that werent corrected fully
    first_date=corrected.index.get_level_values(0).values[0]
    last_date=corrected.index.get_level_values(0).values[-1]

    last_date=datetime.strptime(str(last_date).split(":")[0], '%Y-%m-%dT%H')-timedelta(hours=window)
    first_date=datetime.strptime(str(first_date).split(":")[0], '%Y-%m-%dT%H')+timedelta(hours=window)


    return corrected.sort_index().loc[first_date.strftime("%Y-%m-%d %H:%M:%S"):last_date.strftime("%Y-%m-%d %H:%M:%S"), :, :]


######################################################
def main(t_ObsStationData, t_ModelStationData, t_ModelData, t_Stations_Info, t_ModelStation_Info, t_Models_Info, window=10, radius=10):
 
    if IncF.i_BC==0:
        return {}, {}

    print("")
    print("")
    print(".................Bias Correction................")
    
    ######################
    ## Find variables
    models=[c for c in IncF.c_ModelNames if c]
    obs=[c for c in IncF.c_ObsNetwork if c]

    if len(models)==0:
        print("...No models were read. Skipping this step")
        return {}, {}

    variables=aux.find_equiv_vars(IncF.c_BC_Vars)
    variables=[(i,x) for i in variables["O"] for x in variables["O"][i]]

    mods_asked=[models[int(m.split("_")[1])-1] for m in IncF.c_BC_Models]

    ####################
    new_models={}
    new_info={}
    for freq, window in zip(IncF.c_BC_freq, IncF.i_BC_window):
        new_models[freq]={}
        for mod, new_mod in zip(mods_asked, IncF.c_BC_NewNames):
            corrected_mods=[]
            for net, var in variables:
                print(f"...Correcting model  {mod}  for var  {var}", flush=True)
                data_obs=deepcopy(t_ObsStationData[freq][net])
                data_model=deepcopy(t_ModelData[freq][mod])
                data_model.columns=[c.replace("N_CO", "CO") for c in data_model.columns]
                mesh=t_Models_Info[mod]
                stations_info=t_Stations_Info.loc[net]
                models_info=t_ModelStation_Info[mod].loc[net]
                ##
                New_model=corrected_model(stations_info, models_info, data_obs, data_model, mesh, var, mod, mean_type=mean_type, window=window)
                corrected_mods.append(New_model)
            ##
            corrected=pd.concat(corrected_mods, axis=1)
            corrected.meta.units=t_ModelData[freq][mod].meta.units
            new_models[freq][new_mod]=corrected
            new_info[new_mod]=t_Models_Info[mod]

    ####
    ## Save to file
    print(" ")
    write_netcdf.main(new_models, new_info, {mod:[v[1] for v in variables] for mod in IncF.c_BC_NewNames}, "H,D", direc=IncF.c_BC_Dir)

    

    return new_models,  new_info
