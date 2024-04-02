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
from src.manip import modelstationdata
from src.write  import write_netcdf 
from src import common as aux
from copy import deepcopy
#
import os
os.environ["PYTHONBREAKPOINT"]="pudb.set_trace"
from datetime import datetime, timedelta

mean_type="end"



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

    ##########
    ## Get distance between point and station
    if intersection.empty:
        intersection["distance_geodesic"]=None
        return intersection

    intersection["distance_geodesic"]=intersection.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.Station_coords.y, a.Station_coords.x)).km, axis=1 )
    intersection=intersection.sort_values(by="distance_geodesic")
    intersection["W"]=(radius**2-intersection["distance_geodesic"]**2)/(radius**2+intersection["distance_geodesic"]**2)

    return intersection


######################################################
def get_epsilon(mod, obs, var, type_ep="Cov/Var"):
    def different_maths(obs, mod, type_ep=type_ep):
        ## Turn columns to rows
        mod_rep = mod.reset_index().melt( id_vars='time')
        mod_rep.set_index("time", inplace=True)
        obs_rep = obs.reset_index().melt( id_vars='time')
        obs_rep.set_index("time", inplace=True)


        if len(obs_rep["value"].dropna().drop_duplicates(keep="first").values)<=1:
            return 1

        ## Get arrays
        obs_rep=obs_rep["value"].values
        mod_rep=mod_rep["value"].values



        ###########
        ## Variances
        var_obs=np.nanvar(obs_rep)
        var_mod=np.nanvar(mod_rep)
        var_error=np.nanvar(obs_rep-mod_rep)

        ## Standard Deviations
        std_obs=np.sqrt(var_obs)
        std_mod=np.sqrt(var_mod)

        ## Normalized Covariance
        matrix=np.array([obs_rep, mod_rep])
        cov=np.ma.cov(np.ma.masked_invalid(matrix), bias=1)
        cov.filled(np.nan)
        cov=abs(cov[0][1])

        ## Difference between Obs and Mod
        perc=abs(np.nanmean((obs_rep-mod_rep)/obs_rep))


        if type_ep=="Cov/Var":
            return cov/var_error

        elif type_ep=="Cov":
            return cov

        elif type_ep=="Norm-Cov":
            return cov/(std_obs*std_mod)

        elif type_ep=="Var-Err":
            return var_error/var_obs

        elif type_ep=="Var":
            return var_obs/var_mod

        else:
            return 1

    both=pd.concat([mod, obs], axis=1, keys=["mod", "obs"]) 
    both_R=both.rolling(window=24, min_periods=1, center=True)
   
    both["epsilon"]=list(map(lambda x: different_maths(x["obs"], x["mod"], type_ep=type_ep), both_R) )
    epsilon=both[["epsilon"]].fillna(method='backfill')
    epsilon.columns=[var]

    return  epsilon


######################################################
def interpolate_obs(df):
    def long_nan_series(series):
        # select this series when all values are NaNs
        all_nans = series.isnull().all()

        # and the delta is too long for interpolation
        # note: taking the last value minus the first,
        #       so this is the delta between the last NaN
        #       value and the first NaN value - there's
        #       an hour duration more until the next non-null value
        too_long = series.index[-1] - series.index[0] > pd.Timedelta("3H")
        return too_long & all_nans

    def get_average_value(series, mean_value, date):
        result=np.nan
        days=0       #(3, 6, 9, 12, 15, 18, 21, 24, 27, 30)  
        days_mean=-1 #(1, 3, 5,  7,  9, 11, 13, 15, 17, 19)
        while result!=result:
            days+=3
            ## If nothing is found within a month (15 days before, 15 days after) then stop
            if days>15:
                return np.nan
            ## Get last week
            timedelta=pd.Timedelta("%s days"%days)
            working_data=series.loc[date-timedelta:date+timedelta]
            ## Get only the ones that are the same hour
            working_data=working_data[working_data.index.hour==date.hour]

            result=working_data.mean()
            if result==result:
                return result

            ## Get surrounding three days
            days_mean+=2
            timedelta=pd.Timedelta("%s days"%days_mean)
            working_data=mean_value.loc[date-timedelta:date+timedelta]
            ## Get only the ones that are the same hour
            working_data=working_data[working_data.index.hour==date.hour]

            result=working_data.mean()
            if result==result:
                return result


    #mean_hours=df.groupby(df.index.hour).mean()
    mean_value=df.mean(axis=1)
    for col in df.columns:
        series=df[col]
        df_nan_group_keys = series.isnull().diff().ne(0).cumsum()
        series_long_nans = series.groupby(df_nan_group_keys).transform(long_nan_series)
        
        ## Small gaps 
        series_small_gaps=series[~series_long_nans]  
        series_small_interp = series_small_gaps.interpolate(method="linear")

        ## Only fills 3 hours, rest it leaves as nans
        df[col]=pd.concat([series_small_interp, series[series_long_nans]]).sort_index() 

        ## Interpolates all observations
        """
        ## Long gaps. Groups them by gaps
        series_long_gaps=series[series_long_nans]
        time_dif=series_long_gaps.index.to_series().diff()
        time_dif[time_dif>pd.Timedelta("1H")]=np.nan 
        time_dif=time_dif.replace(pd.Timedelta("1H"), 0).replace(np.nan, 1)  
        time_dif=time_dif.astype(int).cumsum()
        


        ## Retrieve each gap and find the new values
        both=pd.concat([series_long_gaps, time_dif], axis=1)
        both.columns=[col, "group"]
        series_long_interp=[]
        for group, df_group in both.groupby("group"):
            series_long_interp.append(df_group.apply(lambda x: get_average_value(series, mean_value, x.name), axis=1))
        series_long_interp=pd.concat(series_long_interp)

        df[col]=pd.concat([series_small_interp, series_long_interp]).sort_index() 
        """
    return df


######################################################
def SCM_estimate(mini, model, obs, modelstation, epsilon=1):
    stations=mini["ID"].values
    obs=obs[stations]
    var=epsilon.columns[0]

    ## Iteration
    corrected_model=model
    #for i in range(1): 
    ## Sum_esta W*(Obs_esta-Model_esta)
    sum_up=obs.apply(lambda x: np.nansum(mini["W"].values*(x.values-modelstation.loc[x.name].values)), axis=1).to_frame(name=var) 
    #sum_up.columns=model.columns

    ## Sum
    sum_down = obs.apply(lambda x: np.nansum(mini["W"].values*(x.values/x.values)), axis=1).to_frame(name=var) + epsilon*2 #mini["W"].sum()+epsilon*2
    ## New iteration value
    corrected_model=corrected_model+sum_up/sum_down
    

    return corrected_model


######################################################
## Find new corrected var values
def corrected_model(mini, mini_mod, data, model, mesh, var, model_name, modelstation, mean_type="end", window=1):
    ## Obtain points per station that are within radius of window
    corrected=model[[var]]
    for radius in IncF.i_BC_Rads:
        print("   ...Finding grid points inside radius of ", radius, flush=True)
        intersection=grids_within_radius(mini_mod, mesh, radius)
        #intersection[intersection["ID"]=="A01"][["indexes", "lat", "lon", "W", "distance_geodesic", "geometry"]]
        stations_all=intersection["ID"].drop_duplicates(keep="first").values        

        grouped=intersection.groupby(["i_lat", "i_lon"]) 
        for i in range(1):
            for group, df_group in grouped:
                stations=df_group["ID"].values
                ## Model
                mod_at_point = corrected.loc[:, group[0], group[1]]
                ## Obs
                obs_at_obs = data[var]
                obs_at_obs = obs_at_obs[df_group["ID"].values]
                #obs_at_obs = obs_at_obs.loc[mod_at_point.index]
                obs_at_obs = obs_at_obs.reindex(index = mod_at_point.index)
                obs_at_obs.dropna(how="all", axis=1, inplace=True)     
                if obs_at_obs.empty:
                    continue

                ## Checks again with stations have data
                stations = obs_at_obs.columns 

                ## Model Station 
                mod_at_obs = modelstation[var]
                mod_at_obs = mod_at_obs[stations]
                mod_at_obs = mod_at_obs.loc[mod_at_point.index]

                ## Df_group
                df_group=df_group[df_group["ID"].isin(stations)]

                ## Interpolation here
                obs_at_obs=interpolate_obs(obs_at_obs)

                ##
                epsilon=get_epsilon(mod_at_obs, obs_at_obs, var)
                new_model=SCM_estimate(df_group, mod_at_point, obs_at_obs, mod_at_obs, epsilon=epsilon)

                ## Replace Old values by the New values
                corrected.loc[(slice(None), group[0], group[1]), var]=new_model.values


            ## Get new values of model at station's location
            modelstation = modelstationdata.retrieve_model_date(mini_mod[mini_mod["ID"].isin(stations_all)], corrected)


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
    c_BC_NewNames=[c+"_SCM" for c in IncF.c_BC_NewNames]
    for freq, window in zip(IncF.c_BC_freq, IncF.i_BC_window):
        new_models[freq]={}
        for mod, new_mod in zip(mods_asked, c_BC_NewNames):
            corrected_mods=[]
            for net, var in variables:
                print(f"...Correcting model  {mod}  for var  {var}", flush=True)
                data_obs=deepcopy(t_ObsStationData[freq][net])
                data_model=deepcopy(t_ModelData[freq][mod])
                data_model.columns=[c.replace("N_CO", "CO") for c in data_model.columns]
                data_modelstation=deepcopy(t_ModelStationData[freq][mod][net])
                data_modelstation.columns=pd.MultiIndex.from_tuples([(c[0].replace("N_CO", "CO"), c[1]) for c in data_modelstation.columns])

                mesh=t_Models_Info[mod]
                stations_info=t_Stations_Info.loc[net]
                models_info=t_ModelStation_Info[mod].loc[net]
                ##

                New_model=corrected_model(stations_info, models_info, data_obs, data_model, mesh, var, mod, data_modelstation,  mean_type=mean_type, window=window)
                corrected_mods.append(New_model)
            ##
            corrected=pd.concat(corrected_mods, axis=1)
            corrected.meta.units=t_ModelData[freq][mod].meta.units
            new_models[freq][new_mod]=corrected
            new_info[new_mod]=t_Models_Info[mod]

    ####
    ## Save to file
    print(" ")
    write_netcdf.main(new_models, new_info, {mod:[v[1] for v in variables] for mod in c_BC_NewNames}, "H,D", direc=IncF.c_BC_Dir)

    

    return new_models,  new_info
