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
from src.plot import Surface_Maps as maps
from src.write  import write_netcdf 
from src import common as aux
from copy import deepcopy
#
import os
os.environ["PYTHONBREAKPOINT"]="pudb.set_trace"
from datetime import datetime, timedelta
from src.manip import modelstationdata
mean_type="end"
import matplotlib.pyplot as plt 
##TODO: Make units of CO be available
## TODO FIX sum up so that no error happens if dates are not shared 

######################################################
def interpolate_obs(df):
    def long_nan_series(series, window="3H"):
        # select this series when all values are NaNs
        all_nans = series.isnull().all()

        # and the delta is too long for interpolation
        # note: taking the last value minus the first,
        #       so this is the delta between the last NaN
        #       value and the first NaN value - there's
        #       an hour duration more until the next non-null value
        too_long = series.index[-1] - series.index[0] > pd.Timedelta(window)
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

        ## Remove outliers
        df.loc[df[col]>np.nanpercentile(df[col], 98), col] = np.nan

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
        #"""
    return df


######################################################
def get_factor_distance(distance, radius, method="cos"):
    ## Make sure it's one of the two
    if distance>=radius:
        return 0

    if method not in ["r2", "cos"]:
        method="cos"

    ## Return values
    if method =="r2":
        return (radius**2-distance**2)/(radius**2+distance**2)
    elif method == "cos":
        return 0.5-0.5*np.cos(np.pi*((distance-radius)/radius))


######################################################
def get_epsilon(mod, obs, var, type_ep="Cov/Var"):
    def different_maths(obs, mod, type_ep=type_ep):
        ## Turn columns to rows
        mod_rep = mod.reset_index().melt( id_vars='time')
        mod_rep.set_index("time", inplace=True)
        obs_rep = obs.reset_index().melt( id_vars='time').replace(0, np.nan)
        obs_rep.set_index("time", inplace=True)


        if len(obs_rep["value"].dropna().drop_duplicates(keep="first").values)<=1:
            return 1
        if len(mod_rep["value"].dropna().drop_duplicates(keep="first").values)<=1:
            return 1
        if len((mod_rep["value"]+obs_rep["value"]).dropna().drop_duplicates(keep="first").values)<=1:
            return 1

        
        ## Get arrays
        obs_rep=obs_rep["value"].values
        mod_rep=mod_rep["value"].values



        ###########
        ## Variances
        var_obs = np.nanvar(obs_rep)
        var_mod = np.nanvar(mod_rep)
        var_error = np.nanvar(obs_rep-mod_rep)
        var_error = var_error if var_error!=0 else 1

        ## Standard Deviations
        std_obs=np.sqrt(var_obs)
        std_mod=np.sqrt(var_mod)

        ## Normalized Covariance
        matrix=np.array([obs_rep, mod_rep])
        cov=np.ma.cov(np.ma.masked_invalid(matrix), bias=1)
        cov.filled(np.nan)
        cov=abs(cov[0][1])
        cov = cov if cov!=0 else 1

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
    if mod[mod.index.duplicated()].shape[0]!=0:
        
        mod=mod[~mod.index.duplicated()]  
    if obs[obs.index.duplicated()].shape[0]!=0:
        obs=obs[~obs.index.duplicated()]  
        

    both=pd.concat([mod, obs], axis=1, keys=["mod", "obs"]) 
    both_R=both.rolling(window=24, min_periods=1, center=True)
    both["epsilon"]=list(map(lambda x: different_maths(x["obs"], x["mod"], type_ep=type_ep), both_R) )
    epsilon=both[["epsilon"]].fillna(method='backfill')
    epsilon.columns=[var]

    return  epsilon


######################################################
def SCM_estimate(ID, mini, model, obs, modelstation, epsilon=0):
    stations=mini["IDs_within"].values
    obs=obs[stations]
    var=epsilon.columns[0]
    bothh=pd.concat([obs, modelstation], keys=["obs", "model"], axis=1)
    ## Iteration
    corrected_model=model
    #for i in range(1): 
    ## Sum_esta W*(Obs_esta-Model_esta)
    sum_up=bothh.apply(lambda x: np.nansum(mini["distance_factor"].values*mini["W"].values*(x["obs"].values-x["model"].values)), axis=1).to_frame(name=var)

    ## Sum
    sum_down = np.nansum(mini["W"].values) + epsilon*2 
    #sum_down =obs.apply(lambda x: np.nansum(mini["W"].values*(x.values/x.values)), axis=1).to_frame(name=var) + epsilon*2 
    ## New iteration value
    corrected_model=corrected_model+(sum_up/sum_down).replace(np.nan, 0)

    return corrected_model

######################################################
## Creates the bias factor between grid point ad station for a daily period and a window of 10 days of rolling mean.
def bias_station(mini, mod_data, obs_data, mod_at_obs, obs_at_obs, var, ID, mean_type="end", window=10):


    epsilons=get_epsilon(mod_at_obs, obs_at_obs, var)

    new_mod=SCM_estimate(ID, mini, mod_at_obs, obs_data, mod_data, epsilons)

    all_df=pd.concat([new_mod, mod_at_obs, obs_at_obs], keys=["new", "mod", "obs"], axis=1)
    factor_station=all_df["new"]/all_df["mod"]
    factor_station.columns=[ID]
    
    return factor_station.fillna(1)#all_df["new"]/all_df["mod"]#(all_df["new"]-all_df["mod"])/all_df["mo"]


#####################################################
## Obtain bias correction for all stations in mini with closest grid point 
def bias_all_stations(modelstation,  mini,  model, model_name, data, var, mean_type="end", window=10):

    ####################################################
    factors=[]
    ### Get data nearest to stations
    for ID in mini["ID"].drop_duplicates(keep="first"):
        info=mini.loc[mini["ID"]==ID]
        stations_inside=info["IDs_within"].values

        ## Find indexes to extract model data
        nearest=info[["i_lon", "i_lat"]].values[0]

        ## Multiple stations
        mod_data=modelstation[var][stations_inside]
        obs_data=data[var][stations_inside]

        ## At point (which is obs location)
        mod_at_obs=model.loc[:, nearest[1], nearest[0]][[var]]
        obs_at_obs=aux.get_station_data(data, ID)[[var]]
        
        ## Interpolation here
        obs_at_obs=interpolate_obs(obs_at_obs)



        ## Obtain bias factor correction
        factor_station=bias_station(info, mod_data, obs_data, mod_at_obs, obs_at_obs, var, ID, mean_type=mean_type, window=window)
        

        ## Add factor to dictionary
        factors.append(factor_station)

    factors=pd.concat(factors, axis=1)

    if factors.shape != factors.dropna(axis=1, how="all").shape:
        print("Some error that shouldnt happen")
        #breakpoint()

    return factors.dropna(axis=1, how="all").fillna(1)


######################################################
## Find grid points that are within radius of all stations
proj_4326 = pyproj.Proj('+proj=longlat +datum=WGS84')
def grids_within_radius(mini_or, mesh, radius, method_distance="cos"):
    mini=mini_or.copy()
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

            points=transform(project, buf).exterior
            #np.array(transform(project, buf).exterior)
            return np.array(list(zip(np.array(points.xy[0]), np.array(points.xy[1]))))

            ## Returns warning of .coords being deprecated in the future
            #return transform(project, buf).exterior.coords[:]

        points_list=geodesic_point_buffer(lat, lon, km)
        return Polygon(points_list)

    ####################################################
    ##########
    ## Add circular polygon rodenado station as geometry
    mini["geom_model"]=mini.geometry

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
        df["geom_model"]=data.geom_model
        return df

    intersection=[get_points_inside(row) for index, row in mini.iterrows()]
    intersection=pd.concat(intersection).sort_values(by=["i_lat", "i_lon"]).reset_index(drop=True)

    ##########
    ## Get distance between point and station
    if intersection.empty:
        intersection["distance_geodesic"]=None
        return intersection

    ## Get distance
    intersection["distance_geodesic"]=intersection.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.geom_model.y, a.geom_model.x)).km, axis=1 )

    ## Get W factor
    intersection["W"]=intersection.apply(lambda x: get_factor_distance(x.distance_geodesic, radius, method="R2"), axis=1)
    ## Get distance factor
    intersection["distance_factor"]=intersection.apply(lambda x: get_factor_distance(x.distance_geodesic, radius, method=method_distance), axis=1)

    return intersection.sort_values(by="distance_geodesic")



def stations_within_radius(mini_or, radius, method_distance="cos"):
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
            points=transform(project, buf).exterior
            #np.array(transform(project, buf).exterior)
            return np.array(list(zip(np.array(points.xy[0]), np.array(points.xy[1]))))
            #return transform(project, buf).exterior.coords[:]

        points_list=geodesic_point_buffer(lat, lon, km)
        return Polygon(points_list)

    ####################################################
    ##########
    ## Add circular polygon rodenado station as geometry
    mini=mini_or.copy()
    mini=gpd.GeoDataFrame(mini, crs="epsg:4326", geometry=gpd.points_from_xy(mini['Longitud'], mini['Latitud']))
    mini_new=mini.copy()
    geometry=mini.apply(lambda a: get_circle(a.Latitud, a.Longitud, radius), axis=1)
    mini_circle=gpd.GeoDataFrame(mini, crs='epsg:4326', geometry=geometry)
    mini_circle["Points_within"]=[Point(xy) for xy in zip(mini_circle.Longitud, mini_circle.Latitud)]
    mini_circle=mini_circle[["ID", "Points_within", "geometry"]]
    mini_circle.columns=["IDs_within", "Points_within", "geometry"]
    

    intersection=gpd.sjoin(mini_new[["ID", "Nombre", "indexes", "i_lat", "i_lon", "geometry"]], mini_circle, how="left",  predicate="within") 
    ##########
    ## Get distance between point and station
    if intersection.empty:
        intersection["distance_geodesic"]=None
        return intersection

    ## Fill missing ppints
    intersection["Points_within"].fillna(intersection.geometry, inplace=True)
    ## Get distance
    intersection["distance_geodesic"]=intersection.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.Points_within.y, a.Points_within.x)).km, axis=1 )
    ## Get W factor
    intersection["W"]=intersection.apply(lambda x: get_factor_distance(x.distance_geodesic, radius, method="R2"), axis=1) #(radius**2-intersection["distance_geodesic"]**2)/(radius**2+intersection["distance_geodesic"]**2)

    ## Get distance factor
    intersection["distance_factor"]=intersection.apply(lambda x: get_factor_distance(x.distance_geodesic, radius, method="cos"), axis=1)

    return intersection.sort_values(by="distance_geodesic")


######################################################
## Generate grid with all their respective biases. Considers distance, and multiple biases factors from multiple stations. 
def get_factor_grid(intersection, factors_stations, model, var, radius=10, mesh=None, mini=None, data=None):
    factor_info=intersection[["i_lat", "i_lon", "lat", "lon", "distance_geodesic", "ID"]]
    factor_info=factor_info.reset_index(drop=True)

    
    ####################################################  
    ## Model at points within radious      
    if factor_info.empty:
        #corrected=model.reset_index(level=0).loc[intersection.indexes].reset_index().set_index(["time", "i_lat", "i_lon"])[[var]] 
        #corrected[var]=1
        corrected=pd.DataFrame(columns=[var, "time", "i_lat", "i_lon"]).set_index(["time", "i_lat", "i_lon"]) 
        corrected[var]=pd.to_numeric(corrected[var], errors="coerce")
        return corrected.sort_index()
    
    ##########
    ## Find the perdentage with which it will go down. From factor at the top for full strength to 1 when there is no strength (no change in definitive field)
    ## Retuturn multiplier given the distance
    def line_eq(distance, rad=radius):
        distance=distance
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
    factor_info[["multiplier", "sum"]]=factor_info.apply(lambda a: line_eq(a.distance_geodesic), axis=1, result_type="expand")
    factor_info["distance_factor"]=factor_info.apply(lambda x: get_factor_distance(x.distance_geodesic, radius, method="R2"), axis=1)
    factor_info=factor_info.sort_values(by="distance_geodesic")[["i_lat", "i_lon", "lat", "lon", "multiplier",  "sum", "ID", "distance_geodesic", "distance_factor"]]


    ## Add data of the weights coming from each station for the points where more than one station is within range
    sum_factor_distance = factor_info.groupby(["i_lat", "i_lon"])[["distance_factor"]].sum().reset_index()
    sum_factor_distance.columns = ["i_lat", "i_lon", "sum_distance_factor"]

    factor_info = pd.merge(factor_info, sum_factor_distance, on=["i_lat", "i_lon"])
    factor_info["f/sum_f"]= factor_info["distance_factor"]/factor_info["sum_distance_factor"]
    


    ####################################################  
    ## Correct Model at points within radious  
    dfs=[]
    dfs_without=[]
    for name, factors_grid in factor_info.groupby("ID"):
        factor_id=factors_stations[[name]]

        factors_grid.reset_index(drop=True, inplace=True)
        ## change this lat lon for i_lat i_lon
        factors_grid.set_index(['i_lat', 'i_lon'], inplace=True)
        ## Add propotion of grid points that are coming from that specific station
        factors_grid["sum_d"]=factors_grid["sum"]*factors_grid["f/sum_f"]
        factors_grid["multiplier_d"]=factors_grid["multiplier"]*factors_grid["f/sum_f"]

        final_factor=pd.DataFrame((factor_id.values * factors_grid["multiplier_d"].values + factors_grid["sum_d"].values)*1, index=factor_id.index, columns=factors_grid.index).unstack().to_frame(name=var)

        dfs.append(final_factor)




    ############
    final_factor=pd.concat(dfs, axis=0).sort_index()
    final_factor=final_factor.reorder_levels(["time", "i_lat", "i_lon"])
    final_factor[var]=pd.to_numeric(final_factor[var], errors="coerce")

    ## Merge cells that have more than one station
    final_factor=final_factor.groupby(level=[0,1,2]).sum()


    return final_factor.replace(0, np.nan).fillna(1)




######################################################
## Find new final_factor var values
def corrected_model(mini, mini_mod, data, model, mesh, var, model_name, modelstation, mean_type="end", window=1):
    ## Obtain bias correction factor of all stations in mini
    print("   ...Finding bias factor of stations", flush=True)
    stations=stations_within_radius(mini_mod, 5, method_distance="cos")

    boolPM10=False
    ## Obtain points per station that are within radius of window
    corrected=model[[var]]

    initial_dates=corrected.index.get_level_values(0).drop_duplicates(keep="first").strftime("%Y-%m-%d %H:%M:%S")
    for radius in IncF.i_BC_Rads:
        
        ## Get new values of model at station's location
        modelstation = modelstationdata.retrieve_model_date(mini_mod, corrected)
        
        ## Get correction factor of station points
        factors_stations=bias_all_stations(modelstation,  stations, corrected, model_name, data, var, mean_type=mean_type, window=window)

        ##
        print("   ...Finding grid points inside radius of ", radius, flush=True)
        intersection=grids_within_radius(mini_mod, mesh, radius)


        #############
        ## Obtain grid factor for all grid points inside at least one circle of a station. If point has multiple stations in range, then the average of the factor of each one (counting distance) is computed.
        print("      ...Obtaining bias correction factor for grid points", flush=True)
        final_factor=get_factor_grid(intersection, factors_stations, model, var, radius=radius, mesh=mesh, mini=mini, data=data)
        ## Make it so that dates are only those that the model had initially
        final_factor=final_factor.iloc[final_factor.index.get_level_values(0).isin(initial_dates)]


        #############
        ## Obtain final_factor Model data
        print("      ...Obtaining new ", var, " values  ", flush=True)
        
        #final_factor.loc[initial_dates, slice(None), slice(None)]
        corrected=corrected.mul(final_factor, fill_value=1, axis=0)
        if boolPM10:
            corrected_pm10 = corrected_pm10.mul(final_factor.rename(columns={var:var_old}) , fill_value=1, axis=0)

        

    
    if boolPM10:
        corrected=corrected_pm10

    ## Removes the days that werent final_factor fully
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
    if len(variables)==0:
        print("...Requested vars weren't read. Skipping this step", flush=True)
        return {}, {}

    variables=[(i,x) for i in variables["O"] for x in variables["O"][i]]

    mods_asked=[models[int(m.split("_")[1])-1] for m in IncF.c_BC_Models]

    ####################
    new_models={}
    new_info={}
    unique_vars=list(set([g[1] for g in variables]))
    ## Arreglo trucho para poder corregir una variable que no se leyó en ninguna network
    #if "PM10" in unique_vars:
    #    unique_vars.append("PM10ant")
    #    #unique_vars.append("PM25ant")
    #    #unique_vars.append("PM25")
    #if "PM25" in unique_vars:
    #    unique_vars.append("PM25ant")

    #[g[0] for g in variables if g[1]=="PM25"] 
    c_BC_NewNames=[c+"_SCM-Factor" for c in IncF.c_BC_NewNames]
    for freq, window in zip(IncF.c_BC_freq, IncF.i_BC_window):
        new_models[freq]={}
        for mod, new_mod in zip(mods_asked, c_BC_NewNames):
            final_factor_mods=[]
            #for net, var in variables:
            for var in unique_vars:
                #"""
                if var in []:#["PM10", "PM10ant", "PM25ant"]:
                    nets=[g[0] for g in variables if g[1]=="PM25"]
                elif var=="PM10ant":
                    nets=[g[0] for g in variables if g[1]=="PM10"]
                else:
                    nets=[g[0] for g in variables if g[1]==var]
                #"""
                print(f"...Correcting model  {mod}  for var  {var}", flush=True)
                data_obs=[]
                data_modelstation=[]
                stations_info=[]
                models_info=[]
                for net in nets:
                    data_obs.append(deepcopy(t_ObsStationData[freq][net]))
                    data_modelstation.append(deepcopy(t_ModelStationData[freq][mod][net]))
                    stations_info.append(t_Stations_Info.loc[net])
                    models_info.append(t_ModelStation_Info[mod].loc[net])

                data_obs=pd.concat(data_obs, axis=1)
                data_modelstation=pd.concat(data_modelstation, axis=1)
                stations_info=pd.concat(stations_info)
                models_info=pd.concat(models_info)

                data_model=deepcopy(t_ModelData[freq][mod])
                data_model.columns=[c.replace("N_CO", "CO") for c in data_model.columns]  
                data_modelstation.columns=pd.MultiIndex.from_tuples([(c[0].replace("N_CO", "CO"), c[1]) for c in data_modelstation.columns])

                mesh=t_Models_Info[mod]
                New_model=corrected_model(stations_info, models_info, data_obs, data_model, mesh, var, mod, data_modelstation, mean_type=mean_type, window=window)
                final_factor_mods.append(New_model)
            ##
            final_factor=pd.concat(final_factor_mods, axis=1)
            final_factor.meta.units=t_ModelData[freq][mod].meta.units
            final_factor.meta.units={key.replace("N_CO", "CO"):value for key, value in final_factor.meta.units.items()}

            
            ## Add variables that were not asked to correct
            for var in t_ModelData[freq][mod].columns:
                if var not in final_factor.columns:
                    final_factor[var]=t_ModelData[freq][mod][var]

            new_models[freq][new_mod]=final_factor
            new_info[new_mod]=t_Models_Info[mod]

    ####
    ## Save to file
    print(" ")
    #write_netcdf.main(new_models, new_info, {mod:[v[1] for v in variables] for mod in c_BC_NewNames}, "H,D", direc=IncF.c_BC_Dir)
    write_netcdf.main(new_models, new_info, {mod:new_models["H"][mod].columns for mod in c_BC_NewNames}, "H,D", direc=IncF.c_BC_Dir)
    unique_vars

    

    return new_models,  new_info
