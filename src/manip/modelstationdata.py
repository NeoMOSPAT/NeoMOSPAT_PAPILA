# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
import src.common as aux
#
import pandas as pd
import geopandas as gpd
#
from shapely.ops import nearest_points
from shapely import wkt
import geopy.distance
from scipy.spatial import cKDTree
import numpy as np
from src.plot.ticker_formatters import LatDMS, LonDMS



######################################################
## Finds the closest grid points of stations in df to the models
def find_closest_points(df, meshh, t_Model_Filters, neighbor=1):
    ## Find closest grid points of model for the stations
    mini=df.copy()
    mini=gpd.GeoDataFrame(mini, geometry=gpd.points_from_xy(mini.Longitud, mini.Latitud), crs="epsg:4326")

    ## Remove filter rows and vectorize model and obtain average for filters for those rows
    minif=mini.copy()
    minif=minif[minif["ID"].str.contains('^f[0-9]-+')]


    mini=mini[~mini["ID"].isin(minif["ID"])]
    

    #########
    # Get closest point to stations
    mesh=meshh.copy()
    mesh=mesh[mesh["On"]==1]

    def ckd_nearest(gdA, gdB):
        nA = np.array(list(gdA.geometry.apply(lambda x: (x.x, x.y))))
        nB = np.array(list(gdB.geometry.apply(lambda x: (x.x, x.y))))
        btree = cKDTree(nB)
        dist, idx = btree.query(nA, k=1)
        gdB_nearest = gdB.iloc[idx][["geometry"]].reset_index(drop=True)
        gdB_nearest.columns = ["Nearest"]
        gdf = pd.concat([gdA.reset_index(drop=True),  gdB_nearest], axis=1)
        return gdf

    mini = ckd_nearest(mini, mesh)
    mini["Distance"]=mini.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.Nearest.y, a.Nearest.x)).km, axis=1)


    #########
    ## Get indexes information
    mini["str_Nearest"]=mini['Nearest'].apply(wkt.dumps)
    mesh["str_geometry"]=mesh["geometry"].apply(wkt.dumps) 
    
    mesh=pd.DataFrame(mesh).drop(columns=["geometry", "lat", "lon"])
    mini=mini.merge(mesh, how="left", left_on="str_Nearest", right_on="str_geometry").drop(columns=["str_Nearest","str_geometry"])

    #########
    ## Final version
    mini=gpd.GeoDataFrame(pd.DataFrame(mini), geometry=mini["Nearest"]).drop(columns=["Nearest"])

    return mini


######################################################
## Retrieves the information on the nearest grid points specified in df
def retrieve_model_date(df, t_ModelData, t_Model_Filters):
    pairs=df[["i_lat", "i_lon"]].drop_duplicates().values 
    pairs=list([tuple(pair) for pair in pairs])

    ###########
    ## Retrieve data
    ## Find dims that exist
    dims=list(t_ModelData.index.names)
    ## t_dim for 2D data and t dim and height 3D data
    t_dim=[c for c in dims if "t" ==c[0] or "T"==c[0]]
    ## Drop time columns
    Station_Data=t_ModelData.drop(columns=[t_dim[0]] if t_dim[0] in t_ModelData.columns else [])

    ## Output same style if df is empty
    if df.empty:
        return pd.DataFrame(columns=list(Station_Data.columns)+list(Station_Data.index.names)).set_index(list(Station_Data.index.names))
    ## Remove time and height(if 3D) columns from index, then retrieve only lat lons of interest
    Station_Data=Station_Data.reset_index(level=t_dim)#.loc[pairs]


    minif=df.copy()
    minif=minif[minif["ID"].str.contains('^f[0-9]-+')]


    ## Find the grid per station 
    data_list=[]
    keys=[]
    for ID in df[~df["ID"].isin(minif["ID"])]["ID"]:
        df_ID=df[df["ID"]==ID]
        ## Retrieve data for that station only
        pair=tuple(df_ID[["i_lat", "i_lon"]].values[0])
        pair=[pair for i in range(1) if pair in Station_Data.index]
        data=Station_Data.loc[pair]
        ## Reorder index to have same order as before
        data=data.set_index(t_dim).sort_index()
        ## Add time column back again 
        #data[t_dim[0]]=data.index.get_level_values(0)
        
        data_list.append(data)
        keys.append(ID)


    ## Find the grid per filter
    for f in t_Model_Filters:
        print(f)
        filters=t_Model_Filters[f]
        if f=="f-OnlyLand" or f=="f0":
            continue
        for group in filters:
            keys.append(f"{f}-{group}")

            group_partition=t_Model_Filters[f][group]
            print(group_partition)
            total_ratio=group_partition["ratio"].sum()
            print(total_ratio)
            #Station_Data.loc[[pair]]
            indexes=[c for c in group_partition["indexes"].values if c in Station_Data.index]
            data=Station_Data.loc[indexes]
            data["indexes"]=data.index.values 
            data=data.merge(group_partition, on="indexes")
            variables=[c for c in data.columns if c not in ["time", "indexes", "ratio"]]
            
            variables_sum = [v for v in variables if "EMI" in v]
            variables_mean = [v for v in variables if v not in variables_sum]

            data[variables_mean]=data[variables_mean].mul(data["ratio"], axis="index")
            data[variables_sum]=data[variables_sum].mul(data["ratio_grid"], axis="index")
            #ff = {**{v:"sum" for v in variables_sum}, **{v:"mean" for v in variables_mean}}
            data=data[["time"]+variables].groupby("time").sum()#.agg(ff) 
            data[variables_mean]=data[variables_mean]/total_ratio
            data_list.append(data)

    #{'A':['sum','mean'], 'B':['prod'], 'D': lambda g: df.loc[g.index].E.sum()}
    my_data=pd.concat(data_list, axis=1, keys=keys).swaplevel(axis=1).sort_index(axis=1)
    my_data.meta.units=t_ModelData.meta.units

    return my_data


######################################################
## Creates the structure that MOSPAT uses for the data
## Same structure as t_ObsStationData
def main(stations_info, t_ObsStationData, models_info, t_ModelData, t_Model_Filters):
    #######
    if (len([c for c in IncF.c_ObsNetwork if c])==0) or (stations_info.empty):
        #if IncF.i_ModelFilters==0:
        print("...No stations were read. Skipping this step.")
        return {}, {k:{} for k in IncF.c_TimeMeanRes}
    if (len([c for c in IncF.c_ModelNames if c])==0) or (pd.concat(list(models_info.values())).empty):
        print("...No models were read. Skipping this step.")
        return {}, {k:{} for k in IncF.c_TimeMeanRes}
    #if IncF.i_TS==0 and IncF.i_TS_stats==0 and IncF.i_TS_stats==0:
    #    print("...No TS was asked. Skipping this step.")
    #    return {}, {k:{} for k in IncF.c_TimeMeanRes}

    ## In case a model was not found or network had no data given the parameters asked
    models=list(t_ModelData[list(t_ModelData.keys())[0]].keys())
    obs=list(t_ObsStationData[list(t_ObsStationData.keys())[0]].keys())   

    ########
    ## Find nearest grid point
    print("...Finding closest grid points to models")
    nearest={}
    for index, model in enumerate(models):
        mesh=models_info[model]
        ## If next model shares old grid with one obtained previously then use that instead 
        BoolDo=True
        for old_model in list(models_info.keys())[:index]:
            if mesh.equals(models_info[old_model]):
                nearest[model]=nearest[old_model]
                BoolDo=False
                print("   ...Reusing nearest locations from model  %s  for model  %s "%(old_model, model))
                break

        if BoolDo:
            nearest[model]={}
            ## Go through Observation network
            for group in obs:
                ## Find closest grid point to station
                df=find_closest_points(stations_info.loc[group], mesh, t_Model_Filters)
                minif=stations_info.loc[group].copy()
                minif=minif[minif["ID"].str.contains('^f[0-9]-+')]
                df=pd.concat([df, minif])
                df["location"]=df.apply(lambda  x: "%s %s"%(LatDMS(x.geometry.y).to_string(resolution='minutes'), LonDMS(x.geometry.x).to_string(resolution='minutes')), axis=1)
                nearest[model][group]=df
                print("   ...Done for model  %s  and network  %s "%(model, group))

    ## Summary Information
    my_info={model:pd.concat(list(nearest[model].values()), axis=0, keys=list(nearest[model].keys())) for model in models}


    ########
    ## Retrieve data
    print("...Retriving Model Station Data")
    #breakpoint()
    data={}
    for freq in list(t_ModelData.keys()): 
        data[freq]={}
        ## Go through models
        for model in models:
            data[freq][model]={}
            ## Go through Observation network
            for group in obs:
                ## Retrieve Information for grid point
                data[freq][model][group]=retrieve_model_date(nearest[model][group], t_ModelData[freq][model], t_Model_Filters[model])    
                print("   ...Done for model  %s  and network  %s  and freq  %s  "%(model, group, freq))


    return my_info, data
