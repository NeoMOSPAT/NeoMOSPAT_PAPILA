# -*- coding: utf-8 -*-
from os.path import isfile, join

import geopandas as gpd
import numpy as np
import pandas as pd

import src.IncludeFile as IncF
import src.common as aux
from src.read import read_shps
from copy import deepcopy
from shapely.geometry import Polygon
import pygeos

######################################################
def get_intersection_model(shape, mesh_or):
    ## Get bounds and smaller grid
    mesh=mesh_or.copy()
    dx = max(mesh["Cent_lon"].diff().mode().values[0], mesh["Cent_lat"].diff().mode().values[0])
    #
    minx, miny, maxx, maxy = shape.buffer(2*dx).total_bounds
    mesh = mesh.loc[(mesh['Cent_lat'].between(miny, maxy))&(mesh['Cent_lon'].between(minx, maxx))]
    mesh["Area_grid"]=mesh.area
    shape["Area_poly"]=shape.area


    ## Intersects both shapefiles
    intersection=gpd.overlay(shape, mesh, how='intersection')  
    intersection["Area_part"]=intersection.area
    intersection["ratio"]=intersection["Area_part"]/intersection["Area_poly"]
    intersection["ratio_grid"]=intersection["Area_part"]/intersection["Area_grid"]  
    return intersection[["indexes", "ratio", "ratio_grid"]]


######################################################
def transform_to_grid_vec(mesh_or, origin=""):
    mesh=deepcopy(mesh_or)
    ######
    ## Reorder to make sure the resulting ordering will be correct
    mesh.sort_values(by=["i_lat", "i_lon"], inplace=True)
    ## Get max values
    xmax=mesh["i_lon"].max()
    ymax=mesh["i_lat"].max()
    n=xmax*ymax
    ## Setup the indexers 
    left = lower = slice(None, -1)
    upper = right = slice(1, None)
    corners = [[lower, left], [lower, right], [upper, right], [upper, left]]
    

    #######
    ## Get matrix data
    x=np.reshape(mesh["lon"].values, (-1, xmax+1))
    y=np.reshape(mesh["lat"].values, (-1, xmax+1)) 
    ## Allocate output array
    xy = np.empty((n, 4, 2))
    ## Get the cornes of each polygon
    for i, (rows, cols) in enumerate(corners):
        xy[:, i,  0] = x[rows, cols].ravel()
        xy[:, i, 1] = y[rows, cols].ravel()


    ## Retrieve settings to indentify the corner that the information was at.
    if origin=="bottom-left":
        origin, origin_x, origin_y=0, 0, 0
    elif origin=="bottom-right":
        origin, origin_x, origin_y=1, 1, 0
    elif origin=="top-right":
        origin, origin_x, origin_y=2, 1, 1
    elif origin=="top-left":
        origin, origin_x, origin_y=3, 0, 1
    else:
        origin, origin_x, origin_y=0, 0, 0
    ## Get i_lat and i_lon vals
    ilon = np.tile(np.arange(origin_x, origin_x+xmax), ymax)
    ilat = np.repeat(np.arange(origin_y, origin_y+ymax), xmax)

    ## Create geodataframe and plot result
    mesh_geometry = pygeos.creation.polygons(xy) 
    mesh_gdf = gpd.GeoDataFrame({"i_lon":ilon, "i_lat":ilat, "indexes":list(zip(ilat, ilon))}, geometry=mesh_geometry, crs="epsg:4326")
    mesh_gdf[["Cent_lon", "Cent_lat"]] = gpd.GeoDataFrame(geometry = mesh_gdf.boundary).apply(lambda x: x.iloc[0].coords[origin], axis=1, result_type="expand")
    mesh_gdf["Area_model"] = mesh_gdf.area

    return mesh_gdf


############################
## Apply filters of models
def applyfilters_models(t_ModelInfo):
    ## List of integers     
    list_filters=IncF.i_ModelFilters

    ### If filters are off then return all of the information
    t_filters={}
    """
    if not list_filters:
        print("   ...Filters are off. All grid is being stored.")
        tmp_data=deepcopy(t_ModelInfo)
        for model in t_ModelInfo:
            tmp_data_model=tmp_data[model]
            tmp_data_model["On"]=1
            tmp_data_model["ratio"]=1
            t_ModelInfo[model]=tmp_data_model
        return {model:{"f0":{"All":t_ModelInfo[model][["indexes", "ratio"]]}} for model in t_ModelInfo}, t_ModelInfo
    #"""

    BoolChile=False
    mesh_old=None
    for model in t_ModelInfo:
        ## Mesh of model
        mesh=t_ModelInfo[model]

        ## Check if same mesh has been seen before. Could be large improvement in time if meshes are very fine. 
        BoolSameMesh=False
        for model_old in t_filters:
            if mesh.equals(t_ModelInfo[model_old]):
                BoolSameMesh=True
                print("   ...Grid of model  %s  is the same as  %s . Reusing  %s  filters"%(model, model_old, model_old))
                t_filters[model]=deepcopy(t_filters[model_old])
                break
        if BoolSameMesh:
            continue


        ################################
        ## Apply filters
        t_filters[model]={}
        mesh_vec=transform_to_grid_vec(mesh)

        
        #####################   
        if IncF.b_OnlyLand:
            print("   ...Keeping only data in land")
            world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres")).drop(['gdp_md_est'], axis=1)
            xmin, ymin, xmax, ymax = t_ModelInfo[model].total_bounds
            world = world.cx[xmin:xmax, ymin:ymax]
            t_filters[model]["f0"]={}
            inter=get_intersection_model(world, mesh_vec)
            t_filters[model]["f0"]["OnlyLand"]=inter


        ##### 
        if 1 in list_filters:
            t_filters[model]["f1"]={}
            ## Mask filter
            masks = IncF.c_RegionMask # _Mod[index]
            alias = IncF.c_RegionMaskAlias 
            print("   ...Filter  1:   Region Mask   for ", model)
            l_area = []
            for al, area in zip(alias, masks):
                if isinstance(area, str):
                    ## Read geojson
                    gdf = gpd.read_file(area, crs="epsg:4326")
                    
                    ## Get indexes
                    inter=get_intersection_model(gdf, mesh_vec)
                    t_filters[model]["f1"][al]=inter

        ##### 
        if 2 in list_filters:
            t_filters[model]["f2"]={}
            ## Geographical area
            squares=IncF.i_SqrRegions # _Mod
            alias=IncF.c_SqrRegionsAlias # _Mod
            squares=[Polygon([[pol[2], pol[0]], [pol[3], pol[0]], [pol[3], pol[1]], [pol[2], pol[1]]]) for pol in IncF.i_SqrRegions]
            print("   ...Filter  2:   Square geographical areas ", squares, " for model ", model)
            l_squares=[]
            for al, square in zip(alias, squares):
                gdf = gpd.GeoDataFrame({"geometry":[square], "Alias":[al]}, geometry='geometry', crs="epsg:4326")
                
                ## Get indexes
                inter=get_intersection_model(gdf, mesh_vec)
                t_filters[model]["f2"][al]=inter


        ##### 
        if 4 in list_filters:
            t_filters[model]["f4"]={}
            regions=IncF.c_Regions #_Mod
            regions_all=list(set(sum(regions, [])))
            alias=IncF.c_RegionsAlias
            if len(regions)!=0:
                print("   ...Filter  4:   By Chilean region ", regions, " for model ", model)
                BoolChile=True
                ## Get 
                regions_shp=read_shps.read_main("Regiones", "epsg:4326")
                regions_shp["codregion"]=regions_shp["codregion"].astype(str)
                regions_shp=regions_shp[regions_shp["codregion"].isin([str(c) for c in regions_all])]

                for al, reg in zip(alias, regions):

                    ## Get indexes
                    inter=get_intersection_model(regions_shp[regions_shp["codregion"].isin(reg)], mesh_vec)
                    t_filters[model]["f4"][al]=inter

       
        ##### 
        if 5 in list_filters:
            t_filters[model]["f5"]={}
            comunas=IncF.c_Comunas # 
            alias=IncF.c_ComunasAlias
            if len(comunas)!=0:
                print("   ...Filter  5:   By Chilean comuna ", comunas, " for model ", model)
                BoolChile=True
                comunas_all=list(set(sum(comunas, [])))
                comunas_shp=read_shps.read_main("Comunas", "epsg:4326")
                comunas_shp=comunas_shp[comunas_shp["cod_comuna"].isin(comunas_all)]
                for al, com in zip(alias, comunas):
                    ## Get indexes
                    inter=get_intersection_model(comunas_shp[comunas_shp["cod_comuna"].isin(com)], mesh_vec)
                    t_filters[model]["f5"][al]=inter







    ###############
    ### Keep only indexes of interest
    t_ModelInfo_copy=deepcopy(t_ModelInfo)
    ## Check that it wont fail if onlyland is false and no filters are asked
    for model in t_ModelInfo:
        indxs=[list(t_filters[model][f][al]["indexes"].values) for f in t_filters[model] for al in t_filters[model][f]]
        list_indexes=list(set(sum(indxs, [])))
        mesh=t_ModelInfo_copy[model]
        #mesh=mesh[mesh["indexes"].isin(list_indexes)]
        mesh.loc[mesh["indexes"].isin(list_indexes), "On"]=1
        mesh=mesh.replace(np.nan, 0)
        t_ModelInfo_copy[model]=mesh


    if not list_filters and len(t_filters[model].keys())==0:
        print("   ...Filters are off. All grid is being stored.")
        tmp_data=deepcopy(t_ModelInfo)
        for model in t_ModelInfo:
            tmp_data_model=tmp_data[model]
            tmp_data_model["On"]=1
            tmp_data_model["ratio"]=1
            t_ModelInfo[model]=tmp_data_model
        return {model:{"f0":{"All":t_ModelInfo[model][["indexes", "ratio"]]}} for model in t_ModelInfo}, t_ModelInfo




    return t_filters, t_ModelInfo_copy


######################################################
### Applies filters:.
def filter_data(t_ModelInfo, t_ModelData):#, t_Stations_Info):
    dirs_shp=aux.read_yaml("shp")   

    ###########################################
    ## Filter meshes by the general filters
    t_Filters_vs, t_ModelInfo_vs = {}, {}#applyfilters_stations(t_ModelInfo, t_Stations_Info)
    t_Filters_vm, t_ModelInfo_vm = applyfilters_models(t_ModelInfo)
    t_Filters={**t_Filters_vs, **t_Filters_vm}

    if not t_Filters:
        return {}, t_ModelInfo, t_ModelData

    ## mesh only with indexes that will be used
    t_ModelInfo_vm_on={model:t_ModelInfo_vm[model][t_ModelInfo_vm[model]["On"]==1] for model in t_ModelInfo_vm}

    ## No changes asked to be made
    if not IncF.i_ModelFilters:
        return t_Filters, t_ModelInfo_vm, t_ModelData

    ###############
    bool_vs=True if len(t_ModelInfo_vs)==0 else pd.concat(list(t_ModelInfo_vs.values())).empty #pd.concat(list(t_ModelInfo_vs.values())).empty
    bool_vm=pd.concat(list(t_ModelInfo_vm.values())).empty
    if not bool_vs and not bool_vm:
        indexes={k:np.unique(np.concatenate([t_ModelInfo_vs[k]["indexes"].values, t_ModelInfo_vm_on[k]["indexes"].values], axis=0)) for k in t_ModelInfo}
    elif bool_vs and not bool_vm:
        indexes={k:np.unique(t_ModelInfo_vm_on[k]["indexes"].values) for k in t_ModelInfo}
    elif not bool_vs and bool_vm:
        indexes={k:np.unique(t_ModelInfo_vs[k]["indexes"].values) for k in t_ModelInfo}
    else:
        indexes={k:[] for k in t_ModelInfo}

    ## TODO: Compare efficiency of this compared to concating and deleting duplicates
    t_ModelInfo={k:t_ModelInfo[k][t_ModelInfo[k]["indexes"].isin(indexes[k])] for k in t_ModelInfo}

    ## TODO: This needs to be done better. Not catching the empty ones yet
    if pd.concat(list(t_ModelInfo.values())).empty:
        print("...Specific r kequirements for the data to be read result in no data.")
        return {}, t_ModelInfo, t_ModelData


    ###########################################
    ## Get filtered data
    for model in t_ModelInfo:
        mesh=t_ModelInfo[model]
        ##
        data=t_ModelData[model]
        dims=list(data.index.names)
        ## t_dim for 2D data and t dim and height 3D data
        t_dim=dims[:-2]
        data=data.reset_index(level=t_dim).loc[mesh["indexes"]].set_index(t_dim, append=True) 
        data=data.reorder_levels(dims).sort_index()
        data[t_dim[0]]=data.index.get_level_values(0)
        t_ModelData[model]=data


    return t_Filters, t_ModelInfo_vm, t_ModelData

