# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
#
import src.common as aux
from src.read import read_shps
#
from os import makedirs
#
import matplotlib.pyplot as plt
import pygeos
import timeit

list_mappings=[int(c) for c in str(IncF.i_RMap)]
######################################################
def get_intersection_model(models, shape, variables, data_mod, info_mod):
    def my_agg(x, cols):
        for col in cols:
            x[col]=x["ratio"]*x[col]
        x=x[cols].sum()
        return x
    
    result={}
    for count, model in enumerate(models):
        mesh_new=info_mod[model]
        if count!=0 and mesh_new.equals(mesh):
            result[model]=df
            continue
        mesh=mesh_new
        data=data_mod[model].reset_index()
        data["indexes"]=data.apply(lambda x: (x.i_lat, x.i_lon), axis=1)
        var=variables[model]

        ## Intersects both shapefiles
        intersection=gpd.overlay(shape, mesh, how='intersection')  
        intersection["Area_part"]=intersection.area
        idd="ID"#IncF.c_RMap_ID
        

        ## Geometry Shapefile
        if 1 in list_mappings:
            #inter_sjoi=gpd.sjoin(mesh, shape, how='right', op='intersects') 
            intersection["ratio"]=intersection["Area_part"]/intersection["Area_poly"] 
            data_inter=intersection.merge(data, on="indexes")
            df=data_inter.groupby(["time", idd]).apply(my_agg, var).reset_index()

            df=df.merge(shape[[idd, "geometry"]], on=idd) 
            df=gpd.GeoDataFrame(df[[c for c in df.columns if c!="geometry" ]], geometry=df["geometry"], crs="epsg:4326")
            result[model+"_average"]=df


        ## Geometry model
        if 2 in list_mappings:
            #intersection=gpd.overlay(shape, mesh[mesh["indexes"].isin(intersection["indexes"])], how='union').replace(np.nan, 0)
            #intersection["Area_part"]=intersection.area
            intersection["ratio"]=intersection["Area_part"]/intersection["Area_model"]

            data_inter=intersection.merge(data, on="indexes")
            df=data_inter.groupby(["time", idd]).apply(my_agg, var).reset_index()

            df=df.merge(shape[[idd, "geometry"]], on=idd) 
            df=gpd.GeoDataFrame(df[[c for c in df.columns if c!="geometry" ]], geometry=df["geometry"], crs="epsg:4326")

            result[model+"_sum"]=df

        #"""
        ## Plot
        data_geo=data.merge(mesh[["indexes", "geometry"]], on="indexes") 
        data_geo=gpd.GeoDataFrame(data_geo[[c for c in data_geo.columns if c!="geometry" ]], geometry=data_geo["geometry"], crs="epsg:4326")
        data_geo=data_geo[data_geo["indexes"].isin(intersection["indexes"])]

        ax1=data_geo[data_geo["time"]=="2020-07-01"].plot(column="PM25", cmap="viridis", vmin=0, vmax=35)
        shape.boundary.plot(ax=ax1, color="black")
        ax1.axis('equal')
        ax1.axis("off")


        ax2=df[df["time"]=="2020-07-01"].plot(column="PM25", cmap="viridis", vmin=0, vmax=35)
        ax2.axis('equal')
        ax2.axis("off")

        plt.savefig("trial.png")

        #breakpoint()   
        #""" 


    return result


######################################################
def transform_to_grid_vec(mesh, origin=""):
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



######################################################
def main(t_ObsStationData,  t_ModelData, t_Stations_Info, t_Models_Info):
 
    if IncF.i_RMap==0:
        return None

    print("")
    print("")
    print("...................Re Mapping...................")
    
    ######################
    ## Find variables
    models=[c for c in IncF.c_ModelNames if c]
    obs=[c for c in IncF.c_ObsNetwork if c]
    for freq, var in zip(IncF.c_RMap_freq, IncF.c_RMap_data):
        print(f"   ...Re-mapping model  {var}  with freq  {freq}", flush=True)
        variables=aux.find_equiv_vars([var])

        ######################
        ## Retrieve data
        data_obs = pd.DataFrame()
        if "O" in variables:
            data_obs=pd.concat([t_ObsStationData[freq][m][variables["M"][m]] for m in variables["O"]], keys=variables["O"].keys(), axis=1)

        data_mod = pd.DataFrame()
        info_mod = pd.DataFrame()
        if "M" in variables:
            print("   ...Transforming point grid of models requested to polygon grid", flush=True)
            data_mod=pd.concat([t_ModelData[freq][m][variables["M"][m]] for m in variables["M"]], keys=variables["M"].keys(), axis=1)
            info_mod={m:transform_to_grid_vec(t_Models_Info[m], origin="bottom-left") for m in variables["M"]}



        ######################
        ## Read Shapefile
        print(f"   ...Reading shapefile  {IncF.c_RMap_dir}", flush=True)
        shapefile=read_shps.read_shapefile(IncF.c_RMap_dir, IncF.c_RMap_crs, epsg_output=4326)
        if len(IncF.c_RMap_cols)==0:
            shape=shapefile
        else:
            shape=shapefile[IncF.c_RMap_cols]
        ids=shape[["geometry"]].drop_duplicates(keep="first")
        ids["ID"]=range(ids.shape[0])
        shape=shape.merge(ids, on="geometry")
        shape["Area_poly"]=shape.area

       
        ######################
        ## Find intersection
        print("   ...Finding intersections and obtaining averages", flush=True)
        mapped_df=get_intersection_model(list(variables["M"].keys()), shape, variables["M"], data_mod, info_mod)


        ######################
        ## Save models
        for model in mapped_df:
            basename=IncF.c_RMap_dir.split("/")[-1].replace(".shp", "_%s-%s"%(model, freq))
            direc=IncF.c_RMap_savedir

            shp_to_save=mapped_df[model].merge(shapefile, on=[c for c in shapefile.columns if c in mapped_df[model].columns]).to_crs(IncF.c_RMap_crs)
            shp_to_save["time"]=shp_to_save["time"].dt.strftime("%Y-%m-%d %H:%M:%S")
            #######
            if 1 in IncF.c_RMap_savetype:
                ## Add data of original shapefile

                makedirs(f"{direc}/Data/shp/{basename}", exist_ok=True)
                shp_to_save.to_file("%s/Data/shp/%s/%s.shp"%(direc, basename, basename))
                print(f"      ...Saving to: ", f"{direc}/Data/shp/{basename}/Stations.csv", flush=True)
            #######
            if 2 in IncF.c_RMap_savetype:
                for group, df_group in shp_to_save.drop(columns=["geometry"]).groupby(IncF.c_RMap_ID):
                    df_group[mapped_df[model].drop(columns=["ID", "geometry"]).columns].to_csv(f"{direc}/Data/Models-shp/{basename}/ID-{group}_{freq}.csv", index=False)
                    print(f"      ...Saving to: ", f"{direc}/Data/Models-shp/{basename}/ID-{group}_{freq}.csv", flush=True)

                cent=shp_to_save["geometry"].centroid
                shp_to_save["Longitud"]=cent.x
                shp_to_save["Latitud"]=cent.y

                info=shp_to_save[["ID", "Latitud", "Longitud"]].drop_duplicates(keep="first")
                info.columns=["ID", "Latitud", "Longitud"]
                makedirs(f"{direc}/Data/Models-shp/{basename}", exist_ok=True)
                info.to_csv(f"{direc}/Data/Models-shp/{basename}/Stations.csv", index=False)
                print(f"      ...Saving to: ", f"{direc}/Data/Models-shp/{basename}/Stations.csv", flush=True)
            
