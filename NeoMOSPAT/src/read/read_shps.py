# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
import geopandas as gpd
from shapely.geometry import Polygon, Point

from src.common import remove_accents, read_yaml
    
##############################################
def read_shapefile(name, crs_input, epsg_output=4326):
    shp=gpd.read_file(name)
    shp.crs=crs_input
    shp=shp.applymap(remove_accents)
    shp=shp.set_crs(crs_input)
    if isinstance(epsg_output, str):
        epsg_output=epsg_output.replace("epsg:", "")
    shp=shp.to_crs(epsg=epsg_output)
    return shp


##############################################
def read_main(name, crs_out):
    dirs_shp=read_yaml("shp")
    return read_shapefile(dirs_shp["SHP"][name], dirs_shp["CRS"][name], crs_out)


##############################################
def read_shp_maps(crs_out="epsg:4326", n7=IncF.f_Map_Squares, n8=IncF.f_Map_Points):
    list_shp=[int(c) for c in str(IncF.i_Map_Shapefiles)]
    shapefiles={}
    if 1 in list_shp:
        shapefiles["Ciudades"]=read_main("Ciudades", crs_out)
    if 2 in list_shp:
        shapefiles["Comunas"]=read_main("Comunas", crs_out)
    if 3 in list_shp:
        shapefiles["Provincias"]=read_main("Provincias", crs_out)
    if 4 in list_shp:
        shapefiles["Regiones"]=read_main("Regiones", crs_out)
    if 5 in list_shp:
        shapefiles["Chile"]=read_main("Chile", crs_out)
    if 6 in list_shp:
        shapefiles["Topo_250m"]=read_main("Topo_250m", crs_out)
    if 7 in list_shp:
        squares=n7
        squares=[Polygon([[pol[2], pol[0]], [pol[3], pol[0]], [pol[3], pol[1]], [pol[2], pol[1]]]) for pol in squares]
        if len(squares)!=0:
            shapefiles["Squares"] = gpd.GeoDataFrame({"geometry":squares}, geometry = 'geometry', crs = "epsg:4326").to_crs(epsg = crs_out.split(":")[-1])
    if 8 in list_shp:
        points=n8
        points=[Point(p[1], p[0]) for p in points]
        if len(points)!=0:
            shapefiles["Points"] = gpd.GeoDataFrame({"geometry":points}, geometry = 'geometry', crs = "epsg:4326").to_crs(epsg = crs_out.split(":")[-1])
    if 9 in list_shp:
        shapefiles["World"]=gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))#.to_crs(epsg=crs_out)
        
    return shapefiles