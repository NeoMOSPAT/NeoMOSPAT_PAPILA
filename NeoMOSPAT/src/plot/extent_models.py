# -*- coding: utf-8 -*-
import pandas as pd
import geopandas as gpd
import numpy as np
import contextily as ctx
import matplotlib.pyplot as plt
import src.IncludeFile as IncF
from os.path import join
BoolInitialMap=IncF.i_InitialMap



############################################################
def main(t_Models_Info):

    if BoolInitialMap==1:
        polygons=[]
        labels=[]
        for model in t_Models_Info:
            ## Get outline of all polygons
            
            polygons.append(t_Models_Info[model].unary_union.convex_hull)
            #polygons.append(t_Models_Info[model].unary_union.boundary)

            labels.append(model)

        gdf=gpd.GeoDataFrame({"label":labels}, geometry=gpd.GeoSeries(polygons), crs="EPSG:4326")#.to_crs("EPSG:3857")
        ax=gdf.boundary.plot(color="black", figsize=(12, 12))
        ax.set_axis_off()
        ctx.add_basemap(ax=ax, crs="EPSG:4326")
        plt.savefig(join(IncF.c_FigureDir, 'Models-extent.pdf'))
        plt.close()

    