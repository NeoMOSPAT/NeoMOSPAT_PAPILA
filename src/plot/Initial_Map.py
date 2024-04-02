# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
import src.common as aux
#
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from matplotlib import  cm
#
import numpy as np
import pandas as pd
import geopandas as gpd
from os.path import join
#
import contextily as ctx
import src.read.read_shps as shps
from shapely.geometry import Polygon
from adjustText import adjust_text


BoolInitialMap=IncF.i_InitialMap

####################################
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


#  Returns tuple of handles, labels for axis ax, after reordering them to conform to the label order `order`, and if unique is True, after removing entries with duplicate labels.
def reorderLegend(ax=None, order=None, characters=20):
    handles, labels = ax.get_legend_handles_labels()
    info = dict(zip(labels, handles))

    new_handles = [info[l] for l in order]
    #new_labels = [fill(l.replace("\n", "~"), characters).replace("~", "\n") for l in order]
    return new_handles, order


def get_label_stations(info_group):
    mini=info_group.to_crs(epsg=4326)
    label="\n".join(mini[["Nombre", "location"]].apply(lambda x: "(%s)  %s"%(x.location, x.Nombre), axis=1))
    return label


####################################

cmap_squares="Blues"
cmap_regiones="Oranges"
cmap_comunas="RdPu"
cmap_elevation="YlGn"
cmap_stations="copper"

####################################
## Shapefiles
title_commune=""
title_square="" 
title_region=""

## Shapefile to use in case zoom is too big
try:
    regiones=shps.read_main("Regiones", 3857)
    comunas=shps.read_main("Comunas", 3857)
except:
    regiones=pd.DataFrame()
    comunas=pd.DataFrame()

if (2 in IncF.i_ModelFilters) or (2 in IncF.i_StationFilters):
    title_square="\n$\\bf{Square \, Filter}$"
    squares=[Polygon([[pol[2], pol[0]], [pol[3], pol[0]], [pol[3], pol[1]], [pol[2], pol[1]]]) for pol in IncF.i_SqrRegions]
    squares = gpd.GeoDataFrame({"geometry":squares, "Alias":IncF.c_SqrRegionsAlias}, geometry='geometry', crs="epsg:4326").to_crs(epsg=3857)
    cmap_squares = cm.get_cmap(cmap_squares)
    cmap_squares = truncate_colormap(cmap_squares, 0.2, 0.95)

if (4 in IncF.i_ModelFilters) or (4 in IncF.i_StationFilters):
    title_region = "\n$\\bf{Chilean \, Region \, Filter}$"
    regiones=shps.read_main("Regiones", 3857)
    cmap_regiones = cm.get_cmap(cmap_regiones)
    cmap_regiones = truncate_colormap(cmap_regiones, 0.2, 0.95)

if (5 in IncF.i_ModelFilters) or (5 in IncF.i_StationFilters):
    title_commune = "\n$\\bf{Chilean \, Commune \, Filter}$"
    comunas=shps.read_main("Comunas", 3857)
    cmap_comunas = cm.get_cmap(cmap_comunas)
    cmap_comunas = truncate_colormap(cmap_comunas, 0.2, 0.95)

if 6 in IncF.i_StationFilters:
    cmap_elevation = cm.get_cmap(cmap_elevation)
    cmap_elevation = truncate_colormap(cmap_elevation, 0.2, 0.95)  

if 3 in IncF.i_StationFilters:
    cmap_stations = cm.get_cmap(cmap_stations)
    cmap_stations = truncate_colormap(cmap_stations, 0.2, 0.95)    
titles=[title_commune, title_square, title_region]

ms=100
lwe=0.8
lw=3
factor=0.8
markers=["o", "d", ">", "<", "s", "*"]*3

####################################
def model_filters_image():
    #########
    ### Plots
    #fig, ax = plt.subplots(figsize=(10, 10))
    fig = plt.figure(figsize=(12, 10))
    ax=fig.add_axes([0.02, 0.02, 0.55, 0.96])
    ax_legend=fig.add_axes([0.6, 0.02, 0.36, 0.96]) 
    ax.set_axis_off()
    ax_legend.set_axis_off()
    all_labels=[]
    if 2 in IncF.i_ModelFilters:
        ax.add_line(Line2D([], [], color="none", label=title_square, zorder=200)) 
        all_labels.append(title_square)
        ngroups=squares.shape[0]
        for count, row in squares.iterrows():
            row=row.to_frame().T 
            row=gpd.GeoDataFrame(row, geometry=row.geometry)
            label=row['Alias'].values[0]
            all_labels.append(label)
            row.boundary.plot(ax=ax, color='k', linewidth=lw)
            row.boundary.plot(color=cmap_squares((0.5/ngroups + count)/ngroups), ax=ax, zorder=200+count, label=label, linewidth=factor*lw)
    if 4 in IncF.i_ModelFilters:
        ax.add_line(Line2D([], [], color="none", label=title_region, zorder=100))
        all_labels.append(title_region)
        ngroups=len(IncF.c_Regions)
        for count, group in enumerate(IncF.c_Regions):
            label=IncF.c_RegionsAlias[count]
            all_labels.append(label)
            regiones[regiones["codregion"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
            regiones[regiones["codregion"].isin(group)].boundary.plot(color=cmap_regiones((0.5/ngroups + count)/ngroups), ax=ax, zorder=100+count, label=label, linewidth=factor*lw)
    if 5 in IncF.i_ModelFilters:
        ax.add_line(Line2D([], [], color="none", label=title_commune, zorder=150))
        all_labels.append(title_commune)
        ngroups=len(IncF.c_Comunas)
        for count, group in enumerate(IncF.c_Comunas):
            label=IncF.c_ComunasAlias[count]
            all_labels.append(label)
            comunas[comunas["cod_comuna"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
            comunas[comunas["cod_comuna"].isin(group)].boundary.plot(color=cmap_comunas((0.5/ngroups + count)/ngroups), ax=ax, zorder=150+count, label=label, linewidth=factor*lw)
    ########### https://contextily.readthedocs.io/en/latest/providers_deepdive.html
    ## Terrain
    #'Stamen.TerrainBackground'
    ## Similar to google maps
    #'CartoDB.Voyager'
    #ctx.providers.CartoDB.Voyager
    

    ## Add background
    try:
        ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5)
    except:
        ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5, zoom=18)

    ## Legend  
    handles, labels = reorderLegend(ax=ax, order=all_labels)
    leg = ax_legend.legend(handles=handles, labels=labels, fontsize=12, loc='upper left', bbox_to_anchor=(0.05, 1), ncol=1, fancybox=True, framealpha=1, frameon=False)



    for item, label in zip(leg.legendHandles, leg.texts):
        if label._text  in [title_square, title_region, title_commune]:
            width=item.get_window_extent(fig.canvas.get_renderer()).width
            label.set_ha('left')
            label.set_position((-2*width,0))

    name="Model_filters.pdf"
    print("   ...Saving Filter Map of Models to: ", name, flush=True)   
    fig.savefig(join(IncF.c_FigureDir, name), format="pdf")

    extent = ax_legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(join(IncF.c_FigureDir, 'Info_Models.pdf'), bbox_inches=extent)
    plt.close(fig)
    plt.clf()


####################################
def obs_filters_image(t_Stations_Filters, t_Stations_Info):
    #fig, ax = plt.subplots(figsize=(12, 10))
    fig = plt.figure(figsize=(12, 10))
    ax=fig.add_axes([0.02, 0.02, 0.55, 0.96])
    ax_legend=fig.add_axes([0.6, 0.02, 0.36, 0.96]) 
    ax.set_axis_off()
    ax_legend.set_axis_off()

    networks=t_Stations_Info.index.get_level_values(0).drop_duplicates() 
    all_labels=[]
    labels_stations=[]
    
    for index_network, net in enumerate(networks):

        net_info=aux.get_info_group(t_Stations_Info, net)
        filters=t_Stations_Filters[net]
        net_vars=IncF.c_ObsVars[index_network]

        title_network="\n$\\bf{Network \, %s:  %s}$"%(net.replace("NORTE_RM", "NORTE RM"), " , ".join(net_vars))
        titles.append(title_network)

        ax.add_line(Line2D([], [], color="none", label=title_network)) 
        labels_stations.append(title_network)


        #########################
        if 2 in IncF.i_StationFilters:
            stations=filters["f2"]
            if index_network==0:
                ax.add_line(Line2D([], [], color="none", label=title_square, zorder=200)) 
                all_labels.append(title_square)
            ngroups=squares.shape[0]
            for count, row in squares.iterrows():
                if index_network==0:
                    ######### Squares
                    row=row.to_frame().T 
                    row=gpd.GeoDataFrame(row, geometry=row.geometry)
                    label=row['Alias'].values[0]
                    all_labels.append(label)
                    row.boundary.plot(ax=ax, color='k', linewidth=lw)
                    row.boundary.plot(color=cmap_squares((0.5/ngroups + count)/ngroups), ax=ax, zorder=200+count, label=label, linewidth=factor*lw)

                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_squares((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    labels_stations.append(estaciones)


        if 4 in IncF.i_StationFilters:
            stations=filters["f4"]
            if index_network==0:
                ax.add_line(Line2D([], [], color="none", label=title_region, zorder=100))
                all_labels.append(title_region)
            ngroups=len(IncF.c_Regions)
            for count, group in enumerate(IncF.c_Regions):
                if index_network==0:
                    label=IncF.c_RegionsAlias[count]
                    all_labels.append(label)
                    regiones[regiones["codregion"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
                    regiones[regiones["codregion"].isin(group)].boundary.plot(color=cmap_regiones((0.5/ngroups + count)/ngroups), ax=ax, zorder=100+count, label=label, linewidth = factor*lw)

                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_regiones((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    labels_stations.append(estaciones)


        if 5 in IncF.i_StationFilters:
            stations=filters["f5"]
            if index_network==0:
                ax.add_line(Line2D([], [], color="none", label=title_commune, zorder=150))
                all_labels.append(title_commune)
            ngroups=len(IncF.c_Comunas)
            for count, group in enumerate(IncF.c_Comunas):
                if index_network==0:
                    label=IncF.c_ComunasAlias[count]
                    all_labels.append(label)
                    comunas[comunas["cod_comuna"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
                    comunas[comunas["cod_comuna"].isin(group)].boundary.plot(color=cmap_comunas((0.5/ngroups + count)/ngroups), ax=ax, zorder=150+count, label=label, linewidth=factor*lw)

                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_comunas((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    labels_stations.append(estaciones)

        if 6 in IncF.i_StationFilters:
            stations=filters["f6"]
            ngroups=len(IncF.i_Elevation)
            for count, group in enumerate(IncF.i_Elevation):
                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    title_ele="$\\bf{Elevation: \, %s , %s}$"%(group[0], group[1])
                    ax.add_line(Line2D([], [], color="none", label=title_ele, zorder=150))
                    labels_stations.append(title_ele)
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_elevation((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    labels_stations.append(estaciones)
        if 3 in IncF.i_StationFilters:
            stations=filters["f3"]
            ngroups=len(IncF.c_Stations_Obs[index_network])
            for count, group in enumerate(IncF.c_Stations_Obs[index_network]):
                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    title_ele="$\\bf{Stations \, ID:  \, %s}$"%(IncF.c_Stations_Obs_Aliases[index_network][count])
                    ax.add_line(Line2D([], [], color="none", label=title_ele, zorder=150))
                    labels_stations.append(title_ele)
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_stations((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    labels_stations.append(estaciones)
        ## No filter
        if len(IncF.i_StationFilters)==0:
            estaciones=get_label_stations(net_info)
            title_ele="$\\bf{No filters}$"
            ax.add_line(Line2D([], [], color="none", label=title_ele, zorder=150))
            labels_stations.append(title_ele)
            net_info.plot(marker=markers[index_network], ax=ax, label=estaciones, color="orange", edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
            labels_stations.append(estaciones) 

    ## Add background
    try:
        ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5)
    except:
        try:
            ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5, zoom=18)
        except:
            gpd.sjoin(regiones, t_Stations_Info, how="inner").dropna().boundary.plot(ax=ax, edgecolors='black', linewidth=lwe)
            try:
                ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5)
            except:
                pass



    leg_or=ax.legend(fontsize=12, loc='upper left', bbox_to_anchor=(1.05, 1), ncol=1, fancybox=True, framealpha=1, frameon=False)
    renderer=fig.canvas.get_renderer()
    width=leg_or.get_window_extent(renderer).width-leg_or.legendHandles[1].get_window_extent(renderer).width 
    characters=int(width/12)+30
    #fig.savefig(join(IncF.c_FigureDir, "Stations_filters_inbetween.pdf"), format="pdf")

    """
    ## Add labels to ppints
    for x, y, label in zip(t_Stations_Info.geometry.x, t_Stations_Info.geometry.y, t_Stations_Info.Nombre):
        ax.annotate(label, xy=(x, y), xytext=(4, 4), textcoords="offset points")
    #"""
    #########
    ## Add names of stations that have both measurements
    ## Plotting limits
    left, right = plt.xlim()
    texts=[]
    for x, y, name in zip(t_Stations_Info.geometry.x, t_Stations_Info.geometry.y, t_Stations_Info.Nombre):

        texts.append(ax.text(x, y, name, fontsize=10,  verticalalignment='center'))
        """
        ## Saving Data Station name
        if x<(left+right)/2:
            #texts.append(ax.annotate(name, xy=(x, y), xytext=(4, 4)))
            texts.append(ax.text(x, y, name, fontsize=10, horizontalalignment='left', verticalalignment='center'))
        ## To the right if points are near the right edge
        else:
            texts.append(ax.text(x, y, name, fontsize=10, horizontalalignment='right', verticalalignment='center'))
            #texts.append(ax.annotate(name, xy=(x, y), xytext=(4, 4)))
        """
    ## Adjusts text placement
    #adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), expand_points=(1.01, 1.05))#, only_move={'points':'y', 'text':'y'})
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), expand_points=(1.1, 1.1), expand_text=(1.2, 1.2), ax=ax, force_text=(0.5, 0.5)) 
    #fig.savefig("trial.pdf", format="pdf")  
    #breakpoint()
        
    """
    if BoolPlotNames and BoolStationsExist:
        texts=[]
        for index, row in df2.iterrows():
            name="  %s  "%row["Nombre"]
            x=row["X"]
            y=row["Y"]
            idweb=row["ID-Web"]

            ## Saving Data Station name
            if x<(left+right)/2:
                texts.append(plt.text(x, y, name, fontsize=10, horizontalalignment='left',verticalalignment='center'))
            ## To the right if points are near the right edge
            else:
                texts.append(plt.text(x, y, name, fontsize=10, horizontalalignment='right',verticalalignment='center'))
        ## Adjusts text placement
        adjust_text(texts, only_move={'points':'y', 'text':'y'})
        """
    #fig.savefig(join(IncF.c_FigureDir, "Stations_filters_inbetween2.pdf"), format="pdf")
    ## Legend  
    handles, labels = reorderLegend(ax=ax, order=all_labels+labels_stations, characters=characters)

    ax.get_legend().remove()


    leg = ax_legend.legend(handles=handles, labels=labels, fontsize=12, loc='upper left', bbox_to_anchor=(0.05, 1.05), ncol=1, fancybox=True, framealpha=1, frameon=False)


    for item, label in zip(leg.legendHandles, leg.texts):
        if label._text  in titles:
            width=item.get_window_extent(renderer).width
            label.set_ha('left')
            label.set_position( (-2*width, 0))

    name="Stations_filters.pdf"
    print("   ...Saving Filter Map of Networks to: ", name, flush=True) 

    fig.savefig(join(IncF.c_FigureDir, name), format="pdf")

    extent = ax_legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(join(IncF.c_FigureDir, 'Info_Stations.pdf'), bbox_inches=extent)

    plt.close(fig)
    plt.clf()


####################################
def both_filters_image(t_Stations_Filters, t_Stations_Info):
    #########
    ### Plots
    #fig, ax = plt.subplots(figsize=(10, 10))
    networks=t_Stations_Info.index.get_level_values(0).drop_duplicates()
    fig = plt.figure(figsize=(12, 10))
    ax=fig.add_axes([0.02, 0.02, 0.55, 0.96])
    ax_legend=fig.add_axes([0.6, 0.02, 0.36, 0.96]) 
    ax.set_axis_off()
    ax_legend.set_axis_off()
    all_labels=[]
    if 2 in IncF.i_ModelFilters + IncF.i_StationFilters:
        ax.add_line(Line2D([], [], color="none", label=title_square, zorder=200)) 
        all_labels.append(title_square)
        ngroups=squares.shape[0]
        for count, row in squares.iterrows():
            row=row.to_frame().T 
            row=gpd.GeoDataFrame(row, geometry=row.geometry)
            label=row['Alias'].values[0]
            all_labels.append(label)
            row.boundary.plot(ax=ax, color='k', linewidth=lw)
            row.boundary.plot(color=cmap_squares((0.5/ngroups + count)/ngroups), ax=ax, zorder=200+count, label=label, linewidth=factor*lw)
    if 4 in IncF.i_ModelFilters+ IncF.i_StationFilters:
        ax.add_line(Line2D([], [], color="none", label=title_region, zorder=100))
        all_labels.append(title_region)
        ngroups=len(IncF.c_Regions)
        for count, group in enumerate(IncF.c_Regions):
            label=IncF.c_RegionsAlias[count]
            all_labels.append(label)
            regiones[regiones["codregion"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
            regiones[regiones["codregion"].isin(group)].boundary.plot(color=cmap_regiones((0.5/ngroups + count)/ngroups), ax=ax, zorder=100+count, label=label, linewidth=factor*lw)
    if 5 in IncF.i_ModelFilters+ IncF.i_StationFilters:
        ax.add_line(Line2D([], [], color="none", label=title_commune, zorder=150))
        all_labels.append(title_commune)
        ngroups=len(IncF.c_Comunas)
        for count, group in enumerate(IncF.c_Comunas):
            label=IncF.c_ComunasAlias[count]
            all_labels.append(label)
            comunas[comunas["cod_comuna"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
            comunas[comunas["cod_comuna"].isin(group)].boundary.plot(color=cmap_comunas((0.5/ngroups + count)/ngroups), ax=ax, zorder=150+count, label=label, linewidth=factor*lw)
    ########### https://contextily.readthedocs.io/en/latest/providers_deepdive.html
    ## Terrain
    #'Stamen.TerrainBackground'
    ## Similar to google maps
    #'CartoDB.Voyager'
    #ctx.providers.CartoDB.Voyager
    

    
    for index_network, net in enumerate(networks):

        net_info=aux.get_info_group(t_Stations_Info, net)
        filters=t_Stations_Filters[net]
        net_vars=IncF.c_ObsVars[index_network]

        title_network="\n$\\bf{Network \, %s:  %s}$"%(net.replace("NORTE_RM", "NORTE RM"), " , ".join(net_vars))
        titles.append(title_network)

        ax.add_line(Line2D([], [], color="none", label=title_network)) 
        all_labels.append(title_network)


        ##########################
        if 2 in IncF.i_StationFilters:
            stations=filters["f2"]
            if index_network==0:
                ax.add_line(Line2D([], [], color="none", label=title_square, zorder=200)) 
                all_labels.append(title_square)
            ngroups=squares.shape[0]
            for count, row in squares.iterrows():
                if index_network==0:
                    ######### Squares
                    row=row.to_frame().T 
                    row=gpd.GeoDataFrame(row, geometry=row.geometry)
                    label=row['Alias'].values[0]
                    all_labels.append(label)
                    row.boundary.plot(ax=ax, color='k', linewidth=lw)
                    row.boundary.plot(color=cmap_squares((0.5/ngroups + count)/ngroups), ax=ax, zorder=200+count, label=label, linewidth=factor*lw)

                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_squares((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    all_labels.append(estaciones)


        if 4 in IncF.i_StationFilters:
            stations=filters["f4"]
            if index_network==0:
                ax.add_line(Line2D([], [], color="none", label=title_region, zorder=100))
                all_labels.append(title_region)
            ngroups=len(IncF.c_Regions)
            for count, group in enumerate(IncF.c_Regions):
                if index_network==0:
                    label=IncF.c_RegionsAlias[count]
                    all_labels.append(label)
                    regiones[regiones["codregion"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
                    regiones[regiones["codregion"].isin(group)].boundary.plot(color=cmap_regiones((0.5/ngroups + count)/ngroups), ax=ax, zorder=100+count, label=label, linewidth = factor*lw)

                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_regiones((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    all_labels.append(estaciones)


        if 5 in IncF.i_StationFilters:
            stations=filters["f5"]
            if index_network==0:
                ax.add_line(Line2D([], [], color="none", label=title_commune, zorder=150))
                all_labels.append(title_commune)
            ngroups=len(IncF.c_Comunas)
            for count, group in enumerate(IncF.c_Comunas):
                if index_network==0:
                    label=IncF.c_ComunasAlias[count]
                    all_labels.append(label)
                    comunas[comunas["cod_comuna"].isin(group)].boundary.plot(ax=ax, color='k', linewidth=lw)
                    comunas[comunas["cod_comuna"].isin(group)].boundary.plot(color=cmap_comunas((0.5/ngroups + count)/ngroups), ax=ax, zorder=150+count, label=label, linewidth=factor*lw)

                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_comunas((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    all_labels.append(estaciones)

        if 6 in IncF.i_StationFilters:
            stations=filters["f6"]
            ngroups=len(IncF.i_Elevation)
            for count, group in enumerate(IncF.i_Elevation):
                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    title_ele="$\\bf{Elevation: \, %s , %s}$"%(group[0], group[1])
                    ax.add_line(Line2D([], [], color="none", label=title_ele, zorder=150))
                    all_labels.append(title_ele)
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_elevation((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    all_labels.append(estaciones)
        if 3 in IncF.i_StationFilters:
            stations=filters["f3"]
            ngroups=len(IncF.c_Stations_Obs[index_network])
            for count, group in enumerate(IncF.c_Stations_Obs[index_network]):
                ######## Estaciones
                info_group=net_info[net_info["ID"].isin(stations[count])]
                if not info_group.empty:
                    title_ele="$\\bf{Stations \, ID:  \, %s}$"%(IncF.c_Stations_Obs_Aliases[index_network][count])
                    ax.add_line(Line2D([], [], color="none", label=title_ele, zorder=150))
                    all_labels.append(title_ele)
                    estaciones=get_label_stations(info_group)
                    info_group.plot(marker=markers[index_network], ax=ax, label=estaciones, color=cmap_stations((0.5/ngroups + count)/ngroups), edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
                    all_labels.append(estaciones)
        ## No filter
        if len(IncF.i_StationFilters)==0:
            estaciones=get_label_stations(net_info)
            title_ele="$\\bf{No filters}$"
            ax.add_line(Line2D([], [], color="none", label=title_ele, zorder=150))
            all_labels.append(title_ele)
            net_info.plot(marker=markers[index_network], ax=ax, label=estaciones, color="orange", edgecolors='black', markersize=ms, linewidth=lwe, zorder=500)
            all_labels.append(estaciones) 

    ## Add background
    try:
        ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5)
    except:
        try:
            ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5, zoom=18)
        except:
            gpd.sjoin(chile, t_Stations_Info, how="inner").dropna().boundary.plot(ax=ax, edgecolors='black', linewidth=lwe)
            try:
                ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, alpha=0.5)
            except:
                pass


    ## Legend  
    handles, labels = reorderLegend(ax=ax, order=all_labels)
    leg = ax_legend.legend(handles=handles, labels=labels, fontsize=12, loc='upper left', bbox_to_anchor=(0.05, 1), ncol=1, fancybox=True, framealpha=1, frameon=False)

    ## Add labels to ppints
    for x, y, label in zip(t_Stations_Info.geometry.x, t_Stations_Info.geometry.y, t_Stations_Info.Nombre):
        ax.annotate(label, xy=(x, y), xytext=(5, 5), textcoords="offset points", zorder=1000)


    for item, label in zip(leg.legendHandles, leg.texts):
        if label._text  in [title_square, title_region, title_commune]:
            width=item.get_window_extent(fig.canvas.get_renderer()).width
            label.set_ha('left')
            label.set_position((-2*width,0))

    name="Both_filters.pdf"
    print("   ...Saving Filter Map of both to: ", name, flush=True)   
    fig.savefig(join(IncF.c_FigureDir, name), format="pdf")

    extent = ax_legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(join(IncF.c_FigureDir, 'Info_Models.pdf'), bbox_inches=extent)
    plt.close(fig)
    plt.clf()


###############################################################
###############################################################
#####################################
def main(t_Stations_Info, t_Stations_Filters):
    if BoolInitialMap:
        if not t_Stations_Info.empty and len(IncF.i_ModelFilters)!=0:
            both_filters_image(t_Stations_Filters, t_Stations_Info.to_crs(epsg=3857))

        ###### Networks filters
        if not t_Stations_Info.empty:
            obs_filters_image(t_Stations_Filters, t_Stations_Info.to_crs(epsg=3857))
       

        ###### Model Filters
        if len(IncF.i_ModelFilters)!=0:
            model_filters_image()
        else:
            print("...Skipping this step. No filters were applied for models")
    else:
        print("...Skipping this step.")
    
