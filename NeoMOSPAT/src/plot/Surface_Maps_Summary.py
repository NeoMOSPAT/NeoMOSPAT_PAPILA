# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
import src.common as aux
import src.plot.aux_plots as aux_plots
import src.plot.aux_maps as aux_maps
#
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.colors import BoundaryNorm, CenteredNorm
from matplotlib import  cm
#
import numpy as np
import pandas as pd
import geopandas as gpd
#
from os.path import join
from os import makedirs
from datetime import datetime
from copy import copy
from matplotlib import ticker

#plt.rcParams['fig.alpha'] = 0
plt.rcParams['axes.facecolor'] = 'white'
#plt.rcParams['axes.alpha'] = 1
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.edgecolor'] = 'none'


list_maps=[int(c) for c in str(IncF.i_Map_Sum)]
list_addnet=[int(c) for c in str(IncF.c_Map_Sum_AddNet)]

######################################
settings=aux_plots.BaseMapConfig
#cmap=settings.cmap
figarea=settings.image_area
titlesize=settings.titlesize
ticksize=settings.ticksize
labelsize=settings.labelsize
legendsize=settings.legendsize
ms=settings.ms
nlevels=settings.nlevels
markerlist=settings.markerlist
format_fig = settings.format_fig
BoolCropContent=settings.BoolCropContent
pad_inches=settings.pad_inches
center_cmap=True


#############
## Creates name for plot images
def return_name(var, model, vmin, vmax, time_srtf, typle_plot="", add_net="", val="", plot_type=""):
    return "%s-%s__H%s-%s__vmin%.2f-vmax%.2f__%s__%s"%(var, model, val.hour, plot_type, vmin, vmax, time_srtf, "+".join([s for s in [typle_plot, add_net] if s !=""]))



#############
## Plot all types of maps
def plot_map(x, y, z, var, unit, model, time, bounds, obs_data, t_Stations_Info, vmin, vmax, direc="", val="", plot_type="", cmap="", map_type=""):
    ######
    fig = plt.figure(figsize=aux_maps.get_figsize(bounds))
    ax = fig.add_subplot(1, 1, 1)
    ## time
    time_strf = time if isinstance(time, str) else time.strftime("%Y-%m-%d %H:%M:%S").replace("00:00:00", "")

    ## Add background and Shapefile
    ax, cmap = aux_maps.prepare_maps(cmap, bounds, ax=ax)
    ## Color Bar
    extend=aux_maps.find_extend(vmin, vmax, np.nanmin(z), np.nanmax(z))
    

    ## Pixelated plot
    if map_type=="Square":
        cs=ax.pcolor(x, y, z[:-1,:-1], shading="flat", cmap=cmap, vmin=vmin, vmax=vmax, zorder=1, snap=True, edgecolor=None)
        type_plot="GRID"
        ## Format Map (ticks, colorbar, etc)
        aux_maps.format_map(fig, ax, model, var, unit, time_strf, extend, cs, vmin, vmax)

    ## Interpolated plot
    elif map_type=="Interpolated":
        cs=ax.pcolormesh(x, y, z, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, zorder=1, snap=True, edgecolor=None)
        type_plot="SMOOTH"

        ## Format Map (ticks, colorbar, etc)
        aux_maps.format_map(fig, ax, model, var, unit, time_strf, extend, cs, vmin, vmax)

    ## Contour fill plot
    elif map_type=="Contour":
        cs=ax.contourf(x, y, z, cmap=cmap, levels=nlevels, vmin=vmin, vmax=vmax, zorder=1)  
        type_plot="CTRF"
        ## Format Map (ticks, colorbar, etc)
        aux_maps.format_map(fig, ax, model, var, unit, time_strf, extend, cs, vmin, vmax, type_plot=type_plot)




    


    ## Regular without any data station data
    if 0 in list_addnet:
        name=return_name(var, model, vmin, vmax, time_strf, typle_plot=type_plot, add_net="", val=val, plot_type=plot_type)
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)

    ## Add locations only
    if 1 in list_addnet:
        ax=aux_maps.add_stations(obs_data, t_Stations_Info, var, time, vmin, vmax, cmap, ax=ax, case=1)
        name=return_name(var, model, vmin, vmax, time_strf, typle_plot=type_plot, add_net="stations", val=val, plot_type=plot_type)
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)

    fig.clf()
    plt.close(fig)




###############################################################
###############################################################
def direct_to_plot_maps(data_resampled, vars_to_plot, units, model, bounds, model_info, obs_data, t_Stations_Info, plot_freq, colorbar, direc="", filter="", plot_type=""):

    ## Square information in case the data being plotted is not a square
    first_date=data_resampled.index.values[0]
    data_first_date=data_resampled.loc[first_date]
    ilats=data_first_date["i_lat"].drop_duplicates(keep="first")
    ilons=data_first_date["i_lon"].drop_duplicates(keep="first")
    lat_min=ilats.min()
    lon_min=ilons.min()
    lat_max=ilats.max()
    lon_max=ilons.max()
    len_x = lon_max+1 - lon_min
    len_y = lat_max+1 - lat_min
 
    indexes_should_have= [(lat, lon) for lat in range(lat_min, lat_max+1) for lon in range(lon_min, lon_max+1)]
    grid_should_have=model_info.set_index("indexes").reindex(indexes_should_have).reset_index()



    ##############
    if plot_freq=="Whole":
        data_resampled["groupby"]="%s %s"%(data_resampled.index[0].strftime("%Y-%m-%d"), data_resampled.index[-1].strftime("%Y-%m-%d"))
        grouped_data=data_resampled.groupby("groupby")
    else:
        grouped_data=data_resampled.groupby(pd.Grouper(freq=plot_freq))
          
    for time in list(grouped_data.groups.keys()):
        ## Get dataframe
        try:
            group_data=grouped_data.get_group(time)        
            time0=str(group_data.index.values[0]).split("T")[0]
            time1=str(group_data.index.values[-1]).split("T")[0]
            time0=datetime.strptime(str(group_data.index.values[0]), "%Y-%m-%dT%H:%M:%S0.000000000").strftime("%Y-%m-%d %H:%M:%S")
            time1=datetime.strptime(str(group_data.index.values[-1]), "%Y-%m-%dT%H:%M:%S0.000000000").strftime("%Y-%m-%d %H:%M:%S")

            ## Get info of cycles
            if plot_type=="DC":
                group_data["groupby"] = group_data.index.map(lambda t: t.replace(year=2000, month=1, day=1)) 
            elif plot_type=="AC":
                group_data["groupby"] = group_data.index.map(lambda t: t.replace(year=2000)) 

            group_data=group_data.groupby(["groupby", "i_lat", "i_lon"]).mean().reset_index(level=["i_lat", "i_lon"])
        except:
            continue

        data_plot=group_data[vars_to_plot+["lat", "lon", "i_lat", "i_lon"]]

        time="%s - %s"%(time0, time1)

        ## Values to plot
        for val in data_plot.index.unique():
            data=data_plot.loc[val]
            ## Make plots square
            data=data.merge(grid_should_have, on=["lat", "lon", "i_lat", "i_lon"], how="outer")[data.columns]
            data.sort_values(by=["i_lon", "i_lat"], inplace=True)

            x=data["lon"].values.reshape((len_x, len_y)) 
            y=data["lat"].values.reshape((len_x, len_y)) 
            if x.shape[0]<2 or x.shape[1]<2:
                print("      ...Area is too small to plot. Less than two grids in at least one axis")
                continue

            ## Obs data of period
            obs_data_period={k:obs_data[k].mean() for k in obs_data}  
            ## Plot vars requested
            for var in vars_to_plot:
                z=data[var].values.reshape((len_x, len_y))
                z_global=data_plot[var]
                vmin, vmax = aux_maps.find_range_cb(z, z_global, colorbar)
                ## Si variable no tiene valores
                if not np.nanmin(z)!=np.nan:
                    return None

                cmap=aux_maps.get_cmap(vmin, vmax, var)
                if 1 in list_maps:
                    plot_map(x, y, z, var, units[var], model, time, bounds, obs_data_period, t_Stations_Info, vmin, vmax, direc=direc, val=val, plot_type=plot_type, cmap=cmap, map_type="Square")
                if 2 in list_maps:
                    plot_map(x, y, z, var, units[var], model, time, bounds, obs_data_period, t_Stations_Info, vmin, vmax, direc=direc, val=val, plot_type=plot_type, cmap=cmap, map_type="Interpolated")
                if 3 in list_maps:
                    plot_map(x, y, z, var, units[var], model, time, bounds, obs_data_period, t_Stations_Info, vmin, vmax, direc=direc, val=val, plot_type=plot_type, cmap=cmap, map_type="Contour")






###############################################################
###############################################################
def main(t_ObsStationData, t_ModelData, t_Stations_Info, t_Model_Filters, t_Models_Info):
    if IncF.i_Map_Sum==0:
        print("...Plotting was not requested.")
        return None
    if len([c for c in IncF.c_ModelNames if c])==0:
        print("...Plotting was requested but no model was read.")
        return None  

    ## Find vars that need to be plotted
    vars_to_plot=[aux.find_equiv_vars([v]) for v in IncF.c_Map_Sum_Plots]
    if len(vars_to_plot)==0:
        print("...Plotting was requested but variables are not available.")
        return None    


    ########
    path=join(IncF.c_FigureDir, "SurfaceMaps")
    makedirs(path, exist_ok = True)       
    #for freq in IncF.c_Map_Sum_freq:
    for count, plot_freq in enumerate(IncF.c_Map_Sum_freq):
        freq, plot_type, plot_freq = plot_freq.split(",")
        ## Goes through network asked and finds the filters associated to it
        for index_cb, plot in enumerate(vars_to_plot):
            for model in plot["M"]:
                print(f"...Creating plots for  '{IncF.c_Map_Sum_Plots[index_cb]}'  with data frequency  '{freq}'  and plot frequency  '{plot_freq}'")
                variables=vars_to_plot[index_cb]
                colorbar=IncF.c_Map_Sum_ColBar[index_cb]


                path_model=join(path, model)
                makedirs(path_model, exist_ok = True) 
                path_model=join(path_model, f"{freq},{plot_freq}")
                makedirs(path_model, exist_ok = True) 


                ###########
                ## Retrieve model Data needed
                units=t_ModelData[freq][model].meta.units
                data_resampled=t_ModelData[freq][model].reset_index(level=[1,2]).sort_values(by=["i_lat", "i_lon"]) 

                ## Vars that need to be done 
                vars_model=variables["M"][model]
                ## Add lat and lon data to data set
                model_info=t_Models_Info[model]
                data_resampled=data_resampled.merge(model_info, on=["i_lat", "i_lon"]).set_index(data_resampled.index)



                ##########
                ## filters 
                filters=t_Model_Filters[model]
                for f in filters:
                    for group, indexes_group in enumerate(filters[f]):
                            if len(indexes_group)==0:
                                continue
                            ######## Check why this happens
                            ## Without time as column
                            try:
                                data_plot_filter=data_resampled.reset_index().set_index(["i_lat", "i_lon"]).loc[indexes_group]
                            ## With time as column
                            except:
                                data_plot_filter=data_resampled.set_index(["i_lat", "i_lon"]).loc[indexes_group]

                            data_plot_filter=data_plot_filter.reset_index().set_index("time") 

                            ## Find square bounds for the plot
                            if f=="f2":
                                #[Lat_North, Lat_South, Lon_West, Lon_East]
                                bounds = IncF.i_SqrRegions[group]
                                # minx, miny, maxx, maxy 
                                bounds = [bounds[2], bounds[1], bounds[3], bounds[0]] 
                                bounds_models=model_info[model_info["indexes"].isin(indexes_group)].total_bounds
                                if BoolCropContent:
                                    if bounds[0]<bounds_models[0]:
                                        bounds[0]=bounds_models[0]
                                    if bounds[1]<bounds_models[1]:
                                        bounds[1]=bounds_models[1]
                                    if bounds[2]>bounds_models[2]:
                                        bounds[2]=bounds_models[2]
                                    if bounds[3]>bounds_models[3]:
                                        bounds[3]=bounds_models[3]
                            else:
                                bounds=model_info[model_info["indexes"].isin(indexes_group)].total_bounds

                            path_filter=join(path_model, "%s-%s_(%s)"%(f, group, ",".join(["%.2f"%b for b in bounds])))
                            makedirs(path_filter, exist_ok = True) 
                            print("   ...For filter area:  ", "%s-%s"%(f, group))
                            ## Do plots
                            direct_to_plot_maps(data_plot_filter, vars_model, units, model, bounds, model_info, t_ObsStationData[freq], t_Stations_Info, plot_freq, colorbar, direc=path_filter, filter="%s-%s"%(f, group), plot_type=plot_type)
                            plt.close('all')




