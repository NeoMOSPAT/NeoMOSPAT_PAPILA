# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
import src.common as aux
import src.plot.aux_plots as aux_plots
import src.plot.aux_maps as aux_maps
from src.write.statistics import make_stats
#
import matplotlib.pyplot as plt
from matplotlib import  cm
from matplotlib.colors import BoundaryNorm
#
import numpy as np
import pandas as pd
import geopandas as gpd
#
from os.path import join
from os import makedirs
from datetime import datetime
from copy import copy
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable


list_maps=[int(c) for c in str(IncF.i_Map_Stats)]
list_addnet=[int(c) for c in str(IncF.c_Map_Stats_AddNet)]
list_stats=[int(c) for c in str(IncF.c_Map_Stats_type)]
list_types=[int(c) for c in str(IncF.c_Map_Stats_bg)]

######################################
settings=aux_plots.BaseMapConfig
format_fig = settings.format_fig
BoolCropContent=settings.BoolCropContent
pad_inches=settings.pad_inches
nlevels=settings.nlevels
symmetric_centering=settings.symmetric_centering
cmap_stats="winter"

#############
## Creates name for plot images
def return_name(var, model, vmin, vmax, time_srtf, type_plot="", add_net=""):
    return "%s-%s__vmin%.2f-vmax%.2f__%s__%s"%(var, model, vmin, vmax, time_srtf, "+".join([s for s in [type_plot, add_net] if s !=""]))


#############
## All maps
def plot_map(x, y, z, var, unit, model, time, bounds, obs_data, t_Stations_Info, vmin, vmax, direc="", list_addnet=list_addnet, cmap="", map_type="", boolLog=False):
    ######
    fig = plt.figure(figsize=aux_maps.get_figsize(bounds))
    ax = fig.add_subplot(1, 1, 1)
    ## time
    time_strf = time if isinstance(time, str) else time.strftime("%Y-%m-%d %H:%M:%S").replace("00:00:00", "")

    ## Add background and Shapefile
    ax, alpha = aux_maps.prepare_maps(bounds, ax=ax)
    cmap, norm = aux_maps.prepare_colormaps(alpha, vmin, vmax, var, map_type, boolLog=boolLog)
    ## Color Bar 
    extend = aux_maps.find_extend(vmin, vmax, np.nanmin(z), np.nanmax(z))

    ## Pixelated plot
    if map_type=="Square":
        cs=ax.pcolor(x, y, z[:-1,:-1], shading="flat", cmap=cmap, zorder=1, snap=True, rasterized=True, antialiased=True, norm=norm)
        type_plot="GRID"

    ## Interpolated plot
    elif map_type=="Interpolated":
        cs=ax.pcolormesh(x, y, z, shading="gouraud", cmap=cmap, zorder=1, snap=True, linewidth=0, rasterized=True, antialiased=True, norm=norm)
        type_plot="SMOOTH"

    ## Contour fill plot
    elif map_type=="Contour":
        ## Makes it so that the level of the countour plot follows vmin and vmax instead of min and max of data
        if symmetric_centering:
            vm=max(abs(vmin), vmax)
            z[z>vm]=vm+0.001*(vmax-vmin)
            z[z<-vm]=-vm-0.001*(vmax-vmin) 
            extend = aux_maps.find_extend(-vm, vm, np.nanmin(z), np.nanmax(z))           
        else:
            z[z>vmax]=vmax+0.001*(vmax-vmin)
            z[z<vmin]=vmin-0.001*(vmax-vmin)
        cs=ax.contourf(x, y, z, cmap=cmap, zorder=1, levels=nlevels, norm=norm)  
        type_plot="CTRF"

    
    ## Format Map (ticks, colorbar, etc)
    aux_maps.format_map(fig, ax, model, var, unit, time_strf, extend, cs, vmin, vmax, type_plot=type_plot)
    info=t_Stations_Info[list(t_Stations_Info.keys())[0]]
    mini=info[~info["ID"].str.contains('^f[0-9]-+')]
    minif=info[info["ID"].str.contains('^f[0-9]-+')]
    #############################
    ## Regular without any data station data
    if 0 in list_addnet:
        name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="")
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)

    ## Add locations only
    if 1 in list_addnet:
        make_stats(data_models, data_obs, column_obs)
        ax=aux_maps.add_statistics(obs_data, mini, var, time, vmin, vmax, cmap, ax=ax, case=1)
        name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="stations")
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)

    ## Add locations only
    if 3 in list_addnet:
        make_stats(data_models, data_obs, column_obs)
        ax=aux_maps.add_statistics(obs_data, minif, var, time, vmin, vmax, cmap, ax=ax, case=1)
        name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="filters")
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)

    ## Add location with the value of the station
    if 2 in list_addnet:
        ax=aux_maps.add_statistics(obs_data, mini, var, time, vmin, vmax, cmap, ax=ax, case=2)
        name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="sta-vals")
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)


    ## Add location with the value of the station
    if 4 in list_addnet:
        mini_data=pd.concat([obs_data[n].loc[minif["ID"]] for n in obs_data])
        cax = fig.add_axes([ax.get_position().x0-0.165, ax.get_position().y0, 0.035, ax.get_position().height])
        data={n:obs_data[n].loc[minif["ID"]] for n in obs_data}
        for stat in mini_data.columns:
            vmin, vmax = aux_maps.add_statistics(data, minif, stat, time, cmap_stats, ax=ax, cax=cax, case=2, fig=fig)
            name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net=f"filters-{stat}")
            aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)
            cax.clear()
            breakpoint()
    fig.clf()
    plt.close(fig)


#############
def prepare_statistics(obs_data, obs_info, mod_data, model_info, vars_model):
    data={}
    for var in vars_model:
        data[var]={}
        data_network={}
        for network in obs_data:
            data_all=[]
            for ID in obs_info.loc[network]["ID"]:
                df=make_stats(obs_data[network][var][[ID]], mod_data[network][var][[ID]], ID).drop(columns=["Obs"]).T
                data_all.append(df)
            stats_network=pd.concat(data_all, keys=list(obs_info.loc[network]["ID"].values))
            data[var][network]=stats_network
        
    return data


###############################################################
###############################################################
def direct_to_plot_maps(data_resampled, vars_to_plot, units, model, bounds, model_info, statistics, info, plot_freq, colorbar, direc="", filter="", freq="H"):
    def create_agg_dic(vars_to_plot, freq):
        operations=[]
        ## Mean
        if 1 in list_types:
            operations.append("mean")
        ## Max
        if 2 in list_types:
            operations.append("max")
            if freq=="H":
                operations.append(("time-max", lambda x: x.idxmax().hour))
            elif freq=="D":
                operations.append(("time-max", lambda x: x.idxmax().day))
        ## Min
        if 3 in list_types:
            operations.append("min")
            if freq=="H":
                operations.append(("time-min", lambda x: x.idxmin().hour))
            elif freq=="D":
                operations.append(("time-min", lambda x: x.idxmin().day))


        ## Std
        if 4 in list_types:
            operations.append("std")
        
        return operations


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
    grid_should_have.columns = pd.MultiIndex.from_tuples([(c, "") for c in grid_should_have.columns]) 

    ## New
    #breakpoint()
    data_resampled=data_resampled[vars_to_plot+["lat", "lon", "i_lat", "i_lon"]]
    ops=create_agg_dic(vars_to_plot, freq)
    data_all=data_resampled.dropna().groupby(["i_lat", "i_lon", "lat", "lon"]).agg({v:ops for v in vars_to_plot}).reset_index()


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
        except:
            continue
       
        #""" 
        group_data=group_data[vars_to_plot+["lat", "lon", "i_lat", "i_lon"]]
        ops=create_agg_dic(vars_to_plot, freq)
        if len(ops)==0:
            data=group_data.copy()
            data[vars_to_plot]=np.nan
            data.reset_index(inplace=True)
        else:
            data=group_data.dropna().groupby(["i_lat", "i_lon", "lat", "lon"]).agg({v:ops for v in vars_to_plot}).reset_index()

        #"""
        
        time0=str(group_data.index.values[0]).split("T")[0]
        time1=str(group_data.index.values[-1]).split("T")[0]
        time="%s - %s"%(time0, time1)
        time0=datetime.strptime(str(group_data.index.values[0]), "%Y-%m-%dT%H:%M:%S0.000000000").strftime("%Y-%m-%d %H:%M:%S")
        time1=datetime.strptime(str(group_data.index.values[-1]), "%Y-%m-%dT%H:%M:%S0.000000000").strftime("%Y-%m-%d %H:%M:%S")
        

        ## Make plots square   
        try:
            data_plot=data.merge(grid_should_have, on=[("lat", ""), ("lon", ""), ("i_lat", ""), ("i_lon", "")], how="outer")[data.columns]
            data_plot.columns = pd.MultiIndex.from_tuples(data_plot.columns) 
        except:
            data_plot=data.merge(grid_should_have, right_on=[("lat", ""), ("lon", ""), ("i_lat", ""), ("i_lon", "")], left_on=["lat", "lon", "i_lat", "i_lon"], how="outer")[data.columns]

        data_plot.index=np.repeat(time, data_plot.shape[0]) 
        #statistics=prepare_statistics(data_models, data_obs, column_obs)

        ## Process data so it can be plotted
        data_plot.sort_values(by=["i_lon", "i_lat"], inplace=True)
        x=data_plot["lon"].values.reshape((len_x, len_y)) 
        y=data_plot["lat"].values.reshape((len_x, len_y)) 
        if x.shape[0]<2 or x.shape[1]<2:
            print("      ...Area is too small to plot. Less than two grids in at least one axis")
            continue
        ## Plot vars requested
        for var in vars_to_plot:
            for type_plot in data_plot[var].columns:
                z=data_plot[var][type_plot].values.reshape((len_x, len_y))
                z_global=data_all[var][type_plot]
                varcom="%s-%s"%(var, type_plot)
                vmin, vmax, boolLog = aux_maps.find_range_cb(z, z_global, colorbar)
                ## Si variable no tiene valores
                if not np.nanmin(z)!=np.nan:
                    return None

                ## If it's time only do the one without anything
                if type_plot=="time-max":
                    if freq=="H":
                        cmap=settings.cmap_cyclic
                        plot_map(x, y, z, varcom, "h", model, time, bounds, statistics, info, 0, 23, direc=direc, map_type="Square", boolLog=boolLog)
                    elif freq=="D":
                        cmap=settings.cmap_seq
                        plot_map(x, y, z, varcom, "d", model, time, bounds, statistics, info, 0, np.nanmax(z), direc=direc, map_type="Square", boolLog=boolLog)
                    continue

                    
                
                if 1 in list_maps:
                    plot_map(x, y, z, varcom, units[var], model, time, bounds, statistics[var], info, vmin, vmax, direc=direc, map_type="Square")
                if 2 in list_maps:
                    plot_map(x, y, z, varcom, units[var], model, time, bounds, statistics[var], info, vmin, vmax, direc=direc, map_type="Interpolated")
                if 3 in list_maps:
                    plot_map(x, y, z, varcom, units[var], model, time, bounds, statistics[var], info, vmin, vmax, direc=direc, map_type="Contour")
                


###############################################################
###############################################################
def main(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info):
    return None
    if IncF.i_Map_Stats==0:
        print("...Plotting was not requested.")
        return None
    if len([c for c in IncF.c_ModelNames if c])==0:
        print("...Plotting was requested but no model was read.")
        return None  

    ## Find vars that need to be plotted
    vars_to_plot=[aux.find_equiv_vars([v]) for v in IncF.c_Map_Stats_Plots]
    if len(vars_to_plot)==0:
        print("...Plotting was requested but variables are not available.")
        return None    


    ########
    path=join(IncF.c_FigureDir, "SurfaceMaps_Statistics")
    makedirs(path, exist_ok = True)       
    #for freq in IncF.c_Map_freq:
    for count, plot_freq in enumerate(IncF.c_Map_Stats_freq):
        freq, plot_freq =plot_freq.split(",")
        ## Goes through network asked and finds the filters associated to it
        for index_cb, plot in enumerate(vars_to_plot):
            for model in plot["M"]:
                print(f"...Creating plots for  '{IncF.c_Map_Stats_Plots[index_cb]}'  with data frequency  '{freq}'  and plot frequency  '{plot_freq}'")
                variables=vars_to_plot[index_cb]
                colorbar=IncF.c_Map_Stats_ColBar[index_cb]


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
                statistics=prepare_statistics(t_ObsStationData[freq], t_Stations_Info, t_ModelStationData[freq][model],  t_ModelStation_Info, vars_model)
                ## Add lat and lon data to data set
                model_info=t_Models_Info[model]
                data_resampled=data_resampled.merge(model_info, on=["i_lat", "i_lon"]).set_index(data_resampled.index)



                ##########
                ## filters 
                filters=t_Model_Filters[model]
                for f in filters:
                    for group, alias in enumerate(filters[f]):
                            indexes_group=filters[f][alias]["indexes"].values
                            if len(indexes_group)==0:
                                continue
                            indexes_group=pd.MultiIndex.from_tuples(indexes_group, names=("i_lat", "i_lon"))

                            ######## Check why this happens
                            ## Without time as column
                            try:
                                data_resampled = data_resampled.reset_index().set_index(["i_lat", "i_lon"])
                                indexes_group = indexes_group.intersection(data_resampled.index)
                                data_plot_filter = data_resampled.loc[indexes_group]
                            ## With time as column
                            except:
                                data_resampled = data_resampled.set_index(["i_lat", "i_lon"])
                                indexes_group = indexes_group.intersection(data_resampled.index)
                                data_plot_filter = data_resampled.loc[indexes_group]


                            data_plot_filter=data_plot_filter.reset_index().set_index("time") 

                            ## Find square bounds for the plot
                            if f=="f2":
                                #[Lat_North, Lat_South, Lon_West, Lon_East]
                                bounds = IncF.i_SqrRegions[group]
                                # minx, miny, maxx, maxy 
                                bounds = [bounds[2], bounds[1], bounds[3], bounds[0]] 
                            else:
                                bounds=model_info[model_info["indexes"].isin(indexes_group)].total_bounds

                            path_filter=join(path_model, "%s-%s"%(f, alias))
                            makedirs(path_filter, exist_ok = True) 
                            print("   ...For filter area:  ", "%s-%s"%(f, alias))
                            ## Do plots
                            direct_to_plot_maps(data_plot_filter, vars_model, units, model, bounds, model_info, statistics,  t_ModelStation_Info, plot_freq, colorbar, direc=path_filter, filter="%s-%s"%(f, group), freq=freq)
                            plt.close('all')


