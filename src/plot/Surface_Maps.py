# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
import src.common as aux
import src.plot.aux_plots as aux_plots
import src.plot.aux_maps as aux_maps


#
import matplotlib.pyplot as plt
from matplotlib import  cm
from matplotlib.colors import BoundaryNorm, Normalize
#
import numpy as np
import pandas as pd
import geopandas as gpd
#
from os.path import join
from os import makedirs
from datetime import datetime
from copy import copy



list_maps=[int(c) for c in str(IncF.i_Map)]
list_addnet=[int(c) for c in str(IncF.c_Map_AddNet)]
list_types=[int(c) for c in str(IncF.c_Map_type)]


######################################
settings=aux_plots.BaseMapConfig
format_fig = settings.format_fig
BoolCropContent=settings.BoolCropContent
pad_inches=settings.pad_inches
nlevels=settings.nlevels
symmetric_centering=settings.symmetric_centering


#############
## Creates name for plot images
def return_name(var, model, vmin, vmax, time_srtf, type_plot="", add_net=""):
    return "%s-%s__vmin%.2f-vmax%.2f__%s__%s"%(var, model, vmin, vmax, time_srtf, "+".join([s for s in [type_plot, add_net] if s !=""]))


#############
## All maps
def plot_map(x, y, z, var, unit, model, time, bounds, obs_data, t_Stations_Info, vmin, vmax, direc="", list_addnet=list_addnet, cmap="", map_type="", boolLog=False, statistics=None):
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
        cs=ax.pcolormesh(x, y, z, shading="nearest", cmap=cmap, zorder=1, snap=True,  linewidth=0, rasterized=True, antialiased=True, norm=norm)
        type_plot="GRID"

    ## Interpolated plot
    elif map_type=="Interpolated":
        cs=ax.pcolormesh(x, y, z, shading="gouraud", cmap=cmap, zorder=1, snap=True, linewidth=0, rasterized=True, antialiased=True, norm=norm)
        type_plot="SMOOTH"

    ## Contour fill plot
    elif map_type=="Contour":
        ## Makes it so that the level of the countour plot follows vmin and vmax instead of min and max of data
        if symmetric_centering:
            vm=max(abs(vmin), abs(vmax))
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


    #############################
    ## Regular without any data station data
    if 0 in list_addnet:
        name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="")
        aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)

    ## Add locations only
    if 1 in list_addnet:
        if not t_Stations_Info.empty:
            mini = t_Stations_Info[~t_Stations_Info["ID"].str.contains('^f[0-9]-+')]
            ax=aux_maps.add_stations(obs_data, mini, var, time, vmin, vmax, cmap, ax=ax, case=1)
            name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="stations")
            aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)

    ## Add locations only
    if 3 in list_addnet:
        if not t_Stations_Info.empty:
            mini = t_Stations_Info[t_Stations_Info["ID"].str.contains('^f[0-9]-+')]
            ax=aux_maps.add_stations(obs_data, mini, var, time, vmin, vmax, cmap, ax=ax, case=1)
            name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="filters")
            aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)

    ## Add location with the value of the station
    if 2 in list_addnet:
        if not t_Stations_Info.empty:
            mini = t_Stations_Info[~t_Stations_Info["ID"].str.contains('^f[0-9]-+')]
            ax=aux_maps.add_stations(obs_data, mini, var, time, vmin, vmax, cmap, ax=ax, case=2)
            name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="sta-vals")
            aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)

    ## Add location with the value of the filters
    if 4 in list_addnet:
        if not t_Stations_Info.empty:
            mini = t_Stations_Info[t_Stations_Info["ID"].str.contains('^f[0-9]-+')]
            ax=aux_maps.add_stations(obs_data, mini, var, time, vmin, vmax, cmap, ax=ax, case=2)
            name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net="filters-vals")
            aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)


    if IncF.i_Map_Stats!=0:
        direc=direc.replace("SurfaceMaps", "SurfaceMaps_Statistics")
        makedirs(direc, exist_ok = True) 
        ## Add location with the value of the station
        if 2 in list_addnet:
            mini = t_Stations_Info[~t_Stations_Info["ID"].str.contains('^f[0-9]-+')]
            mini_data=pd.concat([statistics[n].loc[mini["ID"]] for n in statistics])
            cax = fig.add_axes([ax.get_position().x0-0.165, ax.get_position().y0, 0.035, ax.get_position().height])
            data={n:statistics[n].loc[minif["ID"]] for n in statistics}
            for stat in mini_data.columns:
                vmin, vmax = aux_maps.add_statistics(data, minif, stat, time,  ax=ax, cax=cax, case=2, fig=fig)
                name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net=f"filters-{stat}")
                aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)
                cax.clear()
                #breakpoint()

        ## Add location with the value of the station
        if 4 in list_addnet:
            minif = t_Stations_Info[t_Stations_Info["ID"].str.contains('^f[0-9]-+')]
            mini_data=pd.concat([statistics[n].loc[minif.set_index("ID").index.intersection(statistics[n].index.get_level_values(0))] for n in statistics])
            cax = fig.add_axes([ax.get_position().x0-0.165, ax.get_position().y0, 0.035, ax.get_position().height])
            data={n:statistics[n].loc[minif.set_index("ID").index.intersection(statistics[n].index.get_level_values(0))] for n in statistics}
            for stat in mini_data.columns:
                vmin, vmax = aux_maps.add_statistics(data, minif, stat, time,  ax=ax, cax=cax, case=2, fig=fig)
                name=return_name(var, model, vmin, vmax, time_strf, type_plot=type_plot, add_net=f"filters-{stat}")
                aux_plots.savefigure(fig, direc, name, pad_inches=pad_inches, format=format_fig)#, transparent=True)
                cax.clear()
                #breakpoint()

    fig.clf()
    plt.close(fig)


###############################################################
###############################################################
def direct_to_plot_maps(data_resampled, vars_to_plot, units, model, bounds, model_info, obs_data, t_Stations_Info, plot_freq, colorbar, direc="", filter="", freq="H", statistics=None):
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
    #first_date=data_resampled
    data_first_date=data_resampled#.loc[first_date]
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
        data=group_data.dropna().groupby(["i_lat", "i_lon", "lat", "lon"]).agg({v:ops for v in vars_to_plot}).reset_index()
        #"""
        
        time0=str(group_data.index.values[0]).split("T")[0]
        time1=str(group_data.index.values[-1]).split("T")[0]
        time="%s - %s"%(time0, time1)
        time0=datetime.strptime(str(group_data.index.values[0]), "%Y-%m-%dT%H:%M:%S0.000000000").strftime("%Y-%m-%d %H:%M:%S")
        time1=datetime.strptime(str(group_data.index.values[-1]), "%Y-%m-%dT%H:%M:%S0.000000000").strftime("%Y-%m-%d %H:%M:%S")
        


        ## Obs data of period
        obs_data_period={k:obs_data[k].loc[time0:time1].agg([c for c in ops if isinstance(c, str)] ) for k in obs_data}

        ## Make plots square
        data_plot=pd.concat([grid_should_have[["i_lat", "i_lon", "lat", "lon"]], data]).drop_duplicates(subset=[("i_lat", ""), ("i_lon", "")], keep="last")[data.columns]#data.merge(grid_should_have.drop(columns=["lat", "lon"]), on=[("i_lat", ""), ("i_lon", "")], how="outer")[data.columns]
        data_plot.columns = pd.MultiIndex.from_tuples(data_plot.columns) 
        data_plot.index=np.repeat(time, data_plot.shape[0]) 


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
                if vmin >= vmax:
                    continue
                ## Si variable no tiene valores
                if np.isnan(np.nanmin(z)):
                    continue

                ## If it's time only do the one without anything
                if type_plot=="time-max":
                    if freq=="H":
                        cmap=settings.cmap_cyclic
                        plot_map(x, y, z, varcom, "h", model, time, bounds, obs_data_period, t_Stations_Info, 0, 23, direc=direc, map_type="Square", boolLog=boolLog)
                    elif freq=="D":
                        cmap=settings.cmap_seq
                        plot_map(x, y, z, varcom, "d", model, time, bounds, obs_data_period, t_Stations_Info, 0, np.nanmax(z), direc=direc, map_type="Square", boolLog=boolLog)
                    continue

                    
                
                if 1 in list_maps:
                    plot_map(x, y, z, varcom, units[var], model, time, bounds, obs_data_period, t_Stations_Info, vmin, vmax, direc=direc, map_type="Square", boolLog=boolLog, statistics=statistics[var])
                if 2 in list_maps:
                    plot_map(x, y, z, varcom, units[var], model, time, bounds, obs_data_period, t_Stations_Info, vmin, vmax, direc=direc, map_type="Interpolated", boolLog=boolLog, statistics=statistics[var])
                if 3 in list_maps:
                    plot_map(x, y, z, varcom, units[var], model, time, bounds, obs_data_period, t_Stations_Info, vmin, vmax, direc=direc, map_type="Contour", boolLog=boolLog, statistics=statistics[var])
                


###############################################################
###############################################################
def mainn(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info):
    if IncF.i_Map==0:
        print("...Plotting was not requested.")
        return None
    if len([c for c in IncF.c_ModelNames if c])==0:
        print("...Plotting was requested but no model was read.")
        return None  

    ## Find vars that need to be plotted
    vars_to_plot=[aux.find_equiv_vars([v]) for v in IncF.c_Map_Plots]
    if len(vars_to_plot)==0:
        print("...Plotting was requested but variables are not available.")
        return None    

    models_all = list(set(list([a for key, val in t_ModelData.items() for a in list(val.keys())])))
    models = list(set(list([a for key, val in t_ModelData.items() for a in list(val.keys())])))

    ######################
    for freq in t_ModelData:
        if "Median" in models_all and "Mean" in models_all:
            models = [c for c in models if c not in ["Median", "Mean"]]
            index_name_or = t_ModelData[freq]["Median"].index.names
            index_name = ["i_lat", "i_lon"]
            df_mad=[]
            df_std=[]
            for stardate, enddate in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                #df={k:t_ModelData[freq][k].loc[stardate:enddate, :, :] for k,v in t_ModelData[freq].items()}
                df={k:t_ModelData[freq][k].reorder_levels(['time','i_lat','i_lon']).sort_index().loc["-".join(stardate.split("-")[::-1]):"-".join(enddate.split("-")[::-1]), :, :]for k,v in t_ModelData[freq].items()}

                time0 = df["Median"].index.get_level_values("time")[0]


                ###
                df_concat = pd.concat([np.abs(df[mod] - df["Median"]) for mod in models])
                df_concat = df_concat.groupby(index_name).median()
                df_concat.index = pd.MultiIndex.from_tuples(df_concat.index)
                df_concat.index.names = index_name
                df_concat["time"] = time0
                #t_ModelData[freq]["MAD"] = df_concat.set_index(['time'], append=True).reorder_levels(index_name_or)
                df_mad.append(df_concat.set_index(['time'], append=True).reorder_levels(index_name_or))

                ###
                try:
                    df_concat = pd.concat([(df[mod].drop(columns=["time"]) - df["Mean"].drop(columns=["time"]))**2 for mod in models])
                except:
                    df_concat = pd.concat([(df[mod] - df["Mean"])**2 for mod in models])

                df_concat = np.sqrt(df_concat.groupby(index_name).sum()/len(models))
                df_concat.index  = pd.MultiIndex.from_tuples(df_concat.index)
                df_concat.index.names = index_name
                df_concat["time"] = time0
                #t_ModelData[freq]["STD"] = df_concat.set_index(['time'], append=True).reorder_levels(index_name_or)
                df_std.append(df_concat.set_index(['time'], append=True).reorder_levels(index_name_or))

            t_ModelData[freq]["STD"] = pd.concat(df_std)
            t_ModelData[freq]["MAD"] = pd.concat(df_mad)

    ########
    path=join(IncF.c_FigureDir, "SurfaceMaps")
    makedirs(path, exist_ok = True)       
    #for freq in IncF.c_Map_freq:
    for count, plot_freq in enumerate(IncF.c_Map_freq):
        freq, plot_freq =plot_freq.split(",")
        ## Goes through network asked and finds the filters associated to it
        for index_cb, plot in enumerate(vars_to_plot):
            for model in ["STD", "MAD"] +  list(plot["M"].keys()):
                if "Median" not in models_all or "Mean" not in models_all:
                    if model in ["STD", "MAD"]:
                        continue
                print(f"...Creating plots for  '{IncF.c_Map_Plots[index_cb]}'  with data frequency  '{freq}'  and plot frequency  '{plot_freq}'")
                variables=vars_to_plot[index_cb]
                colorbar=IncF.c_Map_ColBar[index_cb]

                


                path_model=join(path, model)
                makedirs(path_model, exist_ok = True) 
                path_model=join(path_model, f"{freq},{plot_freq}")
                makedirs(path_model, exist_ok = True) 


                ###########
                ## Retrieve model Data needed
                #units=t_ModelData[freq][model].meta.units
                data_resampled=t_ModelData[freq][model].reset_index(level=[1,2]).sort_values(by=["i_lat", "i_lon"]) 

                ## Vars that need to be done 
                try:
                    vars_model=variables["M"][model]
                    ## Add lat and lon data to data set
                    model_info=t_Models_Info[model]
                    filters=t_Model_Filters[model]
                    units=t_ModelData[freq][model].meta.units
                except:
                    vars_model=variables["M"]["Median"]
                    ## Add lat and lon data to data set
                    model_info=t_Models_Info["Median"]             
                    filters=t_Model_Filters["Median"]
                    units=t_ModelData[freq]["Median"].meta.units
                data_resampled=data_resampled.merge(model_info, on=["i_lat", "i_lon"]).set_index(data_resampled.index)

                statistics={v:None for v in vars_model}
                if False and IncF.i_Map_Stats!=0:
                    statistics=aux_maps.prepare_statistics(t_ObsStationData[freq], t_Stations_Info, t_ModelStationData[freq][model],  t_ModelStation_Info, vars_model)



                ##########
                ## filters 
                #filters=t_Model_Filters[model]
                for f in filters:
                    ## DELETE
                    if f!="f0":
                        continue
                    for group, alias in enumerate(filters[f]):
                            indexes_group=list(filters[f][alias]["indexes"].values)
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
                                if f=="f0":
                                    ##DELETE
                                    minx=-110.0
                                    miny=-58.0
                                    maxx=-33.5
                                    maxy=25.0
                                    bounds = minx, miny, maxx, maxy

                            path_filter=join(path_model, "%s-%s"%(f, alias))
                            makedirs(path_filter, exist_ok = True) 
                            print("   ...For filter area:  ", "%s-%s"%(f, alias))
                            ## Do plots
                            direct_to_plot_maps(data_plot_filter, vars_model, units, model, bounds, model_info, t_ObsStationData[freq], t_Stations_Info, plot_freq, colorbar, direc=path_filter, filter="%s-%s"%(f, group), freq=freq, statistics=statistics)
                            plt.close('all')


