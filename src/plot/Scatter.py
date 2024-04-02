# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import numpy as np
import pandas as pd
#
from scipy.stats import linregress, gaussian_kde
#
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import  cm, legend_handler
from matplotlib.lines import Line2D
#
import src.common as aux
import src.plot.aux_plots as aux_plots
#
from os.path import join
from os import makedirs
#
import geopy.distance
from shapely.ops import nearest_points
#
matplotlib.rc('font', family='Arial')



######################################
settings=aux_plots.SCTConfig
## Settings plot
figsize     = settings.imagesize
titlesize   = settings.titlesize
subtitlesize= settings.subtitlesize
ticksize    = settings.ticksize
labelsize   = settings.labelsize
## Legend Location
legendsize      = settings.legendsize
legend_location = settings.legend_location
## Aesthetics
ms         = settings.markersize
color      = settings.color
edgecolor  = settings.edgecolor
cmap       = settings.cmap
linecolor  = settings.linecolor
## Obs vs Obs Settings
smb_x      = settings.smb_x
smb_y      = settings.smb_y
b_closest  = settings.closest

scatter_plots=[int(plot) for plot in str(IncF.i_SCT)]


######################################################
## Finds the closest grid points of stations in df to the models
def find_closest_points(info1, info2):
    #########
    def near(point, pts=info2.geometry.unary_union):
        nearest = info2.geometry == nearest_points(point, pts)[1]
        return info2[nearest].ID.values[0], info2[nearest].geometry.values[0]

    ## First nearest neighbor informations
    info1[['ID_Nearest', "Nearest"]] = pd.DataFrame(info1.apply(lambda row: near(row.geometry), axis=1).values.tolist(), index=info1.index)
    ## Distance
    info1["Distance"]=info1.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.Nearest.y, a.Nearest.x)), axis=1)

    return info1


######################################
## Function that makes all plots
def make_plots_cases(data_all, info_x, info_y, perc, case, plot_freq, t_Stations_Info, t_ModelStation_Info, t_Stations_Filters,  direc=""):
    def get_distance(df1, df2):
        return geopy.distance.geodesic((df1.geometry.y.values[0], df1.geometry.x.values[0]), (df2.geometry.y.values[0], df2.geometry.x.values[0])).kilometers

    g_x, v_x, var_x, u_x = info_x["g_x"], info_x["v_x"], info_x["var_x"], info_x["units_x"]
    g_y, v_y, var_y, u_y = info_y["g_y"], info_y["v_y"], info_y["var_y"], info_y["units_y"]

    units = {"u_x":u_x, "u_y":u_y}

    for col in data_all:
        if not data_all[col].dropna().empty:
            data_all.loc[data_all[col]>np.nanpercentile(data_all[col], perc), col]=np.nan
    
    if case=="Model-Obs":
        #######
        ## Individual Plots
        ## Check which axis has Obs data
        try:
            sta_info=t_Stations_Info.loc[g_x]
            IDs=sta_info["ID"]
            filters=t_Stations_Filters[g_x]
            mod_info=t_ModelStation_Info[g_y].loc[g_x]
            model=g_y
            network=g_x
            var_network=var_x
            var_model=var_y
            x="Obs"

        except:
            sta_info=t_Stations_Info.loc[g_y]
            IDs=sta_info["ID"]
            filters=t_Stations_Filters[g_y]
            mod_info=t_ModelStation_Info[g_x].loc[g_y]
            model=g_x
            network=g_y
            var_network=var_y
            var_model=var_x
            x="Model"


        ## Check only filters or stations
        sta_x=t_Stations_Info.loc[network]
        if IncF.b_SCT_plotstations and not IncF.b_SCT_plotfilters:
            sta_x=sta_x[~sta_x["ID"].str.contains('^f[0-9]-+')]
        if not IncF.b_SCT_plotstations and IncF.b_SCT_plotfilters:
            sta_x=sta_x[sta_x["ID"].str.contains('^f[0-9]-+')]

        #######
        ## Individual Station Plots
        for ID in sta_x["ID"]:
            ID_mod=mod_info[mod_info["ID"]==ID]
            ID_sta=sta_info[sta_info["ID"]==ID]
            data=data_all[["%s %s"%(var_x, ID), "%s %s"%(var_y, ID)]]
            #data_all[[("%s %s"%(var_x, ID)).replace("Model ", "").replace("Obs ", ""), ("%s %s"%(var_y, ID)).replace("Model ", "").replace("Obs ", "")]]
            data.dropna(inplace=True)
            if data.empty:
                continue

            distance="Distance: %.2f km"%get_distance(ID_mod, ID_sta)
            loc_obs  ="Obs    :  %s    %s"%(ID_sta.location.values[0], ID_sta.Nombre.values[0])
            loc_model="Model :  %s    %s"%(ID_mod.location.values[0], model)

            if x=="Obs":
                v_x="Obs %s"%var_x
                v_y="Model %s"%var_y
                loc_x=loc_obs
                loc_y=loc_model

            elif x=="Model":
                v_x="Model %s"%var_x
                v_y="Obs %s"%var_y
                loc_x=loc_model
                loc_y=loc_obs
                
            locations=[loc_x, loc_y, distance]


            ## Plots with the frequency asked
            if plot_freq=="Whole":
                data["groupby"]=1
                grouped_data=data.groupby("groupby")
            else:
                grouped_data=data.groupby(pd.Grouper(freq=plot_freq))
                  
            #########
            list_dates=list(grouped_data.groups.keys()) 
            for group in list_dates: 
                try:
                    group_data=grouped_data.get_group(group)
                except:
                    group_data=pd.DataFrame()
                if group_data.empty:
                    continue
            
                group_data=group_data[[c for c in group_data.columns if c !="groupby"]]
                timeframe=[group_data.index[0].strftime("%d-%m-%Y"), group_data.index[-1].strftime("%d-%m-%Y")]

            #for timeframe in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                if 1 in scatter_plots:
                    plot_scatter_standard(group_data, ID, v_x, v_y, timeframe, units, direc=direc, locations=locations)   
                if 2 in scatter_plots:
                    plot_scatter_density(group_data, ID, v_x, v_y, timeframe, units, direc=direc, locations=locations)
                if 3 in scatter_plots:
                    plot_scatter_linfit(group_data, ID, v_x, v_y, timeframe, units, direc=direc, locations=locations)
                if 4 in scatter_plots:
                    plot_scatter_densitylinfit(group_data, ID, v_x, v_y, timeframe, units, direc=direc, locations=locations)

        #######
        ## Filter Plots
        network_columns=[c for c in data_all.columns if aux.clean_network(network) in c]
        model_columns=[c for c in data_all.columns if model in c]
        for f in t_Stations_Filters[network]:
            if f in ["f0"]:
                continue
            for index, group in enumerate(t_Stations_Filters[network][f]):
                cols_network=[c for c in network_columns if c.split()[-1] in group]
                cols_model=[c for c in model_columns if c.split()[-1] in group]
                ## Dataframes with averages
                mod_average=data_all[cols_model].mean(axis=1)
                net_average=data_all[cols_network].mean(axis=1)
                data=pd.concat([mod_average, net_average], axis=1, keys=[var_model+" ", var_network+" "])
                data.dropna(inplace=True)
                if data.empty:
                    continue
                
                #########
                ## Generate Texts to add to subplot
                ## Filter information
                if f=="f1":
                    filter_info="Filter by region mask"
                elif f=="f2":
                    filter_info=f"Filter by square region  {IncF.c_SqrRegionsAlias[index]}:  {IncF.i_SqrRegions[index]}"
                elif f=="f3":
                    index_group=aux.get_index_network(network)
                    filter_info=f"Filter by list of stations. Alias: {IncF.c_Stations_Obs_Aliases[index_group][index]}"
                elif f=="f4":
                    filter_info=f'Filter by Chilean regions:  {" , ".join(IncF.c_Regions[index])}'
                elif f=="f5":
                    filter_info=f'Filter by Chilean comunas:  {" , ".join(IncF.c_Comunas[index])}'
                elif f=="f6":
                    filter_info=f"Filter by station's elevation:  {IncF.i_Elevation[index]}"

                ## Stations being averaged
                stations_info=" \t"+"\n \t".join(sta_info[sta_info["ID"].isin(group)][["Nombre", "location"]].apply(lambda  df: "%s    %s "%(df.location, df.Nombre), axis=1).dropna().values)

                gen_info="Model    : %s    \nNetwork : %s\n"%(model, aux.clean_network(network))               
                locations=[gen_info, filter_info, stations_info.replace("\t", "     ")] 
                ## Plots with the frequency asked
                if plot_freq=="Whole":
                    data["groupby"]=1
                    grouped_data=data.groupby("groupby")
                else:
                    grouped_data=data.groupby(pd.Grouper(freq=plot_freq))
                      
                #########
                list_dates=list(grouped_data.groups.keys()) 
                for time in list_dates:
                    try:
                        group_data=grouped_data.get_group(time)
                        group_data=group_data[[c for c in group_data.columns if c !="groupby"]]
                        timeframe=[group_data.index[0].strftime("%d-%m-%Y"), group_data.index[-1].strftime("%d-%m-%Y")]
                    except: 
                        group_data=pd.DataFrame()

                    if group_data.empty:
                        continue

                #for timeframe in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                    if 1 in scatter_plots:
                        plot_scatter_standard(data, "%s_Ave(%s)"%(f,"-".join(group)), v_x, v_y, timeframe, units, direc=direc, locations=locations, ID_x="", ID_y="")   
                    if 2 in scatter_plots:
                        plot_scatter_density(data, "%s_Ave(%s)"%(f,"-".join(group)), v_x, v_y, timeframe, units, direc=direc, locations=locations, ID_x="", ID_y="")
                    if 3 in scatter_plots:
                        plot_scatter_linfit(data, "%s_Ave(%s)"%(f,"-".join(group)), v_x, v_y, timeframe, units, direc=direc, locations=locations, ID_x="", ID_y="")
                    if 4 in scatter_plots:
                        plot_scatter_densitylinfit(data, "%s_Ave(%s)"%(f,"-".join(group)), v_x, v_y, timeframe, units, direc=direc, locations=locations, ID_x="", ID_y="")

    else:#elif case=="Obs-Obs":
        if case=="Obs-Obs":
            sta_x=t_Stations_Info.loc[g_x]
            sta_y=t_Stations_Info.loc[g_y]
        else:
            sta_x=t_ModelStation_Info[g_x]
            sta_y=t_ModelStation_Info[g_y]

        ## Check only filters or stations
        if IncF.b_SCT_plotstations and not IncF.b_SCT_plotfilters:
            sta_x=sta_x[~sta_x["ID"].str.contains('^f[0-9]-+')]
            sta_y=sta_y[~sta_y["ID"].str.contains('^f[0-9]-+')]
        if not IncF.b_SCT_plotstations and IncF.b_SCT_plotfilters:
            sta_x=sta_x[sta_x["ID"].str.contains('^f[0-9]-+')]
            sta_y=sta_y[sta_y["ID"].str.contains('^f[0-9]-+')]

        closest_x=find_closest_points(sta_x, sta_y)
        closest_y=find_closest_points(sta_y, sta_x)
 

        #######
        ## Individual Plots
        #breakpoint()
        for ID_x in sta_x["ID"]:
            for ID_y in sta_y["ID"]:
                ## If its the same station of the same network skip
                if (ID_x==ID_y) and (g_y==g_x) and (var_x==var_y):
                    continue
                ## If they are different networks and I only want the closest station then check if they are the same and skip otherwise
                if b_closest:
                    if ID_y!=closest_x[closest_x["ID"]==ID_x]["ID_Nearest"].values[0]:
                        continue
                try:
                    data=data_all[["%s %s"%(var_x, ID_x), "%s %s"%(var_y, ID_y)]]
                    data.dropna(inplace=True)
                    if data.empty:
                        continue
                except KeyError:
                    continue

                info_x=sta_x[sta_x["ID"]==ID_x]
                info_y=sta_y[sta_y["ID"]==ID_y]

                loc_x=smb_x + u"  %s    %s   %s"%(info_x.location.values[0], var_x, info_x.Nombre.values[0])
                loc_y=smb_y + u"  %s    %s   %s"%(info_y.location.values[0], var_y, info_y.Nombre.values[0])

                distance="Distance: %.2f km"%get_distance(info_x, info_y)

                locations=[loc_x, loc_y, distance]
                for timeframe in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                    if 1 in scatter_plots:
                        plot_scatter_standard(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
                    if 2 in scatter_plots:
                        plot_scatter_density(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
                    if 3 in scatter_plots:
                        plot_scatter_linfit(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
                    if 4 in scatter_plots:
                        plot_scatter_densitylinfit(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)

        """
        #######
        ## Individual Plots
        for ID_y in sta_y["ID"]:
            for ID_x in sta_x["ID"]:
                ## If its the same station of the same network skip
                if (ID_x==ID_y) and (g_y==g_x):
                    continue
                ## If they are different networks and I only want the closest station then check if they are the same and skip otherwise
                if b_closest:
                    if ID_x!=closest_y[closest_y["ID"]==ID_y]["ID_Nearest"].values[0]:
                        continue
                        
                data=data_all[["%s %s"%(var_x, ID_x), "%s %s"%(var_y, ID_y)]]
                data.dropna(inplace=True)
                if data.empty:
                    continue
                info_x=sta_x[sta_x["ID"]==ID_x]
                info_y=sta_y[sta_y["ID"]==ID_y]

                loc_x=smb_x + u"  %s    %s   %s"%(info_x.location.values[0], var_x, info_x.Nombre.values[0])
                loc_y=smb_y + u"  %s    %s   %s"%(info_y.location.values[0], var_y, info_y.Nombre.values[0])

 
                distance="Distance: %.2f km"%get_distance(info_x, info_y)

                locations=[loc_x, loc_y, distance]

                for timeframe in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                    if 1 in scatter_plots:
                        plot_scatter_standard(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
                    if 2 in scatter_plots:
                        plot_scatter_density(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
                    if 3 in scatter_plots:
                        plot_scatter_linfit(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
                    if 4 in scatter_plots:
                        plot_scatter_densitylinfit(data, "%s-%s"%(ID_x, ID_y), var_x, var_y, timeframe, units, direc=direc, ID_x=ID_x, ID_y=ID_y, locations=locations)
        """

######################################
## Find var, network/model and corresponding data to one variable in the string format M_*_* O_*_*
def prepare_data(var, obs, models, freq, t_ObsStationData, t_ModelStationData):
    group_, var_=aux.find_equiv_var(var, obs, models, databases={"O":[], "M":[]}, variables={"O":[], "M":[]})
    group_=aux.clean_dic_lists(group_)
    var_=aux.clean_dic_lists(var_)

    if len(group_)==0 or len(var_)==0:
        return {}, [], [], None, {}

    type_=list(group_.keys())[0]
    groups_=group_[type_]
    varias_=var_[type_]
    if type_=="O":
        data_=t_ObsStationData[freq]
        units_={k:v.meta.units for k,v in t_ObsStationData[freq].items() if k in groups_}
    else:
        data_=t_ModelStationData[freq]
        units_={km:{k:v.meta.units} for km, vm in t_ModelStationData[freq].items() for k, v in t_ModelStationData[freq][km].items() if km in groups_}
        units_={km:{kg:vg.meta.units for kg, vg in t_ModelStationData[freq][km].items()} for km, vm in t_ModelStationData[freq].items() if km in groups_} 

    return data_, groups_, varias_, type_, units_

######################################
## Get name of columns in dataframe that contains data
def get_name_df_cols(ID, ID_x, ID_y, var_x, var_y):
    if (ID_x is None) or (ID_y is None):
        ID_xx=ID
        ID_yy=ID
    else:
        ID_xx=ID_x
        ID_yy=ID_y

    xname="%s %s"%(var_x.replace("Obs ", "").replace("Model ", ""), ID_xx)
    yname="%s %s"%(var_y.replace("Obs ", "").replace("Model ", ""), ID_yy)
    return xname, yname


######################################
## Obtain merged data for the x and y variables
def get_data(data_x, data_y, info_x, info_y, case):
    g_x=info_x["g_x"]
    v_x=info_x["v_x"]
    g_y=info_y["g_y"]
    v_y=info_y["v_y"]

    if case=="Model-Obs":
        try:
            ## Mod Data as Y
            d_y=data_y[g_y][g_x][v_y]
            var_y="%s %s"%(v_y, g_y)
            d_y.columns=["%s %s"%(var_y, c) for c in d_y.columns]
            ## Obs Data as X
            d_x=data_x[g_x][v_x]
            g_x=aux.clean_network(g_x)
            var_x="%s %s"%(v_x, g_x) 
            d_x.columns=["%s %s"%(var_x, c) for c in d_x.columns]
        except:
            ## Mod Data as X
            d_x=data_x[g_x][g_y][v_x]
            var_x="%s %s"%(v_x, g_x)
            d_x.columns=["%s %s"%(var_x, c) for c in d_x.columns]
            ## Obs Data as Y
            d_y=data_y[g_y][v_y]
            g_y=aux.clean_network(g_y)
            var_y="%s %s"%(v_y, g_y)
            d_y.columns=["%s %s"%(var_y, c) for c in d_y.columns]
        string="%s-%s"%(g_x, g_y)
    elif case=="Obs-Obs":
        ## Obs Data as X
        d_x=data_x[g_x][v_x]
        g_x=aux.clean_network(g_x)
        var_x="%s %s"%(v_x, g_x)
        d_x.columns=["%s %s"%(var_x, c) for c in d_x.columns]
        ## Obs Data as Y
        d_y=data_y[g_y][v_y]
        g_y=aux.clean_network(g_y)
        var_y="%s %s"%(v_y, g_y)
        d_y.columns=["%s %s"%(var_y, c) for c in d_y.columns]
        if g_x==g_y:
            string=g_x
        else:
            string="%s-%s"%(g_x, g_y)
    elif case=="Model-Model":
        net=list(data_x[g_x].keys())[0]
        ## Obs Data as X
        d_x=data_x[g_x][net][v_x]
        var_x="%s %s"%(v_x, g_x)
        d_x.columns=["%s %s"%(var_x, c) for c in d_x.columns]
        ## Obs Data as Y
        d_y=data_y[g_y][net][v_y]
        var_y="%s %s"%(v_y, g_y)
        d_y.columns=["%s %s"%(var_y, c) for c in d_y.columns]
        if g_x==g_y:
            string=g_x
        else:
            string="%s-%s"%(g_x, g_y)
    ## If same network wants to be compared save one
    if d_x.equals(d_y):
        data_all=d_x
    ## Else keep all information
    else:
        data_all=d_x.merge(d_y, left_index=True, right_index=True)
    ## Delete
    #data_all=data_all.applymap(lambda l: l if not np.isnan(l) else np.random.choice([1, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5])) 

    return string, data_all, var_x, var_y


###############################################################
###############################################################
def prepare_plot(ID, ID_x, ID_y, var_x, var_y, units, locations, fig=None, ax=None, fig_meta=None, ax_meta=None):
    ## Retrieve Figures
    ## Figure Map
    fig=fig or plt.subplots(1, figsize=figsize)
    ax=ax or plt.gca()

    ## Defines how many ticks to put down
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    
    ## Figure Metadata
    fig_meta=fig_meta or plt.subplots(1, figsize=figsize)
    ax_meta=ax_meta or plt.gca()

    ax_meta.set_axis_off()
    
    ## Titles and subtitles
    subtitle='\n'.join(map(str, locations))
    ax_meta.set_title(subtitle, fontsize=subtitlesize, y=0.975, loc="left", fontdict={'verticalalignment':"top"})

    ## Formatting settings
    ax.xaxis.set_tick_params(labelsize=ticksize)
    ax.yaxis.set_tick_params(labelsize=ticksize)
    ax.xaxis.label.set_size(labelsize)
    ax.yaxis.label.set_size(labelsize)


    ## Add 1:1 line if asked
    if IncF.i_SCT_diag==1:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        gmin=min(xmin, ymin)
        gmax=max(xmax, ymax)
        grange=gmax-gmin
        gmin=gmin-grange*0.04
        gmax=gmax+grange*0.04

        x=np.linspace(gmin-0.1*grange, gmax+0.1*grange, 100) 
        ax.plot(x, x, 'k--', zorder=-100)

        if xmax<=0.07*gmax:
            ax.set_xlim(gmin, 0.125*gmax)
        elif xmax<=0.12*gmax:
            ax.set_xlim(gmin, 0.25*gmax)
        elif xmax<=0.25*gmax:
            ax.set_xlim(gmin, 0.5*gmax)
        else:
            ax.set_xlim(gmin, gmax)

        if ymax<=0.07*gmax:
            ax.set_ylim(gmin, 0.125*gmax)
        elif ymax<=0.12*gmax:
            ax.set_ylim(gmin, 0.25*gmax)
        elif ymax<=0.25*gmax:
            ax.set_ylim(gmin, 0.5*gmax)
        else:
            ax.set_ylim(gmin, gmax)


    ## Axis labels

    ## Model-Obs. ID only shared between model and obs, ID_x and ID_y are nont. 
    if (ID_x is None) or (ID_y is None):
        group=var_x if "Obs" in var_x else var_y
        group=group.split()[-1]
        clean_var_x=var_x.split()[1]
        clean_var_y=var_y.split()[1]
        ## Titles
        if clean_var_x==clean_var_y:
            fig_meta.suptitle('%s   Obs vs Model'%(clean_var_x), fontsize=titlesize, y=1)
        else:
            fig_meta.suptitle('%s vs %s'%(clean_var_x, clean_var_y), fontsize=titlesize, y=1)

        ## Axis Label
        xlabel="%s [%s]"%(" ".join(var_x.split()[:2]), units["u_x"])
        ylabel="%s [%s]"%(" ".join(var_y.split()[:2]), units["u_y"])

    elif ID_x=="" or ID_y=="":
        ## Averages
        if ID[:3] in ["f1_", "f2_", "f3_", "f4_", "f5_", "f6_"]:
            ## Axis Label
            xlabel="%s [%s]"%(var_x, units["u_x"])
            ylabel="%s [%s]"%(var_y, units["u_y"])  
        else:    
            group=[var_x.replace("Model ", "").replace("Obs ", "").split()[-1], var_y.replace("Model ", "").replace("Obs ", "").split()[-1]]
            group=[g for g in group if g in IncF.c_ObsNetwork][0]

            clean_var_x=var_x.split()[1]
            clean_var_y=var_y.split()[1]
            if group in var_x:
                var_x="Obs "+ clean_var_x
                var_y="Model "+ clean_var_y

            elif group in var_y:
                var_y="Obs "+ clean_var_y
                var_x="Model "+ clean_var_x

            ## Titles
            if clean_var_x==clean_var_y:
                fig_meta.suptitle('%s  Obs vs Model'%(clean_var_x), fontsize=titlesize, y=1)
            else:
                fig_meta.suptitle('%s vs %s'%(clean_var_x, clean_var_y), fontsize=titlesize, y=1)

            ## Axis Label
            xlabel="%s [%s]"%(" ".join(var_x.split()[:2]), units["u_x"])
            ylabel="%s [%s]"%(" ".join(var_y.split()[:2]), units["u_y"])      
    else:
        ## Titles 
        fig_meta.suptitle('%s vs %s'%(var_x,var_y), fontsize=titlesize, y=1)

        ## Axis Label
        units_x=aux.read_units(var_x.split()[1])
        units_y=aux.read_units(var_y.split()[1])
        xlabel="%s  %s [%s]"%(smb_x, var_x.split()[0], units["u_x"])
        ylabel="%s  %s [%s]"%(smb_y, var_y.split()[0], units["u_y"])



    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig, ax



#####################################
## Plot scatter
def plot_scatter_standard(data, ID, var_x, var_y, timeframe, units, direc="", ID_x=None, ID_y=None, locations=[]):

    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data.loc[date1:date2].dropna()
    if data.empty:
        return None


    xname, yname=get_name_df_cols(ID, ID_x, ID_y, var_x, var_y)
    name=("SCT__%s___%s_vs_%s___%s-%s"%(ID.strip(), xname.strip(), yname.strip(), date1, date2)).replace(" ", "-")

    #########
    ### Plot
    ########## METADATA
    fig_meta, ax_meta = plt.subplots(1, figsize=figsize) 
    ########## MAIN
    fig, ax = plt.subplots(1, figsize=figsize) 
    ## Plot
    data.plot.scatter(x=xname, y=yname, c=color, ax=ax, s=ms, edgecolors=edgecolor)
    ## Settings
    fig, ax=prepare_plot(ID, ID_x, ID_y, var_x, var_y, units, locations, fig=fig, ax=ax, fig_meta=fig_meta, ax_meta=ax_meta)


    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)

    return None

#####################################
## Plot scatter density
def plot_scatter_density(data, ID, var_x, var_y, timeframe, units, direc="", ID_x=None, ID_y=None, locations=[]):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data.loc[date1:date2].dropna()
    if data.empty:
        return None
 
    xname, yname=get_name_df_cols(ID, ID_x, ID_y, var_x, var_y)
    name=("SCT-density__%s___%s_vs_%s___%s-%s"%(ID.strip(), xname.strip(), yname.strip(), date1, date2)).replace(" ", "-")

    ## Prepare data
    x=data[xname].values
    y=data[yname].values
    if len(x)==1 or len(y)==1:
        return None

    xy = np.vstack([x, y])
    data["Density"] = gaussian_kde(xy)(xy)
    data.sort_values(by="Density", inplace=True)

    #########
    ### Plot
    ########## METADATA
    fig_meta, ax_meta = plt.subplots(1, figsize=figsize) 
    ########## MAIN
    fig, ax = plt.subplots(1, figsize=figsize)         
    data.plot.scatter(x=xname, y=yname, c="Density", cmap=cmap, ax=ax, s=ms, colorbar=False)
    ## Settings
    fig, ax=prepare_plot(ID, ID_x, ID_y, var_x, var_y, units, locations, fig=fig, ax=ax, fig_meta=fig_meta, ax_meta=ax_meta)
    
    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)


    return None


#####################################
## Plot scatter line fit
def plot_scatter_linfit(data, ID, var_x, var_y, timeframe, units, direc="", ID_x=None, ID_y=None, locations=[]):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data.loc[date1:date2].dropna()
    if data.empty:
        return None

    xname, yname=get_name_df_cols(ID, ID_x, ID_y, var_x, var_y)
    name=("SCT-linefit__%s___%s_vs_%s___%s-%s"%(ID.strip(), xname.strip(), yname.strip(), date1, date2)).replace(" ", "-")

    #########
    ### Figure
    ########## METADATA
    fig_meta, ax_meta = plt.subplots(1, figsize=figsize) 
    ########## MAIN
    fig, ax = plt.subplots(1, figsize=figsize)         

    ## Scatter Plot
    data.plot.scatter(x=xname, y=yname, c=color, ax=ax, s=ms, edgecolors=edgecolor)
    xmin, xmax=ax.get_xlim()

    ## Line
    slope, intercept, rvalue, pvalue, stderr = linregress(data[xname], data[yname])
    sns.regplot(x=xname, y=yname, data=data, line_kws={'color': linecolor}, scatter_kws={'s':ms, "color":color, "edgecolors":edgecolor}, ax=ax)


    ## Settings
    fig, ax=prepare_plot(ID, ID_x, ID_y, var_x, var_y, units, locations, fig=fig, ax=ax, fig_meta=fig_meta, ax_meta=ax_meta)  
    fig.legend(labels=['Points', 'Line Fit'], loc=legend_location, title="R:%.2f  P:%.2f  SE:%.2f"%(rvalue, pvalue, stderr), frameon=True, framealpha=0.25, facecolor="lightgray", bbox_to_anchor=(0, 0, 0.875, 0.87), fontsize=int(0.8*subtitlesize), title_fontsize=subtitlesize)

    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)
    return None


#####################################
## Plot scatter line fit + density
def plot_scatter_densitylinfit(data, ID, var_x, var_y, timeframe, units, direc="", ID_x=None, ID_y=None, locations=[]):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data.loc[date1:date2].dropna()
    if data.empty:
        return None

    xname, yname=get_name_df_cols(ID, ID_x, ID_y, var_x, var_y)
    name=("SCT-denlfit__%s___%s_vs_%s___%s-%s"%(ID.strip(), xname.strip(), yname.strip(), date1, date2)).replace(" ", "-")

    ## Prepare data
    x=data[xname].values
    y=data[yname].values
    if len(x)==1 or len(y)==1:
        return None
    xy = np.vstack([x, y])
    data["Density"] = gaussian_kde(xy)(xy)
    data.sort_values(by="Density", inplace=True)

    #########
    ### Figure
    ########## METADATA
    fig_meta, ax_meta = plt.subplots(1, figsize=figsize) 
    ########## MAIN
    fig, ax = plt.subplots(1, figsize=figsize)         
    ## Scatter Plot
    data.plot.scatter(x=xname, y=yname, c="Density", cmap=cmap, ax=ax, s=ms, colorbar=False)
    xmin, xmax=ax.get_xlim()

    ## Line
    slope, intercept, rvalue, pvalue, stderr = linregress(data[xname], data[yname])
    sns.regplot(x=xname, y=yname, data=data, line_kws={'color': linecolor}, scatter=False, ax=ax)
        
    ## Settings
    fig, ax=prepare_plot(ID, ID_x, ID_y, var_x, var_y, units, locations, fig=fig, ax=ax, fig_meta=fig_meta, ax_meta=ax_meta)
    colormap = cm.get_cmap(cmap)
    ms_new=0.12*ms
    handles=[(Line2D([0], [0], marker='o', color=colormap(0.1), markersize=ms_new), Line2D([0], [0], marker='o', color=colormap(0.5), markersize=ms_new), Line2D([0], [0], marker='o', color=colormap(0.9), markersize=ms_new)), Line2D([0], [0], color=linecolor)]
    labels=["Density", "Line Fit"]
    fig.legend(handles=handles, labels=labels, loc=legend_location, title="R:%.2f  P:%.2f  SE:%.2f"%(rvalue, pvalue, stderr), frameon=True, framealpha=0.25, facecolor="lightgray", bbox_to_anchor=(0, 0, 0.875, 0.87), fontsize=int(0.8*subtitlesize), title_fontsize=subtitlesize, scatterpoints=3, handler_map = {tuple: legend_handler.HandlerTuple(None)})

    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)

    return None


###############################################################
###############################################################


#####################################
## TODO: Fix code to include timeseries of 3d models
def main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters):
    if IncF.i_SCT==0:
        print("...Plotting was not requested.")
        return None
    if t_Stations_Info.empty:
        print("...No stations were read. Skipping this step.")
        return None     

    ## Check models and obs that ended up being read
    models=[c for c in IncF.c_ModelNames if c]
    obs=[c for c in IncF.c_ObsNetwork if c]   

    path_base=join(IncF.c_FigureDir, "Scatter")
    makedirs(path_base, exist_ok = True) 
    for count, plot_freq in enumerate(IncF.c_SCT_freq):
        freq, plot_freq =plot_freq.split(",")

        path_freq=join(path_base, f"{freq},{plot_freq}")
        makedirs(path_freq, exist_ok = True) 
        for perc in IncF.i_SCT_perc:
            path=join(path_freq, "Per%d"%(perc))
            makedirs(path, exist_ok = True) 
            for vars_plot in IncF.c_SCT_Plots:
                print("...Plotting:  ", vars_plot, "  for freq   ", freq)
                amount_net=len([var for var in vars_plot if "O" in var])
                if amount_net==1:
                    case="Model-Obs"
                elif amount_net==2:
                    case="Obs-Obs"
                else:
                    case="Model-Model"


                ####################
                ## Finds vars that need to be plotted and total vars requested
                data_x, groups_x, varias_x, type_x, units_x = prepare_data(vars_plot[0], obs, models, freq, t_ObsStationData, t_ModelStationData)
                data_y, groups_y, varias_y, type_y, units_y = prepare_data(vars_plot[1], obs, models, freq, t_ObsStationData, t_ModelStationData)
                
                if len(data_x)==0 or len(data_y)==0:
                    if case=="Model-Obs":
                        print("   ...Skipping. Network requested did not have data.")
                    if case=="Obs-Obs":
                        print("   ...Skipping. At least one of the networks requested did not have data.")
                    if case=="Model-Model":
                        print("   ...Skipping. At least one of the models requested did not have data.")
                    continue         
                               
                for ix, g_x in enumerate(groups_x):
                    for v_x in varias_x[ix]:
                        for iy, g_y in enumerate(groups_y):
                            for v_y in varias_y[iy]:
                                if case=="Model-Model":
                                    u_x=units_x[g_x][list(units_x[g_x].keys())[0]][v_x]
                                    u_y=units_y[g_y][list(units_y[g_y].keys())[0]][v_y]
                                else:
                                    ## Get data for specific variables
                                    try:
                                        ## Axis X is model
                                        u_x=units_x[g_x][g_y][v_x]
                                    except:
                                        ## Axis X is network
                                        u_x=units_x[g_x][v_x] 
                                    
                                    try:
                                        ## Axis Y is model
                                        u_y=units_y[g_y][g_x][v_y] 
                                    except:
                                        ## Axis Y is network
                                        u_y=units_y[g_y][v_y] 

                                info_x={"g_x":g_x, "v_x":v_x, "units_x":aux_plots.get_pretty_units(u_x)}
                                info_y={"g_y":g_y, "v_y":v_y, "units_y":aux_plots.get_pretty_units(u_y)}
                                string, data_all, var_x, var_y=get_data(data_x, data_y, info_x, info_y, case)


                                ## Create Directory
                                direc=join(path, string)
                                makedirs(direc,  exist_ok = True) 
                                makedirs(join(direc, "metadata"),  exist_ok = True) 


                                info_x["var_x"]=var_x
                                info_y["var_y"]=var_y


                                ##############
                                ## make plots             
                                make_plots_cases(data_all, info_x, info_y, perc, case, plot_freq, t_Stations_Info, t_ModelStation_Info, t_Stations_Filters, direc=direc)
