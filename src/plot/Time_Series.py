# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import numpy as np
import pandas as pd
#
import matplotlib.pyplot as plt
#
import src.common as aux
import src.plot.aux_plots as aux_plots
from src.plot.aux_timeseries import *

#
from os.path import join
from os import makedirs
#
import geopy.distance
from shapely.ops import nearest_points



TS_plots=[int(p) for p in str(IncF.i_TS)]
bool_plotstations=IncF.b_TS_plotstations
bool_plotfilters=IncF.b_TS_plotfilters
######################################
settings=aux_plots.TSConfig
##
color1, color2=settings.colors
figsize = settings.imagesize
ms = settings.markersize
ls = settings.linestyle
legend_location = settings.legend_location
nele_legend=settings.nele_legend
linewidth_highlight=settings.linewidth_highlight
linewidth_all=settings.linewidth_all
color_highlight=settings.color_highlight
color_highlight_obs=settings.color_highlight_obs
if IncF.b_TS_hlmodel:
    ls=ls[:len(IncF.c_TS_hlmodel)]+ls
    color2=color_highlight[:len(IncF.c_TS_hlmodel)]+color2
if IncF.b_TS_hlobs:
    ls=ls[:len(IncF.c_TS_hlobs)]+ls
    color1=color_highlight_obs[:len(IncF.c_TS_hlobs)]+color1
    
titlesize = settings.titlesize
ticksize = settings.ticksize
labelsize = settings.labelsize
overlap_percentage = settings.overlap_percentage
overlap_twovar = settings.overlap_percentage
legendsize = settings.legendsize
b_closest = settings.closest
min_perc_window = settings.percwindow
bool_center = settings.center
b_showshared = settings.b_showshared
subtitlesize= settings.subtitlesize



#####################################
###################################
## Plot Two Sides var
def plot_basic(name, data_or, color_dict, timeframe, var1, var2, units, direc="", legend_title="", highlight=""):
    ## Make sure data is for timeframe asked
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data_or.loc[date1:date2]

    ## Check wether it's the same variable that is being plotted
    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    equiv1=[k for  k, v in variables_equivalences.items() if var1 in v][0]
    equiv2=[k for  k, v in variables_equivalences.items() if var1 in v][0]

    if equiv1==equiv2:
        name="%s___%s___%s-%s"%(name, var1, date1, date2)
        ## Prepare data
        case, color1_0, color2_0,  data_1, data_2 = prepare_data_one_var(data_or, color_dict)
        BoolTwoVars=False
        if case is None:
            return None 
    else:
        BoolTwoVars=True
        name="%s___%s-%s___%s-%s"%(name, var1, var2, date1, date2)

        ## Separate datacolor_dict
        data_1=data[[v for v in data.columns if var1 in v ]] 
        data_2=data[[v for v in data.columns if var2 in v ]] 
        color1_0=color_dict[data_1.columns[0]]
        color2_0=color_dict[data_2.columns[0]]




    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  

    y1label, y2label, new_legend = prepare_labels(data.columns, data_1.columns, data_2.columns, var1, var2, units)

    prepare_metadata(fig_meta, ax_meta, legend_title, var1, var2)



    ## Plot
    if not data_1.empty:
        data_1.plot(color=[color_dict[c] for c in data_1.columns], ax=ax1, style=ls, ms=ms, lw=linewidth_all)
    if not data_2.empty:
        data_2.plot(color=[color_dict[c] for c in data_2.columns], ax=ax2, style=ls, ms=ms, lw=linewidth_all)


    ## Legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    ## Both vars werent plotted
    if len(handles2)==0 and len(handles1)==0:
        fig.clf()
        plt.close(fig)
        fig_meta.clf()
        plt.close(fig_meta)
        plt.close('all')
        return None

    ## Save image because everything was plotted
    handles= handles1 + handles2
    labels = labels1  + labels2
    if not data_1.empty:
        ax1.get_legend().remove()
    if not data_2.empty:
        ax2.get_legend().remove()


    align_xaxis(ax1)
    BoolOneAxis = align_yaxis(ax1, ax2, BoolTwoVars=BoolTwoVars, var1=var1, var2=var2)
    if BoolOneAxis:
        ax1.set_ylabel(y1label)
    else:
        ax1.tick_params(axis='y', labelcolor=color1_0)
        ax1.set_ylabel(y1label, color=color1_0)
        ax2.tick_params(axis='y', labelcolor=color2_0)
        ax2.set_ylabel(y2label, color=color2_0)

    if IncF.b_TS_hlmodel:
        if len(labels2)!=0:
            highlight_first_line(ax2, labels2, IncF.c_TS_hlmodel, case="Model")
    if IncF.b_TS_hlobs:
        if len(labels1)!=0:
            highlight_first_line(ax1, labels1, IncF.c_TS_hlobs, case="Obs")

    do_legend(fig, ax1, handles, new_legend, legend_title)

    
    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)
    plt.close('all')
    return None


###################################
## Plot Two Sides var
def plot_rolling(name, data_or, color_dict, timeframe, var1, var2, units, direc="", legend_title="", window=1, highlight=""):
    ## Make sure data is for timeframe asked
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data_or.loc[date1:date2]

    ## Check wether it's the same variable that is being plotted
    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    equiv1=[k for  k, v in variables_equivalences.items() if var1 in v][0]
    equiv2=[k for  k, v in variables_equivalences.items() if var1 in v][0]

    if equiv1==equiv2:
        name="%s___%s_Rolling%s___%s-%s"%(name, var1, window, date1, date2)
        ## Prepare data
        case, color1_0, color2_0,  data_1, data_2 = prepare_data_one_var(data_or, color_dict)
        BoolTwoVars=False
        if case is None:
            return None 
    else:
        BoolTwoVars=True
        name="%s___%s-%s_Rolling%s___%s-%s"%(name, var1, var2, window, date1, date2)
        ## Separate datacolor_dict
        data_1=data[[v for v in data.columns if var1 in v ]] 
        data_2=data[[v for v in data.columns if var2 in v ]] 
        color1_0=color_dict[data_1.columns[0]]
        color2_0=color_dict[data_2.columns[0]]



    data_1R=data_1.rolling(window=window,  center=bool_center, min_periods=int(min_perc_window*window/100)).mean().dropna()
    data_2R=data_2.rolling(window=window,  center=bool_center, min_periods=int(min_perc_window*window/100)).mean().dropna()

    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  


    y1label, y2label, new_legend = prepare_labels(data.columns, data_1.columns, data_2.columns, var1, var2, units)



    prepare_metadata(fig_meta, ax_meta, legend_title, var1, var2)



    ## Plot
    handles1, labels1 = [], []
    if not data_1R.empty:
        data_1R.plot(color=[color_dict[c] for c in data_1R.columns], ax=ax1, style=ls, ms=ms)
        handles1, labels1 = ax1.get_legend_handles_labels()
        ax1.get_legend().remove()

    handles2, labels2 = [], []
    if not data_2R.empty:
        data_2R.plot(color=[color_dict[c] for c in data_2R.columns], ax=ax2, style=ls, ms=ms)
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.get_legend().remove()

    

    ## Both vars werent plotted
    if len(handles2)==0 and len(handles1)==0:
        fig.clf()
        plt.close(fig)
        fig_meta.clf()
        plt.close(fig_meta)
        plt.close('all')
        return None

    ## Save image because everything was plotted
    handles= handles1 + handles2
    labels = labels1  + labels2

    align_xaxis(ax1)
    BoolOneAxis = align_yaxis(ax1, ax2, BoolTwoVars=BoolTwoVars, var1=var1, var2=var2)
    if BoolOneAxis:
        ax1.set_ylabel(y1label)
    else:
        ax1.tick_params(axis='y', labelcolor=color1_0)
        ax1.set_ylabel(y1label, color=color1_0)
        ax2.tick_params(axis='y', labelcolor=color2_0)
        ax2.set_ylabel(y2label, color=color2_0)


    if IncF.b_TS_hlmodel:
        highlight_first_line(ax2, labels2, IncF.c_TS_hlmodel, case="Model")
    if IncF.b_TS_hlobs:
        #print(labels1)
        highlight_first_line(ax1, labels1, IncF.c_TS_hlobs, case="Obs")

    do_legend(fig, ax1, handles, new_legend, legend_title)

    
    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)
    plt.close('all')
    return None


###################################
## Plot Two Sides var
def plot_rolling_all(name, data_or, color_dict, timeframe, var1, var2, units, direc="", legend_title="", window=1, highlight=""):
    ## Make sure data is for timeframe asked
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data_or.loc[date1:date2]

    ## Check wether it's the same variable that is being plotted
    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    equiv1=[k for  k, v in variables_equivalences.items() if var1 in v][0]
    equiv2=[k for  k, v in variables_equivalences.items() if var1 in v][0]

    if equiv1==equiv2:
        name="%s___%s_Raw+Rolling%s___%s-%s"%(name, var1, window, date1, date2)
        ## Prepare data
        case, color1_0, color2_0,  data_1, data_2 = prepare_data_one_var(data_or, color_dict)
        BoolTwoVars=False
        if case is None:
            return None 
    else:
        BoolTwoVars=True
        name="%s___%s-%s_Raw+Rolling%s___%s-%s"%(name, var1, var2, window, date1, date2)
        ## Separate datacolor_dict
        data_1=data[[v for v in data.columns if var1 in v ]] 
        data_2=data[[v for v in data.columns if var2 in v ]] 
        color1_0=color_dict[data_1.columns[0]]
        color2_0=color_dict[data_2.columns[0]]



    data_1R=data_1.rolling(window=window,  center=bool_center, min_periods=int(min_perc_window*window/100)).mean().dropna()
    data_2R=data_2.rolling(window=window,  center=bool_center, min_periods=int(min_perc_window*window/100)).mean().dropna()

    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  
    

    y1label, y2label, new_legend = prepare_labels(data.columns, data_1.columns, data_2.columns, var1, var2, units)



    prepare_metadata(fig_meta, ax_meta, legend_title, var1, var2)



    ## Plot
    handles1, labels1 = [], []
    if not data_1.empty:
        data_1R=data_1.rolling(window=window,  center=bool_center, min_periods=int(min_perc_window*window/100)).mean()
        data_1R.plot(color=[color_dict[c] for c in data_1R.columns], ax=ax1, style=ls, ms=ms)
        data_1.plot(color=[color_dict[c] for c in data_1.columns], ax=ax1, ms=1.5*ms, alpha=0.3, style=".")
        handles1, labels1 = ax1.get_legend_handles_labels()
        ax1.get_legend().remove()
    handles2, labels2 = [], []
    if not data_2.empty:
        data_2R=data_2.rolling(window=window,  center=bool_center, min_periods=int(min_perc_window*window/100)).mean()
        data_2R.plot(color=[color_dict[c] for c in data_2R.columns], ax=ax2, style=ls, ms=ms)
        data_2.plot(color=[color_dict[c] for c in data_2.columns], ax=ax2, ms=1.5*ms, alpha=0.3, style=".")
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.get_legend().remove()

    

    ## Both vars werent plotted
    if len(handles2)==0 and len(handles1)==0:
        fig.clf()
        plt.close(fig)
        fig_meta.clf()
        plt.close(fig_meta)
        plt.close('all')
        return None

    ## Save image because everything was plotted
    handles = handles1[:int(len(handles1)/2)] + handles2[:int(len(handles2)/2)]
    labels = labels1[:int(len(labels1)/2)]  + labels2[:int(len(labels2)/2)]
    align_xaxis(ax1)
    BoolOneAxis = align_yaxis(ax1, ax2, BoolTwoVars=BoolTwoVars, var1=var1, var2=var2)
    if BoolOneAxis:
        ax1.set_ylabel(y1label)
    else:
        ax1.tick_params(axis='y', labelcolor=color1_0)
        ax1.set_ylabel(y1label, color=color1_0)
        ax2.tick_params(axis='y', labelcolor=color2_0)
        ax2.set_ylabel(y2label, color=color2_0)


    if IncF.b_TS_hlmodel:
        highlight_first_line(ax2, labels2, IncF.c_TS_hlmodel, case="Model")
    if IncF.b_TS_hlobs:
        #print(labels1)
        highlight_first_line(ax1, labels1, IncF.c_TS_hlobs, case="Obs")

        
    do_legend(fig, ax1, handles, new_legend, legend_title)

    
    ## Save Figure
    aux_plots.savefigure(fig_meta, join(direc, "metadata"), name+"_metadata", BoolMessage=False)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig_meta.clf()
    plt.close(fig_meta)
    fig.clf()
    plt.close(fig)
    plt.close('all')
    return None



#####################################
## Make Plot. It redirects to the appropriate one depending on what is in vars_to_plot
def make_plots(name, data, vars_to_plot, units, direc="", plot_freq="", case="", legend_title="", window=1, variables=None, highlight=""):

    makedirs(direc, exist_ok = True) 
    makedirs(join(direc, "metadata"), exist_ok = True) 
    ## All variables asked
    if variables is None:
        variables = sorted(list(set(aux.flat(vars_to_plot, out=[]))))


    ## Colors for each line/column
    color_dict=get_color_columns(data.columns)

    #####################
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

        date1=group_data.index[0].strftime("%d-%m-%Y")  
        date2=group_data.index[-1].strftime("%d-%m-%Y")
        timeframe=[date1, date2]


        ## Ommit plot if missing information of one variable
        cols=group_data.columns
        BoolData=True
        for var in variables:
            if group_data[[v for v in cols if var in v]].dropna(how="all").empty:
                BoolData=False
        if b_showshared:
            group_data=group_data.dropna(thresh=2)
            if group_data.empty:
                BoolData=False            

        if not BoolData:
            #print(f"   ...Skipping  {name}. Missing variables.")
            break
        try:
            var1, var2 = variables
        except:
            var1 = var2 = variables[0]

        if 1 in TS_plots:
            plot_basic(name, group_data, color_dict, timeframe, var1, var2, units, direc=direc, legend_title=legend_title, highlight=highlight)
        if 2 in TS_plots:
            plot_rolling(name, group_data, color_dict, timeframe, var1, var2, units, direc=direc, legend_title=legend_title, window=window, highlight=highlight)
        if 3 in TS_plots:
            plot_rolling_all(name, group_data, color_dict, timeframe, var1, var2, units, direc=direc, legend_title=legend_title, window=window, highlight=highlight)


###############################################################
###############################################################
#####################################
## TODO: Fix code to include timeseries of 3d models
def main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData):
    if IncF.i_TS==0:
        print("...Plotting was not requested.")
        return None
    if t_Stations_Info.empty:
        print("...No stations were read. Skipping this step.")
        #return None     

  
    ###########
    ## Start going over specifications for plots 
    path_base=join(IncF.c_FigureDir, "TimeSeries")
    makedirs(path_base, exist_ok = True) 
    for count, plot_freq in enumerate(IncF.c_TS_freq):
        freq, plot_freq =plot_freq.split(",")
        if 2 in TS_plots or 3 in TS_plots:
            window=IncF.i_TS_RM[count]
        else:
            window=1

        path_freq=join(path_base, f"{freq},{plot_freq}")
        makedirs(path_freq, exist_ok = True) 

        ###########
        for vars_plot in IncF.c_TS_Plots:
            vars_to_plot=aux.find_equiv_vars(vars_plot)
            
            print("...Plotting:  ", vars_plot, "  for freq   ", freq, " and time period:  ", plot_freq)

            if len(vars_to_plot)==0:
                print("   ...Skipping. Network requested did not have data.")
                continue
            variables=set(aux.flat(vars_to_plot, out=[]))


            #################
            data, data_obs, data_mod=get_new_data(vars_to_plot, t_ObsStationData[freq], t_ModelStationData[freq], plot_freq)
            data.dropna(how="all", inplace=True)
            data_obs.dropna(how="all", inplace=True)
            data_mod.dropna(how="all", inplace=True)

            units=data.meta.units
            try:
                units["PMC"]=units["DUST"]
            except:
                pass

            networks=sorted(list(set([c[0].split()[1] for c in data.columns if "Obs" in c[0]])))
            models=sorted(list(set([c[0].split()[1] for c in data.columns if "Model" in c[0]])))
            mods="-".join(models)
            

            ## TODO: CHECK THAT THEY ARE THE SAME
            try:
                networks=list(vars_to_plot["O"].keys())
            except:
                networks=[]
            try:
                models=list(vars_to_plot["M"].keys())
                mod_1=models[0]
            except:
                models=[]
                mod_1=""

            ####################
            ## If there are no common dates this will happen
            if len(networks)==0:
                if len(models)==0:
                    continue

                print("   ...Plotting only models for averages found in area filters.")
                for station in minif["ID"]:
                    try:
                        data_station=data.xs(station, axis=1, level=1)
                    except:
                        continue

                    ############
                    name="%s_%s"%(station, mods)
                    if mods=="":
                        name=name[:-1]
                    legend_title="%s:  %s"%(station, mods)
                    locations=[legend_title]
                    make_plots(name, data_station, vars_to_plot,  units, direc=join(path_freq, "Models", station), plot_freq=plot_freq, legend_title=locations, window=window, highlight=mod_1)


                continue           


            if data.empty:
                print("   ...Skipping this step. There is no data to plot for at least one variable. Check timeframe and variables asked to plot.")
                continue
            if len(variables)>2:
                print("   ...Skipping this step. Plot with more than two variables is not implemented.")
                continue
            if len(networks)>2:
                print("   ...Skipping this step. Plot with more than two networks cant be processed yet. The maximum allowed is two. ")
                continue

            ####################
            ## Plot
            direc=join(path_freq, "-".join(["_".join(n.split("_")[:-1]) for n in networks]))
            makedirs(direc, exist_ok = True) 
            makedirs(join(direc, "metadata"), exist_ok = True) 

            if len(networks)==1:
                ##
                network=list(vars_to_plot["O"].keys())[0]
                
                info=aux.get_info_group(t_Stations_Info, network)
                minif=info.copy()
                minif=minif[minif["ID"].str.contains('^f[0-9]-+')]
                mini=info[~info["ID"].isin(minif["ID"])]
                filters=t_Stations_Filters[network]
                mods="-".join(models)
                

                if len(models)!=0:
                    info_mods=pd.concat([aux.get_info_group(t_ModelStation_Info[m], network) for m in models], axis=0, keys=models).reset_index()
                else:
                    info_mods=pd.DataFrame(columns=["ID", "location", "Nombre", "Distance"])

                ##############
                ## Individual plots
                if bool_plotstations:
                    for station in mini["ID"]:
                        try:
                            data_station=data.xs(station, axis=1, level=1)
                        except:
                            continue
                            
                        ID_mod=info_mods[info_mods["ID"]==station]
                        ID_sta=info[info["ID"]==station]
                 
                        loc_obs="\n".join(ID_sta.apply(lambda x:  "(%s)  %s "%(x.location, x.Nombre), axis=1).values)
                        loc_mod="\n".join(ID_mod.apply(lambda x:  "(%s)  %s    Distance: %.3f km"%(x.location, x.level_0, x.Distance), axis=1).values)

                        ############
                        name="%s--%s_%s"%(station, "-".join(info[info["ID"]==station]["Nombre"].values[0].split()).replace(".", "").replace(",", ""), mods)
                        if mods=="":
                            name=name[:-1]
                        legend_title="%s:  %s"%(networks[0], info[info["ID"]==station]["Nombre"].values[0])
                        locations=[legend_title, loc_obs, loc_mod]
                        make_plots(name, data_station, vars_to_plot,  units, direc=direc, plot_freq=plot_freq, legend_title=locations, window=window, highlight=mod_1)


                ##############
                ## filter plots
                if bool_plotfilters:
                    for station in minif["ID"]:
                        try:
                            data_station=data.xs(station, axis=1, level=1)
                        except:
                            continue

                        ############
                        name="%s_%s"%(station, mods)
                        if mods=="":
                            name=name[:-1]
                        legend_title="%s:  %s"%(networks[0], station)
                        locations=[legend_title]
                        make_plots(name, data_station, vars_to_plot,  units, direc=join(direc, station), plot_freq=plot_freq, legend_title=locations, window=window, highlight=mod_1)


            if len(networks)==2:
                network1=list(vars_to_plot["O"].keys())[0]
                network2=list(vars_to_plot["O"].keys())[1]
                
                info1=aux.get_info_group(t_Stations_Info, network1)
                info2=aux.get_info_group(t_Stations_Info, network2)

                closest1=find_closest_points(info1, info2)
                closest2=find_closest_points(info2, info1)
                

                ##############
                ## Individual plots fot network 1
                for station1 in info1["ID"]:
                    data_station1=data[[c for c in data.columns if  networks[0] in c[0] and station1 in c[1]]].xs(station1, axis=1, level=1)
                    for station2 in info2["ID"]:
                        ## If closest station only skip the rest
                        if b_closest:
                            if station2!=closest1[closest1["ID"]==station1]["ID_Nearest"].values[0]:
                                continue
 
                        ########
                        data_station2=data[[c for c in data.columns if networks[1] in c[0] and station2 in c[1]]].xs(station2, axis=1, level=1)
                        ##
                        data_station=pd.concat([data_station1, data_station2], axis=1)

                        ############
                        name="%s--%s_%s--%s"%(station1, info1[info1["ID"]==station1]["Nombre"].values[0], station2, info2[info2["ID"]==station2]["Nombre"].values[0])
                        legend_title=name
                        make_plots(name, data_station, vars_to_plot, units, direc=direc, plot_freq=plot_freq, legend_title=legend_title, window=window)


                ##############
                ## Individual plots fot network 2
                for station2 in info2["ID"]:
                    data_station2=data[[c for c in data.columns if  networks[1] in c[0] and  station2 in c[1]]].xs(station2, axis=1, level=1)
                    for station1 in info1["ID"]:
                        ## If closest station only skip the rest
                        if b_closest: 
                            if station1!=closest2[closest2["ID"]==station2]["ID_Nearest"].values[0]:
                                continue
 
                        ######## 
                        data_station1=data[[c for c in data.columns if networks[0] in c[0] and station1 in c[1]]].xs(station1, axis=1, level=1)
                        ##
                        data_station=pd.concat([data_station1, data_station2], axis=1)

                        ############
                        name="%s--%s_%s--%s"%(station1, info1[info1["ID"]==station1]["Nombre"].values[0], station2, info2[info2["ID"]==station2]["Nombre"].values[0])
                        legend_title=name
                        make_plots(name, data_station, vars_to_plot,  units, direc=direc, plot_freq=plot_freq, legend_title=legend_title, window=window)

