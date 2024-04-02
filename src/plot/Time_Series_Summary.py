# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import numpy as np
import pandas as pd
#
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#
import src.common as aux
import src.plot.aux_plots as aux_plots
#
from os.path import join
from os import makedirs
#
import geopy.distance
from shapely.ops import nearest_points
from copy import deepcopy
from src.plot.aux_timeseries import *
#plt.rcParams['fig.alpha'] = 0
plt.rcParams['axes.facecolor'] = 'white'
#plt.rcParams['axes.alpha'] = 1
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.edgecolor'] = 'none'


TS_plots=[int(p) for p in str(IncF.i_TS_Sum)]

######################################
settings=aux_plots.TSConfig
##
color1, color2=settings.colors
figsize = settings.imagesize
ms = settings.markersize
ls = settings.linestyle
legend_location = settings.legend_location
titlesize = settings.titlesize
ticksize = settings.ticksize
labelsize = settings.labelsize
overlap_percentage = settings.overlap_percentage
overlap_twovar = 0.4
legendsize = settings.legendsize
b_closest = settings.closest
min_perc_window = settings.percwindow
bool_center = settings.center
b_showshared = False #settings.b_showshared
subtitlesize= settings.subtitlesize
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
    


def prepare_metadata(fig_meta, ax_meta, locations, clean_var_x, clean_var_y):
    ## Figure Metadata
    fig_meta=fig_meta or plt.subplots(1, figsize=figsize)
    ax_meta=ax_meta or plt.gca()

    ax_meta.set_axis_off()
    
    fig_meta.suptitle('%s vs %s'%(clean_var_x, clean_var_y), fontsize=titlesize, y=1)
    ## Titles and subtitles
    subtitle='\n'.join(map(str, locations[1:]))
    ax_meta.set_title(subtitle, fontsize=subtitlesize, y=0.975, loc="left", fontdict={'verticalalignment':"top"})



#####################################
def prepare_labels(columns, columns_1, columns_2, var1, var2, units):
    ##
    columns=[c[0] for c in columns][::2]
    columns_1=[c[0] for c in columns_1][::2]
    columns_2=[c[0] for c in columns_2][::2]
    model=list(set([c.split()[1] for c in columns if "Model" in c]))
    network=list(set([c.split()[1] for c in columns if "Obs" in c]))
    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    evar1=[k for  k, v in variables_equivalences.items() if var1 in v] 
    evar1=evar1[0] if evar1 else var1
    evar2=[k for  k, v in variables_equivalences.items() if var1 in v] 
    evar2=evar2[0] if evar2 else var2
    
    ## Dont repeat station because it will always be in reference to it
    if len(network)==1:
        if len(model)!=0:
            ##
            new_legend=[c.replace("Obs ", "").replace("Model ", "") for c in columns_1]+ [c.replace("Obs ", "").replace("Model ", "") for c in columns_2]
            ## Removes var from legend
            new_legend=[" ".join(c.split(" ")[:-1]) for c in new_legend]
        else:
            ## If no models then only show ID of station
            new_legend=[c.split()[-1] for c in columns_1]+ [c.split()[-1] for c in columns_2]

        
        y1label="%s %s"%(evar1, aux_plots.get_pretty_units(units[var1]))
        y2label="%s %s"%(evar2, aux_plots.get_pretty_units(units[var2]))

        return y1label, y2label, new_legend

    elif len(network)==2:
        ##
        new_legend=[c.replace("Obs ", "").replace("Model ", "") for c in columns_1]+ [c.replace("Obs ", "").replace("Model ", "") for c in columns_2]
        ## Removes var from legend
        new_legend=[" ".join(c.split(" ")[:-1]) for c in new_legend]

        units_1=units[var1]#aux.read_units(columns_1[0].split()[1])
        units_2=units[var2]#aux.read_units(columns_2[0].split()[1])
        y1label="%s %s"%(evar1, aux_plots.get_pretty_units(units[var1]))
        y2label="%s %s"%(evar2, aux_plots.get_pretty_units(units[var2]))
        return y1label, y2label, new_legend 

    elif len(network)==0:
        new_legend=[c.replace("Obs ", "").replace("Model ", "") for c in columns_1]+ [c.replace("Obs ", "").replace("Model ", "") for c in columns_2]
        new_legend=[c.replace(var1, "").replace(var2, "") for c in new_legend]

        units=aux.read_units("SINCA")
        y1label="%s %s"%(evar1, aux_plots.get_pretty_units(units[var1]))
        y2label="%s %s"%(evar2, aux_plots.get_pretty_units(units[var2]))
        return y1label, y2label, new_legend     

#####################################
def get_new_data(vars_to_plot, t_ObsData, t_ModData, plot_freq, plot_type):

    BoolObs=False
    BoolMod=False
    ## Merge obs data 
    try:
        obs_data={f"Obs {aux.clean_network(group)} {var}":t_ObsData[group][var] for group in vars_to_plot["O"] for var in vars_to_plot["O"][group]}
        obs_units={var:t_ObsData[group].meta.units[var] for group in vars_to_plot["O"] for var in vars_to_plot["O"][group]}
        #obs_data=pd.concat([t_ObsData[list(vars_to_plot["O"].keys())[0]]["groupby"]]+ list(obs_data.values()), keys=["groupby"]+list(obs_data.keys()), axis=1 )
        obs_data=pd.concat(list(obs_data.values()), keys=list(obs_data.keys()), axis=1)
        obs_data.meta.units=obs_units
        #obs_data.dropna(how="all", axis=1, inplace=True)
        BoolObs=True
    except:
        obs_data=pd.DataFrame()
        obs_units={}
    

    ## Merge model data 
    try:
        mod_data={f"Model {model} {var}":t_ModData[model][group][var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model] for group in vars_to_plot["O"]}
        mod_units={var:t_ModData[model][group].meta.units[var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model] for group in vars_to_plot["O"]}
        mod_data=pd.concat(list(mod_data.values()), keys=list(mod_data.keys()), axis=1)
        mod_data.meta.units=mod_units
        #mod_data.dropna(how="all", axis=1, inplace=True)
        BoolMod=True
    except:
        mod_data=pd.DataFrame()
        mod_units={}

    
    if BoolObs and BoolMod:
        if IncF.i_TS_Sum_SharedXRange==1:
            idx_common = obs_data.index.intersection(mod_data.index)
            mod_data.loc[mod_data.index.difference(idx_common)]=np.nan


        ## Units
        mod_data.meta.units=mod_units
        obs_data.meta.units=obs_units


    ##########
    data=pd.concat([obs_data, mod_data], axis=1)#.dropna()
    ## Get info of cycles
    if plot_type=="DC":
        data_groupby = data.index.map(lambda t: t.replace(year=2000, month=1, day=1)) 
    elif plot_type=="AC":
        data_groupby = data.index.map(lambda t: t.replace(year=2000)) 
    
    try:
        stations = list(set(data.columns.get_level_values(1)))
    except:
        stations = []
    for station in stations:
        data["groupby", station] = data_groupby

    data.meta.units={**obs_units, **mod_units}#{"O":obs_units, "M":mod_units}
    return data, obs_data, mod_data
    

###############################################################
###############################################################
#####################################
## Plot One var
def plot_one_var_basic(name, data_or, color_dict, timeframe, var, units, direc="", legend_title="", plot_type="DC"):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])

    name="%s___%s___%s_%s-%s"%(name, var, plot_type, date1, date2)

    ## Prepare data
    case, color1_0, color2_0,  data_1, data_2 = prepare_data_one_var(data_or, color_dict)
    if case is None:
        return None

    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  
    
    prepare_metadata(fig_meta, ax_meta, legend_title, var, var)

    ## Settings
    y1label, y2label, new_legend = prepare_labels(data_or.columns, data_1.columns, data_2.columns, var, var, units)

    ## Plot
    handles1, labels1=[], []
    if not data_1.empty:
        data_1=data_1.xs("mean", axis=1, level=1)
        data_1.plot(color=[color_dict[c] for c in data_1.columns], ax=ax1, style=ls, ms=ms)
        handles1, labels1 = ax1.get_legend_handles_labels()
        ax1.get_legend().remove()
    handles2, labels2=[], []
    if not data_2.empty:
        data_2=data_2.xs("mean", axis=1, level=1)
        data_2.plot(color=[color_dict[c] for c in data_2.columns], ax=ax2, style=ls, ms=ms)
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.get_legend().remove()


    ################
    handles= handles1 + handles2
    labels = labels1  + labels2
    do_legend(fig, ax1, handles, new_legend, legend_title, plot_type)

    #########
    ## Axis
    BoolOneAxis=align_yaxis(ax1, ax2, BoolTwoVars=False, var1=var, var2=var)
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

    ########
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
def plot_two_var_basic(name, data, color_dict, timeframe, var1, var2, units, direc="", legend_title="", plot_type="DC"):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    var="%s-%s"%(var1, var2)
    name="%s___%s___%s_%s-%s"%(name, var, plot_type, date1, date2)

    ## Separate data

    print(data.head())
    print(data.columns)

    print(var1, var2)
    data_1=data[[v for v in data.columns if var1==v.split()[-1]]] 
    data_2=data[[v for v in data.columns if var2==v.split()[-1]]] 
    color1_0=color_dict[data_1.columns[0]]
    color2_0=color_dict[data_2.columns[0]]

    y1label, y2label, new_legend = prepare_labels(data.columns, data_1.columns, data_2.columns, var1, var2, units)


    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  
    
    prepare_metadata(fig_meta, ax_meta, legend_title, var1, var2)

    ## Settings
    ax1.tick_params(axis='y', labelcolor=color1_0)
    ax1.set_ylabel(y1label, color=color1_0)
    ax2.tick_params(axis='y', labelcolor=color2_0)
    ax2.set_ylabel(y2label, color=color2_0)


    ## Plot
    data_1.plot(color=[color_dict[c] for c in data_1.columns], ax=ax1, style=ls, ms=ms)
    data_2.plot(color=[color_dict[c] for c in data_2.columns], ax=ax2, style=ls, ms=ms)

    ## Legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    if IncF.b_TS_hlmodel:
        if len(labels2)!=0:
            highlight_first_line(ax2, labels2, IncF.c_TS_hlmodel, case="Model")
    if IncF.b_TS_hlobs:
        if len(labels1)!=0:
            highlight_first_line(ax1, labels1, IncF.c_TS_hlobs, case="Obs")

    ## One var wasnt plotted
    if len(handles2)==0 or len(handles1)==0:
        fig.clf()
        plt.close(fig)
        fig_meta.clf()
        plt.close(fig_meta)
        plt.close('all')
        return None

    ## Save image because everything was plotted
    else:
        handles= handles1 + handles2
        labels = labels1  + labels2
        ax1.get_legend().remove()
        ax2.get_legend().remove()
        align_yaxis(ax1, ax2, BoolTwoVars=True, var1=var1, var2=var2)
        do_legend(fig, ax1, handles, new_legend, legend_title, plot_type)
    
        
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
## Plot One var
def plot_one_var_std(name, data_or, color_dict, timeframe, var, units, direc="", legend_title="", plot_type="DC"):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])

    name="%s___%s___%s-std_%s-%s"%(name, var, plot_type, date1, date2)

    ## Prepare data
    case, color1_0, color2_0,  data_1, data_2 = prepare_data_one_var(data_or, color_dict)
    if case is None:
        return None

    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  
    
    prepare_metadata(fig_meta, ax_meta, legend_title, var, var)

    ## Settings
    y1label, y2label, new_legend = prepare_labels(data_or.columns, data_1.columns, data_2.columns, var, var, units)

    ## Plot
    handles1, labels1=[], []
    if not data_1.empty:
        data_mean=data_1.xs("mean", axis=1, level=1)
        data_std=data_1.xs("std", axis=1, level=1)
        up=data_mean+data_std
        down=data_mean-data_std
        len1=len(data_mean.columns)
        data_mean.plot(color=[color_dict[c] for c in data_mean.columns], ax=ax1, style=ls, ms=ms, zorder=100)

        for c in data_mean.columns:
            new_color=aux_plots.lighten_color(color_dict[c], amount=0.5)
            up[c].plot(color=new_color, ax=ax1, style="-", alpha=0.3, zorder=1)
            down[c].plot(color=new_color, ax=ax1, style="-", alpha=0.3, zorder=1)
            ax1.fill_between(data_mean.index, up[c], down[c], color=new_color, alpha=0.3)

        handles1, labels1 = ax1.get_legend_handles_labels()
        ax1.get_legend().remove()
    handles2, labels2=[], []
    if not data_2.empty:
        data_mean=data_2.xs("mean", axis=1, level=1)
        data_std=data_2.xs("std", axis=1, level=1)
        up=data_mean+data_std
        down=data_mean-data_std
        len2=len(data_mean.columns)
        data_mean.plot(color=[color_dict[c] for c in data_mean.columns], ax=ax2, style=ls, ms=ms, zorder=100)

        for c in data_mean.columns:
            new_color=aux_plots.lighten_color(color_dict[c], amount=0.5)
            up[c].plot(color=new_color, ax=ax2, style="-", alpha=0.3, zorder=1)
            down[c].plot(color=new_color, ax=ax2, style="-", alpha=0.3, zorder=1)
            ax2.fill_between(data_mean.index, up[c], down[c], color=new_color, alpha=0.3)

        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.get_legend().remove()

    
    ################
    handles = handles1[:len1] + handles2[:len2]
    labels = labels1[:len1] + labels2[:len2]  
    #handles = handles1 + handles2
    #labels = labels1 + labels2 
    
    do_legend(fig, ax1, handles, new_legend, legend_title, plot_type)


    #########
    ## Axis
    BoolOneAxis=align_yaxis(ax1, ax2, BoolTwoVars=False, var1=var, var2=var)
    if BoolOneAxis:
        ax1.set_ylabel(y1label)
    else:
        ax1.tick_params(axis='y', labelcolor=color1_0)
        ax1.set_ylabel(y1label, color=color1_0)
        ax2.tick_params(axis='y', labelcolor=color2_0)
        ax2.set_ylabel(y2label, color=color2_0)


    ########
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
def plot_two_var_std(name, data_or, color_dict, timeframe, var1, var2, units, direc="", legend_title="", plot_type="DC"):
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    var="%s-%s"%(var1, var2)
    name="%s___%s___%s_%s-%s"%(name, var, plot_type, date1, date2)

    ## Separate data
    data_1=data[[v for v in data.columns if var1==v.split()[-1]]] 
    data_2=data[[v for v in data.columns if var2==v.split()[-1]]] 
    color1_0=color_dict[data_1.columns[0]]
    color2_0=color_dict[data_2.columns[0]]

    y1label, y2label, new_legend = prepare_labels(data.columns, data_1.columns, data_2.columns, var1, var2, units)


    ###############
    ### Plot
    fig_meta, ax_meta = plt.subplots(1) 
    fig, ax1 = plt.subplots(1, figsize=figsize)
    ax2 = ax1.twinx()  
    
    prepare_metadata(fig_meta, ax_meta, legend_title, var1, var2)

    ## Settings
    ax1.tick_params(axis='y', labelcolor=color1_0)
    ax1.set_ylabel(y1label, color=color1_0)
    ax2.tick_params(axis='y', labelcolor=color2_0)
    ax2.set_ylabel(y2label, color=color2_0)


    ## Plot
    data_1.plot(color=[color_dict[c] for c in data_1.columns], ax=ax1, style=ls, ms=ms)
    data_2.plot(color=[color_dict[c] for c in data_2.columns], ax=ax2, style=ls, ms=ms)

    ## Legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    ## One var wasnt plotted
    if len(handles2)==0 or len(handles1)==0:
        fig.clf()
        plt.close(fig)
        fig_meta.clf()
        plt.close(fig_meta)
        plt.close('all')
        return None

    ## Save image because everything was plotted
    else:
        handles= handles1 + handles2
        labels = labels1  + labels2
        ax1.get_legend().remove()
        ax2.get_legend().remove()
        align_yaxis(ax1, ax2, BoolTwoVars=True, var1=var1, var2=var2)
        do_legend(fig, ax1, handles, new_legend, legend_title, plot_type)
    
        
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
def make_plots(name, data, vars_to_plot, units, direc="", plot_freq="", case="", legend_title="", plot_type="DC"):

    ## All variables asked
    variables=sorted(list(set(aux.flat(vars_to_plot, out=[]))))

    ## Colors for each line/column
    color_dict=get_color_columns(data.columns)


    
    #####################
    ## Plots with the frequency asked
    if plot_freq=="Whole":
        data["plot_freq"]=1
        grouped_data=data.groupby("plot_freq")
    else:
        grouped_data=data.groupby(pd.Grouper(freq=plot_freq))
              

    #########
    list_dates=list(grouped_data.groups.keys()) 
    for group in list_dates:
        try:
            group_data=grouped_data.get_group(group)
            date1=group_data.index[0].strftime("%d-%m-%Y")  
            date2=group_data.index[-1].strftime("%d-%m-%Y")
            group_data=group_data.groupby("groupby").agg(["mean", "std"]).dropna(how="all")
        except:
            continue


        group_data=group_data[[c for c in group_data.columns if c[0] !="plot_freq"]]

        timeframe=[date1, date2]


        ## Ommit plot if missing information of one variable
        cols=group_data.columns
        BoolData=True
        for var in variables:
            if group_data[[v for v in cols if var in v[0]]].dropna(how="all").empty:
                BoolData=False
        if b_showshared:
            group_data=group_data.dropna(thresh=2)
            if group_data.empty:
                BoolData=False            

        if not BoolData:
            #print(f"   ...Skipping  {name}. Missing variables.")
            break

        ###########
        
        if len(variables)==2:
            var1, var2 = variables
            if 1 in TS_plots:
                plot_two_var_basic(name, group_data, color_dict, timeframe, var1, var2, units, direc=direc, legend_title=legend_title, plot_type=plot_type)
            if 2 in TS_plots:
                plot_two_var_std(name, group_data, color_dict, timeframe, var1, var2, units, direc=direc, legend_title=legend_title, plot_type=plot_type)


        elif len(variables)==1:
            if 1 in TS_plots:
                plot_one_var_basic(name, group_data, color_dict, timeframe, variables[0], units, direc=direc, legend_title=legend_title, plot_type=plot_type) 
            if 2 in TS_plots:
                plot_one_var_std(name, group_data, color_dict, timeframe, variables[0], units, direc=direc, legend_title=legend_title, plot_type=plot_type) 


def get_model(df, freq, group):
    data = df.reset_index(level=["time"])
    data=data.loc[group].groupby("time").mean()
    var=data.columns[0] 
    return data.resample(freq).mean()[var]


#####################################
## TODO: Fix code to include timeseries of 3d models
def main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData):
    if IncF.i_TS_Sum==0:
        print("...Plotting was not requested.")
        return None
    if t_Stations_Info.empty:
        print("...No stations were read. Skipping this step.")
        return None     

  

    ###########
    ## Start going over specifications for plots 
    path_base=join(IncF.c_FigureDir, "TimeSeries")
    makedirs(path_base, exist_ok = True) 
    for count, plot_freq in enumerate(IncF.c_TS_Sum_freq):
        freq, plot_type, plot_freq =plot_freq.split(",")

        path_freq=join(path_base, f"{plot_type},{plot_freq}")
        makedirs(path_freq, exist_ok = True) 

        ###########
        for vars_plot in IncF.c_TS_Sum_Plots:
            vars_to_plot=aux.find_equiv_vars(vars_plot)
            print("...Plotting:  ", vars_plot, "  for plot type   ", plot_type, " and time period:  ", plot_freq)

            if len(vars_to_plot)==0:
                print("   ...Skipping. Network requested did not have data.")
                continue
            variables=set(aux.flat(vars_to_plot, out=[]))


            ###################
            data, data_obs, data_mod=get_new_data(vars_to_plot, t_ObsStationData[freq], t_ModelStationData[freq], plot_freq, plot_type=plot_type)
            units=data.meta.units

            networks=sorted(list(set([c[0].split()[1] for c in data.columns if "Obs" in c[0]])))
            models=sorted(list(set([c[0].split()[1] for c in data.columns if "Model" in c[0]])))

            ## TODO: CHECK THAT THEY ARE THE SAME
            try:
                networks=list(vars_to_plot["O"].keys())
            except:
                networks=[]
            try:
                models=list(vars_to_plot["M"].keys())
            except:
                models=[]

            ####################
            ## If there are no common dates this will happen
            if len(networks)==0:
                if len(models)==0:
                    continue

                print("   ...Plotting only models for averages found in area filters.")
                ###############
                ## Make plots related to area filters that would be shared between different networks
               
                filters=t_Model_Filters
                mod_direc=join(path_freq, "Models")
  
                
                for f in filters[models[0]]:
                    f_groups=filters[models[0]][f] 
                    filter_data=[]
                    for index, group in enumerate(f_groups):
                        mod_data={f"Model {model} {var}":get_model(t_ModelData[freq][model][var], freq, group) for model in vars_to_plot["M"] for var in vars_to_plot["M"][model]}
                        mod_units={var:t_ModelData[freq][model].meta.units[var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model]}
                        mod_data=pd.concat(list(mod_data.values()), keys=list(mod_data.keys()), axis=1)
                        mod_data.meta.units=mod_units
                        filter_data.append(mod_data)

 
                        direc=join(mod_direc, f"{f}-{index}")
                        makedirs(direc, exist_ok = True) 
                        makedirs(join(direc, "metadata"), exist_ok = True) 
  
                        ############
                        ## Mean for observations and models
                        make_plots(f"{f}-{index}__{'-'.join(models)}", mod_data, vars_to_plot,  mod_units, direc=direc, plot_freq=plot_freq, legend_title=f"Averages {f}-{index}")
                    if f=="f1":
                        data=pd.concat(filter_data, axis=1, keys=range(len(filter_data)))
                        data_mean=data.groupby(pd.Grouper(freq="D")).mean()
                        data_max=data.groupby(pd.Grouper(freq="D")).max()
                        for i in range(2):
                            make_plots(f"{f}-{i}_mean__{'-'.join(models)}", data_mean[i], vars_to_plot,  mod_units, direc=direc, plot_freq=plot_freq, legend_title=f"Averages {f}-{index}")
                            make_plots(f"{f}-{i}_max__{'-'.join(models)}", data_max[i], vars_to_plot,  mod_units, direc=direc, plot_freq=plot_freq, legend_title=f"Averages {f}-{index}")



                        data_mean=data_mean[0]/data_mean[1]
                        data_max=data_max[0]/data_max[1]

                        make_plots(f"{f}-0div1_mean__{'-'.join(models)}", data_mean, vars_to_plot,  mod_units, direc=direc, plot_freq=plot_freq, legend_title=f"Averages {f}-{index}")
                        make_plots(f"{f}-0div1_max__{'-'.join(models)}", data_max, vars_to_plot,  mod_units, direc=direc, plot_freq=plot_freq, legend_title=f"Averages {f}-{index}")
                
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
            direc=join(path_freq, "-".join(networks))
            makedirs(direc, exist_ok = True) 
            makedirs(join(direc, "metadata"), exist_ok = True) 
            
            if len(networks)==1:
                ##
                network=list(vars_to_plot["O"].keys())[0]
                
                info=aux.get_info_group(t_Stations_Info, network)
                filters=t_Stations_Filters[network]
                mods="-".join(models)

                if len(models)!=0:
                    info_mods=pd.concat([aux.get_info_group(t_ModelStation_Info[m], network) for m in models], axis=0, keys=models).reset_index()
                else:
                    info_mods=pd.DataFrame(columns=["ID", "location", "Nombre", "Distance"])

                ##############
                ## Individual plots
                for station in info["ID"]:
                    try:
                        data_station=data.xs(station, axis=1, level=1)
                    except:
                        continue
                        
                    ID_mod=info_mods[info_mods["ID"]==station]
                    ID_sta=info[info["ID"]==station]
             
                    loc_obs="\n".join(ID_sta.apply(lambda x:  "(%s)  %s "%(x.location, x.Nombre), axis=1).values)
                    loc_mod="\n".join(ID_mod.apply(lambda x:  "(%s)  %s    Distance: %.3f km"%(x.location, x.level_0, x.Distance), axis=1).values)

                    ############
                    name="%s--%s_%s"%(station, "-".join(str(info[info["ID"]==station]["Nombre"].values[0]).split()).replace(".", "").replace(",", ""), mods)
                    legend_title="%s:  %s"%(networks[0], info[info["ID"]==station]["Nombre"].values[0])
                    locations=[legend_title, loc_obs, loc_mod]

                    make_plots(name, data_station, vars_to_plot,  units, direc=direc, plot_freq=plot_freq, legend_title=locations, plot_type=plot_type)

