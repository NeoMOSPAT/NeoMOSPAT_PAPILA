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



TS_plots=[int(p) for p in str(IncF.i_TS_stats)]
bool_plotstations=IncF.b_TS_stats_plotstations
bool_plotfilters=IncF.b_TS_stats_plotfilters
######################################
settings=aux_plots.TSConfig
##
color1, color2=settings.colors
figsize = settings.imagesize
ms = settings.markersize
ls = settings.linestyle
linewidth=settings.linewidth_all
legend_location = settings.legend_location
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
dashes=settings.dashes
color_highlight=settings.color_highlight
color_highlight_obs=settings.color_highlight_obs
if IncF.b_TS_hlmodel:
    ls=ls[:len(IncF.c_TS_hlmodel)]+ls
    color2=color_highlight[:len(IncF.c_TS_hlmodel)]+color2
if IncF.b_TS_hlobs:
    ls=ls[:len(IncF.c_TS_hlobs)]+ls
    color1=color_highlight_obs[:len(IncF.c_TS_hlobs)]+color1
    

def make_stats(data_models, data_obs, column_obs, freq_stats):
    ## Substract models with obs
    models_subs=data_models.sub(data_obs[column_obs], axis=0)
    ## Add models with obs
    models_add=data_models.add(data_obs[column_obs], axis=0)
    ## Obs average
    obs_avg=data_obs.copy()
    obs_avg[column_obs] = obs_avg.resample(freq_stats)[column_obs].transform('mean')
    ## Create model Averate
    models_avg=models_add.copy()
    for column in models_add.columns:
        models_avg[column] = models_avg.resample(freq_stats)[column].transform('mean')
    ## Sigamas
    models_sigma=np.sqrt( ( (data_models-models_avg)**2 ).groupby(pd.Grouper(freq=freq_stats)).mean())
    obs_sigma=np.sqrt( ( (data_obs-obs_avg)**2 ).groupby(pd.Grouper(freq=freq_stats)).mean())

    ################
    data_bias=models_subs.groupby(pd.Grouper(freq=freq_stats)).mean()
    data_rmse=(models_subs**2).groupby(pd.Grouper(freq=freq_stats)).mean()
    data_nmbias=(2*models_subs/models_add).groupby(pd.Grouper(freq=freq_stats)).mean()
    data_fge=abs(2*models_subs/models_add).groupby(pd.Grouper(freq=freq_stats)).mean()
    data_r=((data_models-models_avg).mul((data_obs-obs_avg)[column_obs], axis=0)).groupby(pd.Grouper(freq=freq_stats)).mean()/(models_sigma.mul(obs_sigma[column_obs], axis=0))
 
    return data_bias, data_rmse, data_nmbias, data_bias, data_fge, data_r


#####################################
###################################
## Plot Two Sides var
def plot_basic(name, data_or, color_dict, timeframe, var, units, direc="", legend_title="", freq_stats=""):
    ## Make sure data is for timeframe asked
    date1="".join(timeframe[0].split("-")[::-1])
    date2="".join(timeframe[1].split("-")[::-1])
    data=data_or.loc[date1:date2]

    ## Check wether it's the same variable that is being plotted
    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    equiv=[k for  k, v in variables_equivalences.items() if var in v][0]


    
    ## Prepare data
    #breakpoint()
    case, color1_0, color2_0,  data_obs, data_models = prepare_data_one_var(data_or, color_dict)
    if data_obs.dropna().empty:
        return None
    if data_models.dropna(thresh=1).empty:
        return None 
    #data_obs.columns=[" ".join(c.split(" ")[1:]) for c in data_obs.columns]
    #data_models.columns=[" ".join(c.split(" ")[1:]) for c in data_models.columns]
    column_obs=data_obs.columns[0]
    """
    ## Substract models with obs
    models_subs=data_models.sub(data_obs[column_obs], axis=0)
    ## Add models with obs
    models_add=data_models.add(data_obs[column_obs], axis=0)
    ## Obs average
    obs_avg=data_obs.copy()
    obs_avg[column_obs] = obs_avg.resample(freq_stats)[column_obs].transform('mean')
    ## Create model Averate
    models_avg=models_add.copy()
    for column in models_add.columns:
        models_avg[column] = models_avg.resample(freq_stats)[column].transform('mean')
    ## Sigamas
    models_sigma=np.sqrt( ( (data_models-models_avg)**2 ).groupby(pd.Grouper(freq=freq_stats)).mean())
    obs_sigma=np.sqrt( ( (data_obs-obs_avg)**2 ).groupby(pd.Grouper(freq=freq_stats)).mean())

    ################
    data_bias=models_subs.groupby(pd.Grouper(freq=freq_stats)).mean()
    data_rmse=(models_subs**2).groupby(pd.Grouper(freq=freq_stats)).mean()
    data_nmbias=(2*models_subs/models_add).groupby(pd.Grouper(freq=freq_stats)).mean()
    data_fge=abs(2*models_subs/models_add).groupby(pd.Grouper(freq=freq_stats)).mean()
    data_r=((data_models-models_avg).mul((data_obs-obs_avg)[column_obs], axis=0)).groupby(pd.Grouper(freq=freq_stats)).mean()/(models_sigma.mul(obs_sigma[column_obs], axis=0))

    #breakpoint()
    """
    data_bias, data_rmse, data_nmbias, data_bias, data_fge, data_r = make_stats(data_models, data_obs, column_obs, freq_stats)

    ###############
    ### Plot
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(nrows=5, sharex=True,  figsize=(figsize[0], figsize[1]+0.6*figsize[1]*(5-1)))
    fig.subplots_adjust(wspace=0, hspace=0.05)
    #fig.suptitle(f"{var} [{ aux_plots.get_pretty_units(units[var])}]", fontsize=titlesize)
    ax1.set_title(f"{var} ["+aux_plots.get_pretty_units(units[var])+"]\n" , fontsize=titlesize)
    ## Plot
    data_r.plot(color=[color_dict[c] for c in data_r.columns], ax=ax1, style=ls, ms=ms, lw=linewidth)
    data_fge.plot(color=[color_dict[c] for c in data_fge.columns], ax=ax2, style=ls, ms=ms, lw=linewidth)
    data_nmbias.plot(color=[color_dict[c] for c in data_nmbias.columns], ax=ax3, style=ls, ms=ms, lw=linewidth)
    data_bias.plot(color=[color_dict[c] for c in data_bias.columns], ax=ax4, style=ls, ms=ms, lw=linewidth)
    data_rmse.plot(color=[color_dict[c] for c in data_rmse.columns], ax=ax5, style=ls, ms=ms, lw=linewidth)


    labels=["Corr. Coef.", "FGE", "MNBIAS", "BIAS", "RMSE"]
    for ax, label in zip([ax1, ax2, ax3, ax4, ax5], labels):
        handles, labels = ax.get_legend_handles_labels()
        if IncF.b_TS_hl1stmodel:
            highlight_first_line(ax, labels[0])
        #print("labels:", labels)
        labels = [" ".join(l.split(" ")[1:-1]) for l in labels]
        do_legend(fig, ax, handles, labels, "")
        align_yaxis(ax, None, BoolTwoVars=False, var1=var, var2=var, nyticks=3)
        ax.set_ylabel(label)
        #breakpoint()
 
        if IncF.b_TS_hlmodel:
            highlight_first_line(ax, labels, IncF.c_TS_hlmodel, case="Model")

        align_xaxis(ax5)
    fig.align_labels()

    ## Save Figure
    name="%s___%s-Statistics___%s-%s"%(name, var, date1, date2)
    aux_plots.savefigure(fig, direc, name)

    ## Close Figures
    fig.clf()
    plt.close(fig)
    plt.close('all')
    return None



#####################################
## Make Plot. It redirects to the appropriate one depending on what is in vars_to_plot
def make_plots(name, data, vars_to_plot, units, direc="", plot_freq="", case="", legend_title="", variables=None, freq_stats=""):

    makedirs(direc, exist_ok = True) 
    makedirs(join(direc, "metadata"), exist_ok = True) 

    ## Get variable that was asked to plot
    var = sorted(list(set(aux.flat(vars_to_plot, out=[]))))[0]
    ## Colors for each line/column
    color_dict=get_color_columns(data.columns)

    #####################
    ## Plots with the frequency asked
    if plot_freq=="Whole":
        data["groupby"]=1
        grouped_data=data.groupby("groupby")
    else:
        grouped_data=data.groupby(pd.Grouper(freq=plot_freq))

    ## Go through frequencies asked, called group
    for group in list(grouped_data.groups.keys()) :
        ## Retrieve group data
        group_data=grouped_data.get_group(group)
        ## Remove times that are not shared if it's requested to do so.
        if b_showshared:
            group_data=group_data.dropna(thresh=2)
        ## If data is empty then skip plotting things
        if group_data.empty:
            continue

        ## Remove groupby column
        group_data=group_data[[c for c in group_data.columns if c !="groupby"]]

        ## Time frame
        date1=group_data.index[0].strftime("%d-%m-%Y")  
        date2=group_data.index[-1].strftime("%d-%m-%Y")
        timeframe=[date1, date2]

        ## Plot
        plot_basic(name, group_data, color_dict, timeframe, var, units, direc=direc, legend_title=legend_title, freq_stats=freq_stats)



###############################################################
###############################################################
#####################################
## TODO: Fix code to include timeseries of 3d models
def main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData):
    if IncF.i_TS_stats==0:
        print("...Plotting was not requested.")
        return None
    if t_Stations_Info.empty:
        print("...No stations were read. Skipping this step.")
        return None     
    #breakpoint()
  
    ###########
    ## Start going over specifications for plots 
    path_base=join(IncF.c_FigureDir, "TimeSeries-Statistics")
    makedirs(path_base, exist_ok = True) 
    for count, plot_freq in enumerate(IncF.c_TS_stats_freq):
        freq, freq_stats, plot_freq =plot_freq.split(",")
        path_freq=join(path_base, f"{freq},{freq_stats},{plot_freq}")
        makedirs(path_freq, exist_ok = True) 

        ###########
        for vars_plot in IncF.c_TS_stats_Plots:
            #breakpoint()
            vars_to_plot=aux.find_equiv_vars(vars_plot)
            print("...Plotting:  ", vars_plot, "  for freq   ", freq, " and time period:  ", plot_freq)
            #breakpoint()
            if len(vars_to_plot)==0:
                print("   ...Skipping. Network requested did not have data.")
                continue


            #################
            data, data_obs, data_mod = get_new_data(vars_to_plot, t_ObsStationData[freq], t_ModelStationData[freq], plot_freq)
            units=data.meta.units
            
            variables=set(aux.flat(vars_to_plot, out=[]))
            networks=sorted(list(set([c[0].split()[1] for c in data.columns if "Obs" in c[0]])))
            network=list(vars_to_plot["O"].keys())[0]
            models=sorted(list(set([c[0].split()[1] for c in data.columns if "Model" in c[0]])))
            mods="-".join(models)


            if data_obs.empty:
                print("   ...Skipping this step. No observation data was found. ")
                continue
            if data_mod.empty:
                print("   ...Skipping this step. No model data was found. ")
                continue
            if len(variables)>=2:
                print("   ...Skipping this step. Plot of statistics with more than one variables is not implemented.")
                continue
            if len(networks)>=2:
                print("   ...Skipping this step. Plot of statistics with more than one network cant be processed yet.")
                continue




            ####################
            ## Plot
            direc=join(path_freq, "-".join(networks))
            makedirs(direc, exist_ok = True) 
            makedirs(join(direc, "metadata"), exist_ok = True) 


            ##    
            info=aux.get_info_group(t_Stations_Info, network)
            minif=info.copy()
            minif=minif[minif["ID"].str.contains('^f[0-9]-+')]
            mini=info[~info["ID"].isin(minif["ID"])]
            filters=t_Stations_Filters[network]
            mods="-".join(models)


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
                    loc_mod="\n".join(ID_mod.apply(lambda x:  "(%s)  %s    Distance: %.3f km"%(x.location, x.level_0, x.Distance.km), axis=1).values)

                    ############
                    name="%s--%s_%s"%(station, "-".join(info[info["ID"]==station]["Nombre"].values[0].split()).replace(".", "").replace(",", ""), mods)
                    if mods=="":
                        name=name[:-1]
                    legend_title="%s:  %s"%(networks[0], info[info["ID"]==station]["Nombre"].values[0])
                    locations=[legend_title, loc_obs, loc_mod]
                    make_plots(name, data_station, vars_to_plot,  units, direc=direc, plot_freq=plot_freq, legend_title=locations, freq_stats=freq_stats)


            ##############
            ## filter plots
            if bool_plotfilters:
                for station in minif["ID"]:
                    try:
                        data_station=data.xs(station, axis=1, level=1)
                    except:
                        continue
                    if data_station.dropna(thresh=1).empty:
                        continue

                    ############
                    name="%s_%s"%(station, mods)
                    if mods=="":
                        name=name[:-1]
                    legend_title="%s:  %s"%(networks[0], station)
                    locations=[legend_title]
                    make_plots(name, data_station, vars_to_plot,  units, direc=join(direc, station), plot_freq=plot_freq, legend_title=locations, freq_stats=freq_stats)


