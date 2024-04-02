# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import numpy as np
import pandas as pd
#
import matplotlib.pyplot as plt
#
import src.common as aux
from src.plot.aux_timeseries import *
from src.plot.Time_Series import *
import src.plot.aux_plots as aux_plots
#
from os.path import join
from os import makedirs
#

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



###############################################################
def apply_filter_obs(data_obs, t_Stations_Filters, filt):
    columns = data_obs.columns
    cols_network = {net:[ c for c in columns  if net in c[0]] for net in t_Stations_Filters}
    
    dfs=[]
    nfilts=[]
    cols=[]
    for net, cols_net in cols_network.items():
        filters=t_Stations_Filters[net][filt]
        columns_names = [filt+"-"+c for c in filters]
        nfilts.append(len(filters))
        df=data_obs[cols_net]
        select = df.columns.get_level_values(1).isin(columns_names)
        cols=cols+columns_names
        df = df.loc[:, select]#.xs(columns_names, level=1, axis=1) 
        dfs.append(df)

    return pd.concat(dfs), max(nfilts), list(set(cols))


###############################################################
def get_vars_obs(data_obs, t_Stations_Filters, variables):
    return data_obs, len(variables)

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
    path_base=join(IncF.c_FigureDir, "TimeSeries_Mosaic")
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
        for vars_plot in IncF.c_TS_mosaic:

            vars_to_plot=aux.find_equiv_vars(vars_plot[:-1])

            
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

            #########################################
            mplot, mcol, mrow = vars_plot[-1].split(";")
            type_plot=int(mplot.split(":")[-1])
            type_col=mcol.split(":")[-1]
            type_row=mrow.split(":")[-1]

            
            if "f" in type_col:
                data_obs, ncol1, cols = apply_filter_obs(data_obs, t_Stations_Filters, type_col)
                ncol2=ncol1
                #data_mod, ncol2 = apply_filter_mod(data_mod, t_Model_Filters, type_col)
                ncol = max(ncol1, ncol2)

            if "f" in type_row:
                data_obs, nrow1, cols = apply_filter_obs(data_obs, t_Stations_Filters, type_row)
                nrow2=nrow1
                #data_mod, ncol2 = apply_filter_mod(data_mod, t_Model_Filters, type_col)
                nrow = max(nrow1, nrow2)

            if type_row=="Vars":
                nrow=len(variables)
                rows=list(variables)
                #data_obs, nrow1 = get_vars_obs(data_obs, t_Stations_Filters, variables)
                #nrow2=nrow1
                #data_mod, ncol2 = get_vars_mod(data_mod, t_Model_Filters)
                #nrow = max(nrow1, nrow2)
            ## Colors for each line/column
            color_dict=get_color_columns(data_obs.columns)

            #####################
            ## Plots with the frequency asked
            if plot_freq=="Whole":
                data_obs["groupby"]=1
                grouped_data=data_obs.groupby("groupby")
            else:
                grouped_data=data_obs.groupby(pd.Grouper(freq=plot_freq))
                      

            #########
            list_dates=list(grouped_data.groups.keys()) 
            for group in list_dates:
                try:
                    group_data=grouped_data.get_group(group)
                except:
                    group_data=pd.DataFrame()

                if group_data.empty:
                    continue

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, sharey='row', figsize=(19,12))
                for i in range(nrow):
                    row = rows[i]
                    for j in range(ncol):
                        col = cols[j]
                        if type_plot==1:
                            data = data_obs[[c for c in data_obs.columns if row in c[0]]].xs(col, level=1, axis=1)
                           

                            dataRM=data.rolling(window=30,  center=True, min_periods=int(30*0.7)).mean()
                            dataRM.plot(color=aux_plots.lighten_color(color1, amount=1.2), ax=axes[i,j], style="-", ms=ms, zorder=10, xlabel='', lw=3)
                            data.plot(color=color1, ax=axes[i,j], ms=1.5*ms, alpha=0.3, style=".", zorder=1, markeredgecolor="None", markeredgewidth=0, xlabel='')
                            #brealpoint()
                            #data.plot(ax=axes[i,j], style=ls, ms=ms, color=color1, xlabel='')
                            axes[i,j].legend([f"{col.replace('f1-', '')} {row}"], frameon=False)
                        if j==0:
                            axes[i,j].set_ylabel(f"{row} [{aux_plots.get_pretty_units(units[row])}]" )
                            #continue
                        #elif i==nrow-1:
                        #    continue
                        #elif j==0:
                        #    continue                   
                        #else:
                        #    axes[i,j].set_xticklabels([])
                        #    axes[i,j].set_yticklabels([])   

                fig.subplots_adjust(wspace=0, hspace=0)
                aux_plots.savefigure(fig, path_freq, "trial")
                breakpoint()

             