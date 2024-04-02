# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import numpy as np
import pandas as pd
import geopandas as gpd
#
import matplotlib.pyplot as plt
import matplotlib as mpl
#
import src.common as aux
import src.plot.aux_plots as aux_plots
import src.plot.aux_maps as aux_maps
from datetime import datetime
from src.plot.ticker_formatters import LatDMS, LonDMS
from src.plot.Surface_Maps import *
#
from os.path import join
from os import makedirs
from  matplotlib.colors import LogNorm, Normalize, BoundaryNorm
from matplotlib.ticker import LogLocator, LogFormatterSciNotation
#


# "ENS-MedianPython"
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

world=gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
minx=-110.0
miny=-58.0
maxx=-33.5
maxy=25.0
bounds = minx, miny, maxx, maxy

###############
## Change ticks():
def format_ticks(ax=None):
    ax=ax or plt.gca()

    nticksx = len(ax.get_xticks())
    nticksy = len(ax.get_yticks())
    nticks=min(nticksx, nticksy)

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    range_y=ymax-ymin
    range_x=xmax-xmin
    prop=range_y/range_x

    
    ## Verticla rectangle
    if prop>=5/3:
        ax.xaxis.set_major_locator(plt.MaxNLocator(3)) 
        ax.yaxis.set_major_locator(plt.MaxNLocator(5)) 
    ## Horizontal rectangle
    elif prop<=3/5:
        ax.yaxis.set_major_locator(plt.MaxNLocator(3)) 
        ax.xaxis.set_major_locator(plt.MaxNLocator(5)) 
    ## Squareish
    else:
        if nticks >=5:
            nbins=5
        else:
            nbins=3
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins)) 
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins)) 
    

    ## Set ticks
    xticks = ax.get_xticks().tolist()[1:-1]
    yticks = ax.get_yticks().tolist()[1:-1]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)

    ## Modify labels
    xlabels = [LonDMS(t).to_string(resolution='minutes').replace("00'", "")  for t in xticks]
    ylabels = [LatDMS(t).to_string(resolution='minutes').replace("00'", "") for t in yticks]
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)

    return ax


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
def main(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info):
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
  
    ###########
    ## Start going over specifications for plots 
    path_base=join(IncF.c_FigureDir, "Maps_Mosaic")
    makedirs(path_base, exist_ok = True) 
    for count, plot_freq in enumerate(IncF.c_Map_freq):
        freq, plot_freq =plot_freq.split(",")

        path_freq=join(path_base, f"{freq},{plot_freq}")
        makedirs(path_freq, exist_ok = True) 

        ###########
        for vars_plot in IncF.c_Map_mosaic:
            print("...Plotting:  ", vars_plot, "  for freq   ", freq, " and time period:  ", plot_freq)
            vars_to_plot=aux.find_equiv_vars(vars_plot[:-1])
            variables=set(aux.flat(vars_to_plot, out=[]))


            #########################################
            mplot, mcol, mrow = vars_plot[-1].split(";")
            type_plot=int(mplot.split(":")[-1])
            type_col=mcol.split(":")[-1]
            type_row=mrow.split(":")[-1]
            c=0
            for start, end in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                c=c+1
                d_start = datetime(int(start.split("-")[-1]), int(start.split("-")[1]), int(start.split("-")[0]))
                d_end = datetime(int(end.split("-")[-1]), int(end.split("-")[1]), int(end.split("-")[0]))
                timeframe = f"{d_start.strftime('%Y-%m-%d')}-{d_end.strftime('%Y-%m-%d')}"
                nrow=2
                ncol=3
                for type_plot in ["Square", "CTRF"]:
                    for scale in ["Normal", "Log"]:
                        row = list(variables)[0]
                        vals_cont=pd.concat([t_ModelData[freq][model].reset_index().set_index("time").loc[d_start.strftime("%Y-%m-%d"):d_end.strftime("%Y-%m-%d")][row] for model in t_ModelData[freq].keys()]).values
                        if type_plot=="Square":
                            if scale=="Log":
                                vmax=np.nanpercentile(vals_cont, 99)
                                vmin = np.nanmin(vals_cont[np.nonzero(vals_cont)])
                                norm=LogNorm(vmin=vmin, vmax=vmax)
                            else:
                                vmax=np.nanpercentile(vals_cont, 97)
                                vmin=0
                                norm=Normalize(vmin=vmin, vmax=vmax)
                            
                            #norm = Normalize(vmin=vmin, vmax=vmax)
                        if type_plot=="CTRF":
                            nlevels=13
                            cmap_ctrf = "plasma"
                            if scale=="Log":
                                vmax=np.nanpercentile(vals_cont, 99)
                                
                                levels = np.logspace(np.log2(vmax)-8, np.log2(vmax), num=nlevels, base=2.0)
                                #levels[0] = 0
                                levels[~np.isfinite(levels)] = 0
                                if 0 not in list(levels):
                                    levels = [0]+list(levels)
                                levels = sorted(list(set(levels)))
                                nlevels=len(levels)
                                vmin = 0
                                cmap_ctrf = cm.get_cmap('plasma', nlevels)
                                norm = BoundaryNorm(levels, nlevels-1)
                            else:
                                vmax = np.nanpercentile(vals_cont, 97)
                                vmin=0
                                levels = np.linspace(vmin, vmax, nlevels)
                                norm = BoundaryNorm(levels, 256, extend="max")
                            

                           
                        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, sharey='row', figsize=(14/2*nrow+1,3*ncol))
                        i=0
                        j=0
                        for model in vars_to_plot["M"]:
                            col = model
                            ## Vars that need to be done 
                            vars_model=vars_to_plot["M"][model]
                            ## Add lat and lon data to data set
                            model_info=t_Models_Info[model].sort_values(by=["i_lon", "i_lat"])
                            units=t_ModelData[freq][model].meta.units
                            data_resampled=t_ModelData[freq][model].reset_index(level=[0,1]).loc[d_start.strftime("%Y-%m-%d"):d_end.strftime("%Y-%m-%d")]
                            data_resampled=model_info.merge(data_resampled, on=["i_lat", "i_lon"], how="outer")#.set_index(data_resampled.index)
                            filters=t_Model_Filters[model]
                            ## Grid
                            ilats=model_info["i_lat"].drop_duplicates(keep="first")
                            ilons=model_info["i_lon"].drop_duplicates(keep="first")
                            lat_min=ilats.min()
                            lon_min=ilons.min()
                            lat_max=ilats.max()
                            lon_max=ilons.max()
                            len_x = int(lon_max+1 - lon_min)
                            len_y = int(lat_max+1 - lat_min)


                            data_resampled.sort_values(by=["i_lon", "i_lat"], inplace=True)
                            x=model_info["lon"].values.reshape((len_x, len_y)) 
                            y=model_info["lat"].values.reshape((len_x, len_y)) 

                            data = data_resampled

                            data = data.sort_values(by=["i_lon", "i_lat"])
                            data = data.groupby(["i_lon", "i_lat"]).mean()
                            z = data[row].values.reshape((len_x, len_y)) 
                            unit=units[row]
                            if type_plot=="Square":
                                cs=axes[i,j].pcolormesh(x, y, z, shading="nearest", cmap="plasma", zorder=1, snap=True,  linewidth=0, rasterized=True, antialiased=True, norm=norm)
                            if type_plot=="CTRF":
                                cs = axes[i,j].contourf(x, y, z, cmap=cmap_ctrf, zorder=1, antialiased=True, levels=levels, extend="max", norm=norm)
                            world.boundary.plot(ax=axes[i,j], color="black")
                            axes[i,j].set_ylim(bounds[1], bounds[3])
                            axes[i,j].set_xlim(bounds[0], bounds[2])
                            axes[i,j].axes.set_aspect('equal')
                            axes[i,j].set_title("%s  %s [%s]"%(model, row, aux_plots.get_pretty_units(unit)), pad=10)
                            format_ticks(ax=axes[i,j])
                            j = j+1
                            if j>=3:
                                i=i+1
                                j=0
                                if i>=2:
                                    i=0
                                    j=0
                        if timeframe=="2015-01-01-2015-01-31":
                            timeframe="Jan"
                        if timeframe=="2015-07-01-2015-07-31":
                            timeframe="Jul"                            
                        fig.subplots_adjust(wspace=0, hspace=0.2)
                        #ax, kw = mpl.colorbar.make_axes(axes[i,:])
                        #cax = fig.add_axes([axes[i,ncol-1].get_position().x1+0.02, axes[i,ncol-1].get_position().y0, 0.035, axes[i,ncol-1].get_position().height]) 
                        
                        if scale=="Log" and type_plot=="Square":
                            locator = LogLocator(base=2)
                            formatter=LogFormatterSciNotation(base=2)
                            fig.colorbar(cs, ax =axes.ravel().tolist(), extend="max", spacing='uniform', ticks=locator, format=formatter)
                        else:
                            fig.colorbar(cs, ax=axes.ravel().tolist(), extend="max", spacing='uniform')
                        #breakpoint()
                        aux_plots.savefigure(fig, path_freq, f"{row}_{type_plot}_{scale}_grid_models_{timeframe}.png")

                        #breakpoint()

     