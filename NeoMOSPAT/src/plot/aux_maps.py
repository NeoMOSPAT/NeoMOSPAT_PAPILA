# -*- coding: utf-8 -*-
import numpy as np
from numpy import arange, floor, ceil
#
import src.IncludeFile as IncF
import src.common as aux
import src.plot.aux_plots as aux_plots
import src.read.read_shps as shps
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from src.write.statistics import make_stats
#
from src.plot.ticker_formatters import LatDMS, LonDMS
#
import contextily as ctx
#
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, CenteredNorm, Normalize, TwoSlopeNorm, LogNorm

###############
## Add shapefiles to plot and keep limit of interest
settings_shp=aux_plots.ShapefileMapConfig
shapefiles=shps.read_shp_maps(crs_out="epsg:4326")
settings_shp=aux_plots.ShapefileMapConfig
###############
## Map settings
settings = aux_plots.BaseMapConfig
figarea = settings.image_area
symmetric_centering=settings.symmetric_centering
cmap_cyclic = settings.cmap_cyclic
cmap_divergent = settings.cmap_divergent
cmap_seq = settings.cmap_seq
cmap_stats = settings.cmap_stats
ticks_cb = settings.ticks_cb
## Satelite background confg
alpha_min = settings.alpha_min
alpha_max = settings.alpha_max
## Stations plot
ms = settings.ms
markerlist = settings.markerlist


###############
## Creates figure size
def get_figsize(bounds, figarea=figarea):
    ratio=(bounds[1]-bounds[3])/(bounds[0]-bounds[2])
    widen_ratio=1
    widen_constant=1.5
    if ratio>1:
        ratio=1/ratio
        fig_size=(widen_ratio*np.sqrt(figarea)*ratio, np.sqrt(figarea)/ratio/widen_ratio)
    else:
        fig_size=(widen_ratio*np.sqrt(figarea)/ratio, np.sqrt(figarea)*ratio/widen_ratio)

    fig_size_new=(widen_ratio*fig_size[0]+widen_constant, fig_size[1]*fig_size[0]/(fig_size[0]*widen_ratio+widen_constant))

    return fig_size_new


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


##############
def clippedcolorbar(fig, CS, **kwargs):
    vmin = CS.get_clim()[0]
    vmax = CS.get_clim()[1]
    m = ScalarMappable(cmap=CS.get_cmap())
    m.set_array(CS.get_array())
    m.set_clim(CS.get_clim())
    step = CS.levels[1] - CS.levels[0]
    cliplower = CS.zmin<vmin
    clipupper = CS.zmax>vmax
    noextend = 'extend' in kwargs.keys() and kwargs['extend']=='neither'
    # set the colorbar boundaries
    #boundaries = arange((floor(vmin/step)-1+1*(cliplower and noextend))*step, (ceil(vmax/step)+1-1*(clipupper and noextend))*step, step)
    boundaries = arange((floor(vmin/step)+1*(cliplower and noextend))*step, (ceil(vmax/step)-1*(clipupper and noextend))*step, step) 
    kwargs['boundaries'] = boundaries
    # if the z-values are outside the colorbar range, add extend marker(s)
    # This behavior can be disabled by providing extend='neither' to the function call
    if not('extend' in kwargs.keys()) or kwargs['extend'] in ['min','max']:
        extend_min = cliplower or ( 'extend' in kwargs.keys() and kwargs['extend']=='min' )
        extend_max = clipupper or ( 'extend' in kwargs.keys() and kwargs['extend']=='max' )
        if extend_min and extend_max:
            kwargs['extend'] = 'both'
        elif extend_min:
            kwargs['extend'] = 'min'
        elif extend_max:
            kwargs['extend'] = 'max'
    return fig.colorbar(m, **kwargs)


###############
def add_colorbar(fig, cax, vmin, vmax, norm, extend):
    bool_cero = vmin<0 and vmax>0
    norm.set_clim(vmin, vmax) 
    if bool_cero:
        try:
            if symmetric_centering:
                cb=clippedcolorbar(fig, norm, cax=cax, norm=CenteredNorm(vcenter=0, halfrange=max(abs(vmin), vmax)))
            else:
                cb=clippedcolorbar(fig, norm, cax=cax, norm=TwoSlopeNorm(0, vmin=vmin, vmax=vmax))
        except:
            if symmetric_centering:
                cb=fig.colorbar(norm, cax=cax, extend=extend, norm=CenteredNorm(vcenter=0, halfrange=max(abs(vmin), vmax)))
            else:
                cb=fig.colorbar(norm, cax=cax, extend=extend, norm=TwoSlopeNorm(0, vmin=vmin, vmax=vmax))
    else:
        try:
            cb=clippedcolorbar(fig, norm, cax=cax)
        except:
            cb=fig.colorbar(norm, cax=cax, extend=extend)


###############
## Format maps (title, colorbar, etc)
def format_map(fig, ax, model, var, unit, time, extend, norm, vmin, vmax, type_plot=""):
    fig.patch.set_alpha(0)
    ax.set_facecolor('white')
    ax.axes.set_aspect('equal')
    norm.get_cmap()._rgba_bad = (1.0, 1.0, 1.0, 0)

    ## Title
    ax.set_title("%s\n%s [%s]    %s"%(model, var, aux_plots.get_pretty_units(unit),  time), pad=10)

    ## Format ticks
    format_ticks(ax=ax)

    ## Color Bar Axis
    cax=fig.add_axes([ax.get_position().x1+0.02,ax.get_position().y0,0.035,ax.get_position().height])  
    #divider = make_axes_locatable(ax)
    #cax = divider.new_horizontal(size="5%", pad=0.2)
    ## Add Color bar
    add_colorbar(fig, cax, vmin, vmax, norm, extend)
    return None

###############
def prepare_colormaps(alpha, vmin, vmax, var, map_type, boolLog=False):
    if not boolLog:
        ## Get type of cmap 
        if var in ["WDIR", "WDIR20m", "V0", "VD", "time-max"]: 
            cmap = cmap_cyclic
            norm = Normalize(vmin=vmin, vmax=vmax)
        ## Divergent
        if vmin<0 and vmax>0:
            cmap = cmap_divergent
            if symmetric_centering:
                norm = CenteredNorm(vcenter=0, halfrange=max(abs(vmin), vmax))
            else:
                norm = TwoSlopeNorm(0, vmin=vmin, vmax=vmax)
        ## Sequential
        else:
            cmap = cmap_seq
            norm = Normalize(vmin=vmin, vmax=vmax)

        ##########
        ## If background was added then add alpha channel
        if alpha:
            cmap = custom_cmap(cmap)
    else:
        ## Get type of cmap 
        if var in ["WDIR", "WDIR20m", "V0", "VD"]: 
            cmap = cmap_cyclic
            norm = Normalize(vmin=vmin, vmax=vmax)
        ## Divergent
        else:
            cmap = cmap_seq
            norm = LogNorm(vmin=vmin, vmax=vmax)
    return cmap, norm


###############
## Prepare Maps (Adds background and shapefiles)
def prepare_maps(bounds, ax=None, background=IncF.i_Map_Background):
    ax=ax or plt.gca()

    ## Change this eventually for it to have more specifications regarding maps
    for shp in shapefiles:
        shapefiles[shp].plot(ax=ax, **getattr(settings_shp, shp), zorder=200)

    ## Add Shapefiles that will be added
    ax.set_ylim(bounds[1], bounds[3])
    ax.set_xlim(bounds[0], bounds[2])

    ## Add background
    alpha=True
    if background==1:
        ctx.add_basemap(ax=ax, source=ctx.providers.Esri.WorldImagery, crs='epsg:4326')
    elif background==2:
        ctx.add_basemap(ax=ax, source=ctx.providers.OpenStreetMap.Mapnik, crs='epsg:4326')
    elif background==3:
        ctx.add_basemap(ax=ax, source=ctx.providers.Stamen.TerrainBackground, crs='epsg:4326')
    else:
        alpha=False
        
    return ax, alpha


###############
## Create map with alpha  
def custom_cmap(cmap, ncolors=256, alpha_min=alpha_min, alpha_max=alpha_max):
    ## Get colormap 
    color_array = plt.get_cmap(cmap)(range(ncolors))
    ## Change alpha values 
    color_array[:,-1] = np.linspace(alpha_min, alpha_max, ncolors)**2

    ## Create colormap object
    map_object = LinearSegmentedColormap.from_list(name=f'{cmap}_alpha', colors=color_array)

    return map_object


###############
## Finds vmin and vmax values depending on what was asked on IncludeFile 
def find_range_cb(z, z_global, cb):
    #breakpoint()
    ## Standard of maximum per day
    vmin, vmax = np.nanmin(z), np.nanmax(z)
    if cb is not None:
        ## Range given
        if len(cb)==2:
            vmin, vmax = cb
            boolLog=False
        elif len(cb)==3:
            vmin, vmax, boolLog = cb
        else:
            print("      ...Invalid defintion of cb: ", cb)
            vmin = vmax = 0
            BoolLog = False
        ## Percentile option for minimum value
        if isinstance(vmin, (str, bytes)):
            minper=float(vmin.split(":")[1])
            type_min=vmin.split(":")[0].split("_")[1]
            if type_min=="all":
                vmin = np.nanpercentile(z_global, minper)
            elif type_min=="freq":
                vmin = np.nanpercentile(z, minper)
        ## Percentile option for maximum value
        if isinstance(vmax, (str, bytes)):
            maxper=float(vmax.split(":")[1])
            type_max=vmax.split(":")[0].split("_")[1]
            if type_max=="all":
                vmax = np.nanpercentile(z_global, maxper)
            elif type_max=="freq":
                vmax = np.nanpercentile(z, maxper)
    return vmin, vmax, boolLog
 

###############
## Finds the appropriate option for variable "extend" in fig colorbar
def find_extend(vmin, vmax, datamin, datamax):
    if datamin >= vmin:
        if datamax <= vmax:
            extend="neither"
        else:
            extend="max"
    else:
        if datamax <= vmax:
            extend="min"
        else:
            extend="both"
    return extend


#############
def prepare_statistics(obs_data, obs_info, mod_data, model_info, vars_model):
    data={}
    for var in vars_model:
        data[var]={}
        data_network={}
        for network in obs_data:
            data_all=[]
            for ID in obs_info.loc[network]["ID"]:
                if ID in mod_data[network][var].columns and ID in obs_data[network][var].columns:
                    df=make_stats(obs_data[network][var][[ID]], mod_data[network][var][[ID]], ID).drop(columns=["Obs"]).T
                    data_all.append(df)
            stats_network=pd.concat(data_all, keys=list(obs_info.loc[network]["ID"].values))
            data[var][network]=stats_network
        
    return data

###############
## Add statistics to maps
def add_statistics(t_ObsStationData, t_Stations_Info, var, day, ax=None, cax=None, case=0, levels=8, extend="", fig=None):

    ax=ax or plt.gca()
    #cax=cax or plt.gca()
    #cax = fig.add_axes([ax.get_position().x0, ax.get_position().y1, ax.get_position().width, 0.035])
    
    if t_ObsStationData=={}:
        print("      ...Requested to plot stations but no stations were read")
        return None

    var_name=var
    sinfo=t_Stations_Info.reset_index(level=1, drop=True)   
    summary = pd.concat([t_ObsStationData[g][var] for g in t_ObsStationData]).describe()
    vmin = summary["min"]
    vmax = summary["max"]
    for index, group in enumerate(sinfo.index.drop_duplicates()):
        data=t_ObsStationData[group]
        variables_model=aux.equiv_vars(data.columns.get_level_values(0))
        data=data[[var_name]].reset_index()
        data=data[["level_0", var_name]]
        data.columns=["ID", var_name]
        data=sinfo.loc[[group]][["ID", "Latitud", "Longitud", "geometry"]].merge(data, on="ID").dropna()

        if case==2:
            #psta=data.plot(ax=ax, cax=cax, column=var_name, marker="o", markersize=int(1.4*ms), cmap=cmap, edgecolor="black", zorder=200)
            ps=ax.scatter(data["Longitud"].values, data["Latitud"].values, c=data[var_name].values, s=int(1.4*ms), cmap=cmap_stats, vmin=vmin, vmax=vmax, edgecolor="black", zorder=200)
            fig.colorbar(ps, cax=cax)
        elif case==3:   
            ps=ax.scatter(data["Longitud"].values, data["Latitud"].values, c=data[var_name], s=int(1.4*ms), cmap=cmap_stats, vmin=vmin, vmax=vmax, edgecolor="black", norm=BoundaryNorm(np.linspace(vmin, vmax, levels+1), ncolors=colormap.N), zorder=200)
            fig.colorbar(ps, cax=cax)

    ax.axes.set_aspect('equal')
    return vmin, vmax



###############
## Add stations to maps
def add_stations(t_ObsStationData, t_Stations_Info, var, day, vmin, vmax, cmap, ax=None, case=0, levels=8):

    ax=ax or plt.gca()
    if t_ObsStationData=={}:
        print("      ...Requested to plot stations but no stations were read")
        return None
    ## Do stuff   
    try:
        var_name, var_op = var.split("-")
        time=False
    except:
        var_name, time, var_op = var.split("-")
        time=True

    sinfo=t_Stations_Info.reset_index(level=1, drop=True)    
    for index, group in enumerate(sinfo.index.drop_duplicates()):
        data=t_ObsStationData[group]
        variables_model=aux.equiv_vars(data.columns.get_level_values(0))
        if var_name in variables_model:
            var_obs=[v for v in set(data.columns.get_level_values(0)) if var_name in aux.equiv_vars([v])][0]

        
            data=data[var_name].loc[var_op].reset_index()
            data.columns=["ID", var_name]
            data=sinfo.loc[[group]][["ID", "Latitud", "Longitud", "geometry"]].merge(data, on="ID").dropna()
    
            if case==1 :
                data.plot(ax=ax, color="black", marker=markerlist[index], markersize=ms, zorder=200)
            elif case==2 and time is False:
                data.plot(ax=ax, column=var_name, marker="o", markersize=int(1.4*ms), cmap=cmap, vmin=vmin, vmax=vmax, edgecolor="black", zorder=200)
            elif case==3 and time is False:   
                data.plot(ax=ax, column=var_name, marker="o", markersize=int(1.4*ms), cmap=cmap, edgecolor="black", norm=BoundaryNorm(np.linspace(vmin, vmax, levels+1), ncolors=colormap.N), vmin=vmin, vmax=vmax, zorder=200)
                
    ax.axes.set_aspect('equal')
    return ax
