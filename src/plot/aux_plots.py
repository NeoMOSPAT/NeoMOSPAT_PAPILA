# -*- coding: utf-8 -*-
from typing import List

import matplotlib.pyplot as plt
#
from matplotlib.colors import cnames, to_rgb, LinearSegmentedColormap
import colorsys
from os.path import join
import yaml 
###########################################################
###########################################################
##################      SETTINGS      #####################
list_colors1=[
    "black",
    "tab:orange", 
    "goldenrod", 
    "indianred", 
    "darkorchid", 
    "crimson", 
    "orangered", 
    "y", 
    "xkcd:light red",
    "xkcd:dull red", 
    "xkcd:rusty red", 
    "xkcd:rusty orange", 
    "xkcd:tangerine", 
    "xkcd:orange yellow", 
    "xkcd:goldenrod", 
    "xkcd:sun yellow", 
    "xkcd:off yellow", 
    "xkcd:chartreuse"]*5


list_colors2=[
    "tab:red",
    "tab:orange",
    "tab:green",
    "tab:cyan",
    "tab:blue", 
    "tab:purple",
    "indianred", 
    "xkcd:tea",
    "xkcd:grey pink",
    "xkcd:desert",
    "xkcd:velvet",
    "tab:olive",
    "xkcd:dark blue grey",
    "xkcd:muted purple",
    "xkcd:slate green",
    "xkcd:orchid",
    "xkcd:perwinkle blue",
    "xkcd:light sage",
    "tab:cyan",
    "xkcd:ocher", 
    "greenish beige",
    "ocean blue",
    "light grey blue",
    "dark mauve",
    "purple grey",
    "lavender blue",
    "muted pink",
    "coffee",
    "lightseagreen",
    "mediumslateblue",
    "blueviolet",
    "cornflowerblue",
    "mediumorchid",
    "dodgerblue",
    "xkcd:light indigo", 
    "xkcd:dusky blue", 
    "xkcd:light navy", 
    "xkcd:lavender blue", 
    "cadetblue", 
    "xkcd:tealish", 
    "xkcd:dark cyan", 
    "xkcd:twilight", 
    "xkcd:cornflower", 
    "xkcd:deep sky blue", 
    "xkcd:dusty blue", 
    "xkcd:blue grey", 
    "xkcd:light blue", 
    "xkcd:bright sky blue", 
    "xkcd:bright aqua"]*5




########################################
### Definitions for all plots. This for homogeneity
class BaseFigureConfig:
    # TODO: add xticksize and yticksize
    labelsize = 19
    legendsize = 17
    ticksize = 17
    titlesize = 20
    subtitlesize = 16

    # TODO: agregar configuración de escala de colormap que el usuario pueda elegir

########################################
##  Save format kwargs for fig.savefig
class SaveFigureConfig:
    ## Format to save figure
    format="png" 
    ## Removes extra white padding
    bbox_inches="tight"
    ## Transparent background with white background for plot
    facecolor="white"
    edgecolor="none"
    #transparent=True 
    ## If png is used it sets the resolution. 
    ## If vectorized image is used (pdf, svg), it sets the resolution of the imported images to the figure only
    dpi=300


#######################################
###########   BaseMaps   ##############
class BaseMapConfig(BaseFigureConfig):
    format_fig = "png"
    #########
    ## Map Settings
    ## Area that image needs to have. It reshapes its size depending on the ratio of what is being plotted
    image_area=60
    ## Color scheme to use
    cmap_divergent = "coolwarm"
    cmap_seq = "plasma"
    cmap_cyclic = "twilight_shifted"
    cmap_stats = "winter"

    symmetric_centering = True
    ticks_cb = 16
    ## For contourf plot, defines how many intervals to have. 
    nlevels=7
    ## Crop Square of map for only content of map
    BoolCropContent=False
    ## Padding for edge of figure. Adds very small border
    pad_inches=0.01

    ########
    ## Network information
    ## Defines the size of the marker
    ms=55
    ## Defines types of marker for the different networks
    markerlist=["o", "s", "D", "<", ">", "d", "*", "p"]

    ########
    ## Background Contextily
    alpha_min=0.45
    alpha_max=0.65


#######################################
######   How to plot Shapefile   ###### 
class ShapefileMapConfig(BaseFigureConfig):
    ## Defines all the kwargs to give to shapefile.plot() to aesthetics of maps. The variable needs to have the exact name as the shapefile  key value in shp_directories.yaml
    ## "facecolor":"none" doesnt fill the polygon 

    ## Ciudades
    Ciudades = {"facecolor":"none", "edgecolor":"black", "linewidth":1}

    ## Comunas
    Comunas = {"facecolor":"none", "edgecolor":"black", "linewidth":1.1}

    ## Provincias
    Provincias = {"facecolor":"none", "edgecolor":"black", "linewidth":1.2}

    ## Regiones
    Regiones = {"facecolor":"none", "edgecolor":"black", "linewidth":1.2}

    ## Chile
    Chile = {"facecolor":"none", "edgecolor":"black", "linewidth":1.4}

    ## Topo lines every 250m
    Topo_250m = {"facecolor":"none", "edgecolor":"black", "linewidth":0.5, "alpha":0.6}

    ## Lakes
    Lakes = {"facecolor":"dodgerblue", "edgecolor":"dodgerblue", "linewidth":1}

    ## Squares
    Squares = {"facecolor":"none", "edgecolor":"black", "linewidth":1.2}
    
    ## Points
    Points = {"color":"black", "markersize":55}

    ## World
    World = {"facecolor":"none", "edgecolor":"black", "linewidth":1.1}

    

#######################################
###########   Scatter   ###############
class SCTConfig(BaseFigureConfig):
    ## Shared settings
    imagesize=(7, 6)
    markersize=60
    labelsize=17
    legend_location='upper right'
    ## Standard Scatter
    color='tab:orange'
    edgecolor='black'
    ## Density colormap used
    cmap="plasma"
    ## Line Fit color
    linecolor="black"
    ##### Obs vs Obs 
    ## Choose symbol to differentiate network
    smb_x=u'\u2666' # Diamond
    smb_y=u'\u25cf' # Circle
    ## If to pick closest station to station or do all available for the other network
    closest=True



#######################################
########### Time Series ###############
class TSConfig(BaseFigureConfig):
    #########
    ## General Settings
    colors = [list_colors1, list_colors2]
    imagesize = (16, 5.5)
    legend_location = "upper center"
    nele_legend = 8
    linestyle = ["-o"]*80#["--o", ":o"]*80
    linewidth_all=3
    linewidth_highlight=3.5

    ## -, -- , -. , -..
    dashes=[(0, (7, 3)), (0, (7,3,1,3)), (0, (7,3,1,3,1,3)), (0, (7,3,2,5)), (0, (3,1,3,2)), (0, (3,6,3,6,3,18))]*10
    color_highlight=["black", "darkblue"]
    color_highlight_obs=["black", "red"]
    ## Percentaeg of plotting area for which the maximum should be plotted
    perc_top=80
    markersize = 5
    ## Only show shared points 
    b_showshared = False

    ########  One var Obs-Model or var1-var2
    ## Percentage of overlap range so that if plotting same var between to use same axis
    overlap_percentage = 0.4

    ######## Obs-Obs
    ## When asking TS with two networks, specifies whether to plot only closest station to each other or not
    closest = True

    ######
    ## Rolling plots
    ## Minimum percentage of data to have in window to consider the mean of it
    percwindow = 75
    ## Center rolling mean
    center = True



#######################################
########### cross sections ############
class XSSConfig(BaseFigureConfig):
    cmap = 'inferno'
    figheight = 7
    norm = 'LogNorm' # Exact name of the method for the argument norm
    xlabel = 'Coordinates'
    ylabel = 'Height [m]'
    representation = ['pcolormesh']
    # TODO: not implemented yet
    # 'regular_pressure': raw data
    # 'regular_height': data at given heights, maybe from h_start to h_end using a step
    # choose one or both
    data_distribution = ['regular_pressure', 'regular_height']





###########################################################
###########################################################
#################      HELPER FUNC     ####################

########################################
## Lighten or Darken Color
def lighten_color(color, amount=0.5):
    def lighten_single_color(color, amount=0.5):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.
        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        try:
            c = cnames[color]
        except KeyError:
            c = color
        c = colorsys.rgb_to_hls(*to_rgb(c))
        return [colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])]

    if isinstance(color, list):
        new_colors = [lighten_single_color(c, amount=amount)[0] for c in color]
        return new_colors
    else:
        return lighten_single_color(color, amount=amount)



######################################
## Get attrs based on class 
attrs = {k:v for k,v in vars(SaveFigureConfig).items() if "__" not in k}
def savefigure(fig, direc, name, BoolMessage=True, **kwargs):
    ## Append format if its not defined
    try:
        figformat="."+kwargs["format"]
    except:
        figformat="."+attrs["format"]
    name = name if figformat in name else name+figformat

    fig.patch.set_alpha(0)
    ## If kwargs overwrites an option it ke eps kwargs value
    options={**attrs, **kwargs} 


    name=name.replace(" ", "-")
    if BoolMessage:
        print("      ...Saved to: ", name)

    try:
        ## Save figure
        fig.savefig(join(direc, name.replace(" ", "-")), **options)
    except OSError as exc:
        ## Filename too long
        if (exc.errno == 63) or (exc.errno == 36):
            ## Save figure
            name_or=name.split("__")[0].split("/")[-1]
            new_name="LongName_S%s"%len(name_or)
            fig.savefig(join(direc, name.replace(name_or, new_name)), **options)
            print(f"         ...Renamed from  {name}  to  {name.replace(name_or, new_name)}")
        else:
            raise 


######################################
## Get pretty units
## Add units that have mathematical symbols in latex format 
pretty_units={  "ug/m3" :r"$\mu$g/m${}^3$", 
                "W/m2"  :r"W/m${}^2$", 
                "deg.C" :r"$\degree$ C", 
                "m s-1" :"m/s",
                "Deg.M" :r"$\degree$"}
pretty_units = yaml.safe_load(open("config/pretty_units.yaml"))
def get_pretty_units(unit):
    try:
        return pretty_units[unit]
    except:
        return unit


