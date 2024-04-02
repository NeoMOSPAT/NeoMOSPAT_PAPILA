# -*- coding: utf-8 -*-

################################################################
## All parameters should be defined here with default values.
## Their usage should be here and in the README files.
## This files need to be continuously updated to contain all flags needed to run NeoMOSPAT
################################################################
## TODO Iniital Map
i_InitialMap = 0

################################################################
################################################################
###############         General   SET-UP        ################
## Freq for models and obs
## Takes as input all recognized frequencies of pandas (google as Offset aliases)
## https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
## Most commonly used:
##        T     :  Minutes
##        H     :  Hourly
##        D     :  Daily 
##        M, MS :  Monthly
##        W     :  Weekly
##        SM    :  15th and end of month
##        Q, QS :  Quarterly
##        Y, YS :  Yearly
##   Num combin :  3H (every 3 hours), 12H (every 12 hours) or 2D (every 2 days)
c_TimeMeanRes = []

## Time period to read. Format: %d-%m-%Y
c_Start_Date = []
c_Last_Date = []

## Directory used to save figures (doesn't have to exist)
c_FigureDir  = './Figures'
## Dont change this. In case c_RunName is specified it will create folder ./Figures/c_RunName and save everything there. 
## This is recommended to keep figures tidy for different researches or analysis
b_CreateNameSubFolder = True
c_RunName = ""


################################################################
################################################################
###############        Models     SET-UP        ################
## Model Directory folder 
c_ModelDir = "./examples/data/Mod_DB"

## Model name. It has to be a folder within c_ModelDir
c_ModelNames = []
## Alisases. Need to be kept less than 15 characters
c_ModelAlias = []
## Variables to read for each model. Models share variables read.
c_ModelVars = []


## Time Zone (UTC zone)
i_TimeZone = -4  # NOTE: Not used yet
## How often is the model data stored in 
c_ModelTimeFreq = 'DD'
## How to fill empty values 
f_FillValue = "nan" # NOTE: Not used yet


## Whether to only keep model data that is "in land" and discard what is above the ocean
b_OnlyLand = False      #NOTE: Works but it is very slow. Needs to be revised

################################################################
################################################################
###############       Observation    SET-UP     ################
## Observations networks to read. Observations available can be found in config/data_directories.yaml
c_ObsNetwork = []
## c_ObsVars  =  list of lists, one for each network. 
##    Each network list contains the variables to read for each network.
c_ObsVars = []

## Height wanted for MET data, in meters
i_height = 20

## Margin of error allowed for height
i_height_me = 20

## Unkown heights in SINCA are marked as zero. Do you want to consider them if no closer height was found.
b_include_zero = True

## datetime to consider
## 1 : "date_utc"
## 2 : "date_local"
## Default : 2
i_date_column = 1 




################################################################
################################################################
########################      Filters        ###################
### All Stations and area of models defined within the tags will be read depending if they are being asked to do so.


####### NOT YET
### F1  : By Region Mask (defined in mospat_inc_filters.py) and given shapefile
c_RegionMask = []
c_RegionMaskAlias = []

####### 
### F2  : By region provided by user
## i_SqrRegions  =  List of regions. 
##       region  =  [Lat_North, Lat_South, Lon_West, Lon_East]
i_SqrRegions = []
c_SqrRegionsAlias = []


#######
### F4 : By Chilean region 
## Regions are in numbers: 
       # Arica y Pari :     15 
       # Tarapaca     :      1 
       # Antofagasta  :      2 
       # Atacama      :      3 
       # Coquimbo     :      4 
       # Valparaiso   :      5 
       # Metropolitana:     13 
       # O'Higgins    :      6 
       # Maule        :      7 
       # Bio-Bio      :      8 
       # Araucania    :      9 
       # Los Rios     :     14 
       # Los Lagos    :     10 
       # Aysen        :     11 
       # Magallanes   :     12     
## c_Regiones  =  List of groups (list). Each group can be constituted by multiple regiones and they are grouped because it makes sense to study the components together.
c_Regions = []
c_RegionsAlias = []

####### 
### F5 : By Chilean comuna. See config/comunas_codes.csv  
## c_Comunas  =  List of groups (list). Each group can be constituted by multiple comunas and they are grouped because it makes sense to study the components together.
c_Comunas = []
c_ComunasAlias = []


####### 
### F6 : By Station's Elevation or terrain elevation
## i_Elevation = List of tuples. Tuples have the information like [height_min, height_max]
i_Elevation = []
 


################################################################
###############             Filters             ################
##         =  0 Reads everything
##         =  1 Filter by region mask
##         =  2 Filter by square region
##         =  3 Filter by a list of stations' ID
##         =  4 Filter by Chilean region (number)
##         =  5 Filter by Chilean comuna (comuna code)
##         =  6 Filter by terrain elevation (Works for stations only for now)
######  Combinations as concatenated numbers. Examples:  
##  =  13  combination of 1 + 3
##  =  245 combination of 2 + 4 + 5

###### Model Filters.
## Specifies which filters are turned on for models
i_ModelFilters = 0


###### Station Filters. One additional option
## Specifies which filters are turned on for stations
i_StationFilters = 0     

## F3: Filter Number 3: By list of stations
## c_Stations_Obs  =  List of lists, one for each network read. 
##    Each network list is a list of groups.
##        Each group is a list of stations' ID that want to be studied together.
c_Stations_Obs = []
c_Stations_Obs_Aliases = []




################################################################
################### Statistics
## Turns on or off whether to create statistics
## 1 : Statistics per stations and as a whole
## 2 : Statistics per filters and by filters type and all together
## 12: All of the above
i_Stats = 0
i_Stats_freq = []



################################################################
################################################################
###############              PLOTS              ################
## TODO: Explanation of O_1_1 

################################################################
################### Cross sections        
# Cross sections are vertical views for 3D model data. These are defined by two
# (lat, lon) decimal coordinates, usually inside the model grid. There are three
# supported methods to create this view:
# - 'closest_point'
# - 'two_points_interpolation'
# - 'four_points_interpolation' (default)
# Go to 'src/manip/cross_section.py' to learn about the differences or
# run example notebook 'examples/planes_on_model_data/vertical_planes.ipynb'.
#
# Configuration example:
#
# t_cross_sections_config  =  [
#        {'name': 'My cross section',
#         'path': [(lat1, lon1), (lat2, lon2), (lat3, lon3)],
#         'method': 'closest_point'
#         'models': [0, 1]},
#        {'name': 'Another cross section',
#         'path': [(lat3, lon3), (lat4, lon4)]},
#        {'path': [(lat6, lon6), (lat7, lon7)]}]
#
# 'name' is used to identify cross section on plots. If isn't provided it will
# be named 'CROSS_SECTION_<INDEX_OF_CROSS_SECTION>'
# 'path' is a collection of coordinates (max 2 for this version of NeoMOSPAT)
# 'method' is the method used to create the cross section. If isn't provided,
# then 'four_points_interpolation' is used.
# 'models': List of indexes. Indicates on which models the system should get the
# cross section. If isn't defined applies for all models in c_ModelNames.
#
# Cross sections are calculated for all 3D variables found in the selected
# models.
t_cross_sections_config = [{
        "path": []
    }]

# Possible values ['pcolormesh'] or ['pcolormesh', 'scatter']
t_cross_sections_representation = ['pcolormesh']


################################################################
################### Scatter
## 0: None
## 1: Standard
## 2: Density
## 3: Linfit
## 4: Linfit + Density
i_SCT = 0

### Time frequency of data to use for plots
c_SCT_freq = []

### Defines which variable combination to plot
## To Plot every possible combination use 
##        [["O_*_*", "O_*_*"], ["O_*_*", " M_*_*"]]
##     Note: This will create a stupidly amount of images.
c_SCT_Plots = []


i_SCT_diag=1
## Percentile of data to use. 100 means all data, 98 considers up to 98 percent of data.
i_SCT_perc = [100, 98] 
## Plot TS for individual stations
b_SCT_plotstations=False
## Plot TS for filters
b_SCT_plotfilters=True


################################################################
################### Time Series
## Can plot up to two variables and as many models as given.
## 0: None
## 1: Standard Time Series
i_TS = 0    
### Defines which variable combination to plot
c_TS_Plots = []
## First value: Frequency of data to use
## Second value: How often to do the plots. Can use same values as Time freq of whole MOSPAT
##      H,D: makes   daily plots with  hourly data
##      D,W: makes  weekly plots with   daily data 
c_TS_freq = []
## Rolling mean window for corresponding c_TS_freq
i_TS_RM=[]
## Defines if all range of points is shown or if only the 99.5 percentile. (Made to remove outliers in case they are too high and prevent from looking at the general behavior)
i_TS_ShowAll = 0   
## Defines whether to force same range for variables asked to plot
i_TS_SharedYRange = 1
## Wether to erase model data if obs data is not there
i_TS_SharedXRange = 1
## Plot TS for individual stations
b_TS_plotstations=False
## Plot TS for filters
b_TS_plotfilters=True
## Highlight first model
b_TS_hlmodel=False
## Indexes of models to highlight
c_TS_hlmodel=[0]
## Highlight first model
b_TS_hlobs=True
c_TS_hlobs=[0]

################################################################
################### Time Series Summary
## 0: None
## 1: Raw data
## 2: Raw data + std
i_TS_Sum = 0
### Defines which variable combination to plot
c_TS_Sum_Plots = []
## First value: Type of cycle to create
## Second value: How often to do the plots. Can use same values as Time freq of whole MOSPAT
##           H,DC,M: makes  daily cycle  for every month that was read with hourly data
##       D,AC,Whole: makes  anual cycle  for data from the whole period that was read with daily data points
c_TS_Sum_freq = []
## Rolling mean window for corresponding c_TS_freq
i_TS_Sum_RM=[]
## Defines if all range of points is shown or if only the 99.5 percentile. (Made to remove outliers in case they are too high and prevent from looking at the general behavior)
i_TS_Sum_ShowAll = 0
## Defines whether to force same range for variables asked to plot
i_TS_Sum_SharedYRange = 1
## Wether to erase model data if obs data is not there
i_TS_Sum_SharedXRange = 1

################################################################
################### Time Series Statistics
## Can plot up to two variables and as many models as given.
## 0: None
## 1: Standard Time Series
i_TS_stats = 0    
### Defines which variable combination to plot
c_TS_stats_Plots = []
## First value: Frequency of data to use
## Second value: How often to do the plots. Can use same values as Time freq of whole MOSPAT
##      H,D: makes   daily plots with  hourly data
##      D,W: makes  weekly plots with   daily data 
c_TS_stats_freq = []
## Defines if all range of points is shown or if only the 99.5 percentile. (Made to remove outliers in case they are too high and prevent from looking at the general behavior)
i_TS_stats_ShowAll = 0   
## Defines whether to force same range for variables asked to plot
i_TS_stats_SharedYRange = 1
## Wether to erase model data if obs data is not there
i_TS_stats_SharedXRange = 1
## Plot TS for individual stations
b_TS_stats_plotstations=True
## Plot TS for filters
b_TS_stats_plotfilters=True



################################################################
################### Maps
### Only available for models for now
## 0: None
## 1: Map, model grid visible
## 2: Map, values interpolated
## 3: Map, countur fill
i_Map = 0

### Defines which variable combination to plot
## For every possible combination use 
##        ["M_*_*"]
##     Note: This will create a stupidly amount of images.
c_Map_Plots = []

### How often to do the maps. Takes values defined in c_TimeMeanRes
c_Map_freq = []

### Type of plot to do
## 1: mean
## 2: max
## 3: std
c_Map_type = 1

### Value range for colorbar. It is one per plot in c_Map_Plots
## Written as list of [vmin, vmax]
## They can take values such as:
##    Number    : Equal that number
##    PCTL_all  : Percentile value for the whole data of timeframe asked
##    PCTL_freq : Percentile value for the data of the map in question
## Examples:     ["PCTL_all:0", "PCTL_all:98"], ["PCTL_freq:0", "PCTL_freq:98"], [0, 50]
c_Map_ColBar = []

### Add points of network in the map (You can ask all options)
## 0: No network points in the maps
## 1: Points of networks
## 2: Points of networks with values of stations for equivalent variables
## 3: Points of filters
## 4: Points of filters with average value for equivalent variables
c_Map_AddNet = 0

### Add Shapefile boundaries 
## Aesthetics of how they are plotted can be changed in  src.plot.aux_plots class ShapefileMapConfig
## 0: None
## 1: Ciudades 
## 2: Comunas
## 3: Provincias
## 4: Regiones
## 5: Chile
## 6: Topo_250m
## 7: Squares
## 8: Points
i_Map_Shapefiles = 0

### Add background (Data is plotted with 50% alpha)
## Aesthetics of how they are plotted can be changed in  src.plot.aux_plots class ShapefileMapConfig
## 0: No background
## 1: Satelite       (ctx.providers.Esri.WorldImagery)
## 2: Google Maps    (ctx.providers.OpenStreetMap.Mapnik)
## 3: Simple Terrain (ctx.providers.Stamen.TerrainBackground)
i_Map_Background = 0


## Squares to plot on N7 i_Map_Shapefiles
f_Map_Squares=[]

## Points to plot on N8 i_Map_Shapefiles
f_Map_Points=[]
c_Map_mosaic=[]

################################################################
################### Maps Summary
## 0: None
## 1: Map, model grid visible
## 2: Map, values interpolated
## 3: Map, countur fill  
#########
## To Plot every possible combination use 
##     ["M_*_*"]
##     Note: Only available for models for now
i_Map_Sum = 0
c_Map_Sum_Plots = []
c_Map_Sum_ColBar = []
## How often to do the maps
## First value: Frequency of data to use 
## Second value: How often to do the plots. Can use same values as Time freq of whole MOSPAT
##            H,D: makes    daily maps  with  average from hourly data
##            D,M: makes  monthly maps  with  average from daily data 
##        H,Whole: makes  one map of whole period with the average from  the  hourly data 
##           H,DC,M: makes  daily cycle  for every month that was read with hourly data
##       D,AC,Whole: makes  anual cycle  for data from the whole period that was read with daily data points
c_Map_Sum_freq = []
## Value range for colorbar. It is one per plot in c_Map_Plots
## Written as list of [vmin, vmax]
## They can take values such as:
##      Number    : Equal that number
##      PCTL_all  : Percentile value for the whole data of timeframe asked
##      PCTL_freq : Percentile value for the data of the map in question
#c_Map_ColBar = [["PCTL_all:0", "PCTL_all:98"], ["PCTL_freq:0", "PCTL_freq:98"]]
######### Add points of network in the map (You can ask all options)
## 0: No network points in the maps
## 1: Points of networks
## 2: Points of networks with values of stations for equivalent variables
c_Map_Sum_AddNet = 0
######### Add Shapefile boundaries 
## 0: None
## 1: Ciudades 
## 2: Comunas
## 3: Provincias
## 4: Regiones
## 5: Chile
## 6: Topo_250m
## 7: Squares
i_Map_Sum_Shapefiles = 4
f_Map_Sum_Squares = []

## 0: None
## 1: Map, model grid visible
## 2: Map, values interpolated
## 3: Map, countur fill
i_Map_Stats=0
c_Map_Stats_Plots=[]
c_Map_Stats_ColBar=[]
c_Map_Stats_freq=[]
c_Map_Stats_AddNet=0
## 0 : All
## 1 : Bias / NMBias / STD
## 2 : FGE / RMSE 
## 3 : R / FGE
c_Map_Stats_type=0
c_Map_Stats_bg=0





################################################################
################################################################
###############              Manip              ################



################################################################
################### Model data to shapefile polygons
#### Model data onto shapefile/netcdf in c_RMap_dir (keeps geometry of shapefile)
## 0: None 
## 1: Aggregates variables as  average  in geometry
## 2: Aggregates variables as  sum  in geometry
i_RMap = 0

## Full path of shapfile to read
c_RMap_dir = ""
## Projection of shp. Format "epsg:<projnumber>". 
## Examples of typical <projnumber>
##     4326   =   lat lon (in degrees) 
##     3857   =   Pseudo-Mercator (typical of maps, in meters)
c_RMap_crs = "epsg:4326"
c_Rmap_ID="ID"
## How to save the new data
## 1: As shapefile
## 2: As network equivalent 
c_RMap_savetype = 1
## Directory where to save mapped information
c_RMap_savedir = ""
## Data to map to shapefile. It can be models or observations
c_RMap_data = []
## Data frequency to save. Specified frequency 
## Example: "D" saves remapped daily data of model or observation 
c_RMap_freq = []



c_TS_mosaic = []
c_Map_mosaic = []
################################################################
################### Bias Correction
i_BC = 0
## Data frequency to correct
c_BC_freq = []
## Window to obtain statistical parameters
i_BC_window = []
## Radius to use
i_BC_Rads = []
## Observational vars used to correct. Only "O_i_j"
c_BC_Vars = []
## Models to correct. Only "M_i"
c_BC_Models = []
## New names for Model
c_BC_NewNames = ["%s_BC_R%s"%(m.replace(" ", "_"), "-".join([str(r) for r in i_BC_Rads])) for m in c_ModelAlias]
## Dir where to save corrected data
c_BC_Dir = "."




###################################################################
######## Write 
## 0 : Nothing is saved
## 1 : Models are saved
## 2 : Observations are saved
## 3 : Model at observations locations are saved
i_Write = 0
## How to save network files
## 1 : Save individual stations data
## 2 : Save filter data 
i_Write_net = 12
