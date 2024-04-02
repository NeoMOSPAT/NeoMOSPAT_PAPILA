# NeoMOSPAT

New version of MOSPAT (Modelling and Observation System and Analysis Tool) project. Adopts newer technology in the python ecosystem.
    
Allows the analysis and comparison of surface measurements (observations) with that of simulations (models) through time series, scatter plots, statistical analysis and others.  
    
### Current finished features:
- Allows the reading of multiple timeframes, shared between observations and models.
- Reading of `DMC`, `SINCA` and `Explorador (NORTE_RM)` databases for surface variables.
- Reading 2D and 3D models from NetCDF4 files.
- Models can be outputed as `xarray.Dataset` or `Pandas.Dataframe` with a separate dataframe for coordinates labels.
- Obtains the closest grid point of stations to model grid.
- Posses filters that allow to group and only view and perform the functions to the information of interest
- "Cross sections" figures. A map that shows the points used and a 'levels' view.


### Examples and showcase of data structures 
Files in `examples/data` are a small version of their original counterparts. NetCDF4 files (extension `.nc`) has been reduced using `ncks`. Their domains has been reduced to adjust to a bounding box close to north and central zones of Chile using the following coordinates:

    - north-west -17.498, -72.044
    - north-east -17.498, -66.800
    - south-west -34.243, -72.044
    - south-east -34.243, -66.800
    
![models_bbox](images/example_data_bbox.png)

### Directory tree description

- `config`: Defines databases' locations (paths) and other needed "permanent" information.
- `examples`:    
    - `data`: Stores data files used in examples scripts and notebooks.
    - `modelreader`: Examples that use `modelreader.py` API.
- `src`: Folder with source code.
    - `manip`: Holds any kind of data manipulation (reformating, reshaping, merging, etc)
    - `plot`: Codes that generate graphs, maps, plots, etc.
    - `read`: Routines that read data from databases and consolidates them into unified data structures.
    - `write`: Routines that write text files and CSVs. 
    
    
___________________________________________________
##     Include File Options
### General 
```
## Freq in which to store data for models and obs
c_TimeMeanRes = ["DD", 'HH'] 
                                ##  TT: Every minute
                                ##  HH: Hourly
                                ##  3H: Every 3 hours
                                ##  DD: Daily
                                #### T: Minute, H:Hour, D:Day, M:Month, Y:Year
                                #### {number} {letter} == Every {number} {letter} 

## Timeframes desired
c_Start_Date = ['01-01-2017', '01-01-2015'] 
c_Last_Date =  ['30-12-2017', '30-12-2015']

## Images savepath
c_FigureDir = your_full_path

```
### Observations 

NeoMOSPAT is created so that it looks for stations that have all variables requested, if no stations are found then it skips it. As such, one can repeat the name of the network in `c_ObsNetwork` and call it with different variables.

```
## Networks to read. Options: "SINCA", "NORTE_RM", "DMC". 
c_ObsNetwork = [net1, net2, net3]   
## Vars to read. Vars available for each network can be found in config/data_directories.yaml under UNITS
c_ObsVars = [list_vars_net1, list_vars_net2, list_vars_net3]

## Methereological data requested height  
i_height  ## Height wanted for MET data, in meters  
i_height_me  ## Margin of error allowed for height  
b_include_zero=True/False  ## Unkown heights in SINCA are marked as zero. Read those if no closer height was found to be available  
```

### Models

If no model is found then it skips it.

```
## Directory
c_ModelDir = path_where_models_are
## Model name (folder name)
c_ModelNames = [model1, model2]
## Short name that will appear in Figures
c_ModelAlias = [alias_model1, alias_model2]
## Vars all models will read
c_ModelVars  = [var1, var2]
```
### Filters

Filters are meant as a way to reduce the amount of data stored when a certain area of interest in known. It also serves the purpose of creating groups of stations/areas for which the same codes can be applied as a whole

```
#################################
##          = 0 Turned Off: Keeps the whole grid
##  NOTYET  = 1 Filter by region mask (.nc or shapefile)
##          = 2 Filter by square region
##          = 3 Filter by a list of stations' ID
##          = 4 Filter by Chilean region (number)
##          = 5 Filter by Chilean comuna (number)
##  NOTYET  = 6 Filter by terrain elevation or stations elevation
######  Combinations as concatenated numbers. Examples:  
## = 13  combination of 1 + 3
## = 245 combination of 2 + 4 + 5

###### Model Filters.
i_ModelFilters = your_combination_models

###### Station Filters. One additional option
i_StationFilters = your_combination_networks 


#################

### !!!!!! Only exclusive for networks          
### F3: Filter Number 3: By list of stations
c_Stations_Obs=[[list1_stations_net1, list2_stations_net1], [list1_stations_net2, list2_stations_net2]]


### F1  : By Region Mask (defined in mospat_inc_filters.py) and given shapefile
c_RegionMask=[]


### F2  : By region provided by user
### [Lat_North, Lat_South, Lon_West, Lon_East]
i_SqrRegions=[square1, square2]
c_SqrRegionsAlias=["Region1", "Region2"]


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
c_Regions=[[4], [5]]


### F5 : By Chilean comuna. See config/comunas_codes.csv  
c_Comunas=[list_comunas_1, list_comunas_2, list_comunas_3]


### F6 : By Station's Elevation or terrain elevation: range_elev = [elev_min, elev_max]
i_Elevation=[range_elev_1, range_elev_2]
 


```
