___________________________________________________
## Directory tree description

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
##     Instructions to adding new functionalities
### Considerations to have when writing code
- Add leading text to variables to specify what kind of type it is
> t_ : Dictionary or DataFrame  
> i_ : Int or list of Ints  
> c_ : String or list of Strings  
> f_ : Float or list of Floats  

- Add `# -*- coding: utf-8 -*-` as the first line of every `.py`.

- Comment what you are doing throughout the code so that it can be understood by an outsider. 

- Create help message for all functions by adding a comment with `"""` after the def of each function (after first line).

- Use 4 white spaces to indent. Your text editor should be able to transform tabs to 4 white spaces.

### Adding functionality to main
- Write one or multiple stand alone files and save them in the appropriate `read`, `write`, `plot`, and/or `manip` folders. Your code should be able to be called from one package only.

- In `main.py` call the package as
> from folder import name_package  
>  .  
>  .  
>  .  
> name_package.main(** args)  

- Add the appropriate flags in `IncludeFile.py` with a brief explanation of how they change the functionality of your code

- Add text explaining `IncludeFile.py` flags in `README.md` and explain the purpose of your new package

___________________________________________________
##     Adding new databases

If it's a new type of database you'll have to create a new corresponding read routine in the `read` folder. The routine must read it into an appropriate and homologous structure so it can be added to the code and keep the existing functionalities. For 2D variables see data structure of surface observations and for 3D variables see that of models (look at the examples folder).

### Surface observations
If another network needs to be added first you need to create the extraction routine of the data which must have certain conditions:
- Create folder `Data` where the extraction code is. The data and general information about the stations will be saved here.
- Save one file per station with historic information for all variables in folder `Data`.  Filename must be `ID-{IDstation}--{typedata}_{freqdata}.csv` with:
     - `IDstation` identifier of station, usually provided in some way or another by the owner.
     - `typedata` is either `Met` or `Cal` referencing to Metereological data or air quality data.
    - `freqdata` is the frequency of measurement of the data, usually HH (hourly) but can be smaller. Use HH for hourly, DD for daily, 10m for every 10 minutes, 3H for every 3 hours, etc.  
- Data in files `ID-{IDstation}--{typedata}_{freqdata}.csv` must fulfill:
     - Column with datetime must be saved as a string with format `%Y-%m-%d %H:%M:%S` and be the first column in file.
     - All numeric columns must not have strings and must account for errors, i.e. no negative or infinite numbers if it doesn't make sense.
- Create a `Stations.csv` with at least the following columns in `Data` folder:
    - `ID`: identifier of stations
    - `lat`, `lon`: coordinates of stations
    - `cod_comuna`, `Comuna`: code and name of stations' comuna, if in Chile
    - `cod_region`, `Region`: number and name of stations' region, if in Chile
- Add path to `Data` in `config/data_directories.yaml` under `Observations`. Name your new network in the proccess.

If all steps were followed, you should be able to use your network by adding to the `IncludeFile.py` the name used in the last step in `c_ObsNetwork`, as per usual.