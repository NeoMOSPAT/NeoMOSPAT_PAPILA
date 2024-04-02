# -*- coding: utf-8 -*-
from datetime import datetime
from typing import List, Type, Union, Tuple
from os import listdir
from os.path import isdir

import geopandas as gpd
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import src.IncludeFile as IncF

# ModelType: A model can be any of these data structures
ModelType = Union[Type[xr.Dataset], Type[pd.DataFrame]]

# Possible keys used to refer time arrays
c_time_keys = {'Time', 'time', 'TIME', "times", "Times"}
# Possible keys used to refer latitude, longitude and height
c_lat_keys = {'Lat', 'lat', 'LAT', 'Latitude', 'latitude', 'LATITUDE'}
c_lon_keys = {'Lon', 'lon', 'LON', 'Longitude', 'longitude', 'LONGITUDE'}
c_height_keys = {'Height', 'height', 'HEIGHT'}


#####################################
def nc_reader(c_file: str, c_variables: List[str]) -> Type[xr.Dataset]:
    """
    Reads NetCDF4 files containing model simulations.

    It can retrieve multiple variables from a single file. Since I/O operations
    are costly this function only works with one file. If you want to read multiple
    files call this function multiple times.

    Data is usually spatial (planar or volumetric) with a time dimension.

    Args:
        c_file (str): Full path of NetCDF4 file
        c_variables (List[str]): Variables to be extracted from the file.
        TODO: add arguments for time/spatial filtering

    Returns:
        xarray object that consists in a NetCDF-in-memory like structure.

    Raises:
        FileNotFound.
        IOError: Usually if a variable doesn't exist or is spelled incorrectly.
    """
    c_variables = [c_variables] if isinstance(c_variables, str) else c_variables
    # Open file
    with nc.Dataset(filename=c_file) as ds:
        # NOTE: to open multiple files dataset exist xr.open_mfdataset
        
        # Get variables
        c_vars_set = set(c_variables)
        c_file_vars = list(ds.variables)
        nonexistent_vars = c_vars_set.difference(c_file_vars)
        existent_vars = c_vars_set.difference(nonexistent_vars)

        ## Gets units of variables asked
        units={var:ds.variables[var].units  for var in c_vars_set if var not in nonexistent_vars}

        #if nonexistent_vars:
        #    print("         ...Omitting variables", " ".join(nonexistent_vars),  "not in file.")

        #f_vars_data = {c_var_name:ds.variables[c_var_name][:] for c_var_name in existent_vars}
        #breakpoint()
        #f_vars_data = {c_var_name:(ds.variables[c_var_name].astype(casting="safe")[:] if c_var_name!="optdaero" else ds.variables[c_var_name][:][:, 0, :, :]) for c_var_name in existent_vars}
        f_vars_data = {}
        for var in existent_vars:
            arr=ds.variables[var][:]
            fill_value=arr.fill_value
            arr=np.array(arr)
            arr[np.where(arr==fill_value)]=np.nan
            arr=arr.astype("float32")
            f_vars_data[var]=arr
        # Get time, latitude, longitude and height

        # time key is not the same for all files
        try:
            c_time_key = c_time_keys.intersection(c_file_vars).pop()
            f_time = ds.variables[c_time_key]
            #f_time_units = f_time.units
            #f_calendar = f_time.calenda
            d_time  = nc.num2date(f_time[:], units=f_time.units, calendar=f_time.calendar)
            d_time = pd.to_datetime([datetime(cftime.year,cftime.month,cftime.day,cftime.hour,cftime.minute) for cftime in d_time])
            if sum(n < 0 for n in f_time[:])>0:
                ## Happens when order variables (time, lat, lon) are not available at variables to get the values inside. Why? Idk
                time_len=ds.dimensions[c_time_key].size
                date=c_file.split("_")[-1].replace(".nc", "")
                time_object=datetime.strptime(date, "%Y%m%d00")
                d_time = [time_object + timedelta(hours=n) for n in range(time_len)]

        except:
            ## Happens when order variables (time, lat, lon) are not available at variables to get the values inside. Why? Idk
            c_time_key=c_time_keys.intersection(ds.dimensions).pop()
            time_len=ds.dimensions[c_time_key].size
            date=c_file.split("_")[-1].replace(".nc", "")
            time_object=datetime.strptime(date, "%Y%m%d00")
            d_time = [time_object + timedelta(hours=n) for n in range(time_len)]

            #return None
        # Primero se usa num2date pues las fechas vienen en formato "Days since X".
        # Esto las convierte en objetos compatibles con datetime.datetime
        #d_time  = nc.num2date(f_time[:], units=f_time.units, calendar=f_time.calendar)
        # Luego se convierte a un objeto DateTimeIndex
        #d_time = pd.to_datetime([datetime(cftime.year,cftime.month,cftime.day,cftime.hour,cftime.minute) for cftime in d_time])
        # TODO: Detect lat lon models to only carry a 1-dimenssional array
        # lat and lon have len(shape) == 2,  dims = (lat, lon)
        #breakpoint()
        ### When there is no lat lon info
        try:
            c_lat_key = c_lat_keys.intersection(c_file_vars).pop()
            c_lon_key = c_lon_keys.intersection(c_file_vars).pop()
            f_lat = ds.variables[c_lat_key][:]
            f_lon = ds.variables[c_lon_key][:]
            f_lat = np.float32(f_lat)
            f_lon = np.float32(f_lon)
        except KeyError:
            return None

        ## Fix if netcdf has a 1D latitudes and longitudes
        if f_lat.ndim==1:
            size_lon=f_lon.shape
            f_lat=np.tile(f_lat, [size_lon[0], 1]).T

        if f_lon.ndim==1:
            size_lat=f_lat.shape
            f_lon=np.tile(f_lon, [size_lat[0], 1])



        # height has len(shape) == 4, dims = (time, height, lat, lon)
        c_height_key = c_height_keys.intersection(c_file_vars)
        f_height = ds.variables[c_height_key.pop()][:] if c_height_key else None
        is3D = not (f_height is None)


    # Build new xarray Dataset with the extracted variables. TODO: An improvement
    # can be made for files that hold only one variable or when param "variables"
    # is equal to ds.variables.keys(), then the whole content can be returned.
    # A xarray.Dataset is a collection of xarray.DataArray. Each one containing a
    # variable with declared dimensions and coordinates. All xarray.DataArrays
    # share the coordinates and dimmensions.

    if is3D:
        dims = ('time', 'z', 'y', 'x')
        coords = {
            'time': (('time',), d_time),
            'height': (dims, f_height),
            'lat': (('y', 'x'), f_lat),
            'lon': (('y', 'x'), f_lon)
        }

    else:
        dims = ('time', 'y', 'x')
        coords = {
            'time': (('time',), d_time),
            'lat': (('y', 'x'), f_lat),
            'lon': (('y', 'x'), f_lon)
        }
    t_dataArrays = {}
    for var in existent_vars:
        # TODO: try to optimize by adding dims only once
        # Adds entry to dictionary {var_name: dataArray}
        try:
            t_dataArrays[var] = xr.DataArray(
                data=f_vars_data[var],
                coords=coords,
                dims=dims,
                attrs=None # Units can be added here as a dictionary
            )
        except:
            t_dataArrays[var] = xr.DataArray(
                data=f_vars_data[var][:, 0, :, :],
                coords=coords,
                dims=dims,
                attrs=None # Units can be added here as a dictionary
            )            

    out_ds = xr.Dataset(data_vars=t_dataArrays, attrs={'model_name': c_file.split('/')[-1].split('_')[1], "units":units})
    del f_vars_data
    del t_dataArrays

    #out_ds["lat"]=out_ds["lat"].astype(np.float32)
    #out_ds["lon"]=out_ds["lon"].astype(np.float32)

    return out_ds


#####################################
def model_to_pandas(f_model: Type[xr.Dataset]) -> Tuple[Type[pd.DataFrame], Type[pd.DataFrame]]:
    """
    Transform model data stored in a xarray.Dataset to a MultiIndex Pandas.DataFrame.
    To save space coordinates values are returned as a separate Pandas.DataFrame.

    # NOTE: not GeoPandas.DataFrame yet
    Pandas.DataFrame is useful when you need to extract data when you possess k
    the geometry you need. Usually you get such a geometry when using Shapefiles.
    In NeoMOSPAT this format is useful to work model data alongside observations.j
    When there are large areas of NaN in model simulations this format can help you
    to save RAM, in that case you must make sure the previous Dataset is deleted
    to free memory.

    Args:
        f_model: Type[xarray.Dataset] Model Dataset with at least (time, lan, lon)
        dimenssions.
    Returns:
        MultiIndex Pandas.DataFrame object with index as described above.
        Pandas.DataFrame containing values of coordinates lat, lon and height when applies.

    Raises:
        Nothing

    Example:
        > df, coords = model_to_pandas(f_model=my_model)
        > df.head()
                                                    PM25        PM10
        time                       i_lat i_lon
        datetime(2017-07-01 00:00)   0   0      3.758323    3.922479
        datetime(2017-07-01 01:00)   1   0      2.214582    2.387613
        datetime(2017-07-01 02:00)   2   0      2.594995    2.764024
        datetime(2017-07-01 03:00)   3   0      2.755236    2.961900
        datetime(2017-07-01 04:00)   4   0      4.090765    4.343379

        > coords.head()
            lat         lon         indexes
        0   -37.385757  -71.526520  (0, 0)
        1   -37.116779  -71.523193  (1, 0)
        2   -36.847893  -71.519897  (2, 0)
        3   -36.579121  -71.516632  (3, 0)
        4   -36.310455  -71.513367  (4, 0)
    """
    is3D = len(f_model.dims) == 4
    var_list = list(f_model.variables)
    ## If no data was read
    if len(var_list)==0:
        return pd.DataFrame(), pd.DataFrame()

    var_list.remove('lat')
    var_list.remove('lon')
    # Create a df for model variables and other for the grid
    data_df = f_model[var_list].to_dataframe()
    data_df.drop(columns=['lat', 'lon'], inplace=True)
    # 'height' is time dependant so isn't include here, but in the 'data_df'
    coords_df = f_model[['lat', 'lon']].to_dataframe()

    if is3D:
        # height, lat and lon are z, y and x
        # data_df.index names are (time, x, y, z)
        data_df.index.rename(['i_lon', 'i_lat', 'i_height'], level=[1, 2, 3], inplace=True)
        # this order makes more sense
        data_df = data_df.reorder_levels(['time','i_height', 'i_lat','i_lon'])

    else:
        # lat and lon are y and x
        data_df.index.rename(['i_lon', 'i_lat'], level=[1,2], inplace=True)
        # this order makes more sense
        data_df = data_df.reorder_levels(['time','i_lat','i_lon'])

    # Save indexes as tuples:
    # (0,0), (0,1), ... , (0, nlon-1), ..., (nlat-1, nlon-1) 
    ## Reorders so that it is lat, lon and saves indexes as such
    ## Important because order of index will always be i_lat, i_lon
    coords_df['indexes'] = coords_df.swaplevel(i=0, j=1, axis=0).index.values 
    coords_df.index = list(range(len(coords_df))) 
    coords_df[['i_lat', 'i_lon']] = pd.DataFrame(coords_df['indexes'].tolist())    
    data_df.meta.units=f_model.units

    del f_model
    del var_list

    return data_df, coords_df


#####################################
def main():#t_Stations_Info):
    """
    Reads models as asked in IncludeFile.

    Args:
        None. Information needed is in IncludeFile as
        Models to read      : IncF.c_ModelNames
        Variables to read   : IncF.c_ModelVars
        Frequencies to have : IncF.c_TimeMeanRes
        Filter data as with : IncF.i_ModelFilter (and their corresponding other variables)
    Returns:
        Two dictionaries. One with the grid information as a Geopandas.DataFrame. The other with the data as a MultiIndex Pandas.DataFrame object.

    Raises:
        Nothing

    Example:
        ## Usage
        > t_ModelInfo, t_ModelData = model_to_pandas(f_model=my_model)

        ## Structure
        > t_ModelInfo.keys() = [model1, model2]
        > t_ModelData.keys() = [freq1, freq2]
        > t_ModelData[freq1].keys()==t_ModelData[freq2].keys() = [model1, model2]

        ## Content
        > t_ModelData[freq1][model1].head()
                                                    PM25        PM10
        time                       i_lat i_lon
        datetime(2017-07-01 00:00)   0   0      3.758323    3.922479
                                     1   0      2.214582    2.387613
                                     2   0      2.594995    2.764024
                                     3   0      2.755236    2.961900
                                     4   0      4.090765    4.343379

        > t_ModelInfo[model1].head()
            lat         lon         indexes   geometry
        0   -37.385757  -71.526520  (0, 0)   Point(-71.526520, -37.385757)
        1   -37.116779  -71.523193  (1, 0)   Point(-71.523193, -37.116779)
        2   -36.847893  -71.519897  (2, 0)   Point(-71.519897, -36.847893)
        3   -36.579121  -71.516632  (3, 0)   Point(-71.516632, -36.579121)
        4   -36.310455  -71.513367  (4, 0)   Point(-71.513367, -36.310455)
    """

    import src.common as aux
    from src.manip import select_models
    from src.manip import cross_section

    if len(IncF.c_ModelNames)==0:
        print("...No models were asked to be read.")
        return {}, {}, {}, []

    ####################################################
    ####################################################
    ### Part 1:  Read raw information from netcdf files and filter to store less data
    data = {}
    model_info = {}
    models_x_sections = []
    for index, model in enumerate(IncF.c_ModelNames):
        print("...Reading model", model)
        intervals = []
        ###########
        ## Find all files for model
        string_in_common = f"{IncF.c_ModelTimeFreq}[-_]{model}"
        string_in_common = f".*{model}"
        mod_dir=IncF.c_ModelDir+model
        
        ## TODO: Debiese estar en el checker
        if not isdir(mod_dir):
            print("   ...No directory ", mod_dir, ". Skipping model.")
            continue

        c_file_paths = aux.find_files_with_pattern(parent_dir=mod_dir, pattern=string_in_common + r".*nc")

        d_files = {}
        # NOTE(xZevalx): this only works for hourly files, isn't?
        for c_file_path in c_file_paths:
            #breakpoint()
            preposition=c_file_path.split("/")[-1][:2]
            ## Try with period files
            if preposition=="PP":#IncF.c_ModelTimeFreq=="PP":
                file_date1, file_date2 = c_file_path.split("_")[-1].split(".")[0].split("-")
                file_date1 = datetime.strptime(file_date1, '%Y%m%d')
                file_date2 = datetime.strptime(file_date2, '%Y%m%d')
                day_count=(file_date2-file_date1).days+1
                for file_date in (file_date1 + timedelta(n) for n in range(day_count)):
                    if not file_date in d_files:
                        d_files[file_date] = []
                    d_files[file_date].append(c_file_path) 
            ## Daily files
            if preposition=="DD":#IncF.c_ModelTimeFreq=="DD":
                file_date = aux.get_datetime_from_within(c_file_path, datetime_format='%Y%m%d%H')
                # Multiple files with the same date can exist, usually for 3D case.
                if not file_date in d_files:
                    d_files[file_date] = []
                d_files[file_date].append(c_file_path)




        ## Find keys for initial and end date to read
        dates = sorted(list(d_files.keys()))

        model_dfs = []
        model_3d_dfs = []

        # For every period read a bunch of files (?)
        total_dates_range=[]
        for start_date, end_date in zip(IncF.d_Start_Date, IncF.d_Last_Date):
            min_date = min(dates, key=lambda d: abs(d - start_date))
            max_date = min(dates, key=lambda d: abs(d - end_date))

            ## If model is not of the year asked
            if min_date.year != start_date.year:
                if max_date.year != end_date.year:
                    continue
            intervals.append({"Start": min_date, "End": max_date})

            ###########
            if start_date != min_date:
                print(f"   ...No first date {start_date} in files. Used {min_date} instead")
            if end_date != max_date:
                print(f"   ...No last date {end_date} in files. Used {max_date} instead")


            ###########
            ## Reading netcdf for the required dates
            total_dates_range = total_dates_range + [min_date + timedelta(days=x) for x in range((max_date-min_date).days)] 

        # Create a list of all files to read
        files_for_current_period = []
        for date in sorted(list(set(total_dates_range))):
            try:
                files_for_current_period += d_files[date]
            except:
                print("   ...No data found for ", date)
        if len(files_for_current_period)==0:
             print("   ...No data found for period")
             return {}, {}, {}, {}

        print("   ...Importing variables", " ".join(IncF.c_ModelVars))

        model_mesh_old=pd.DataFrame()
        for c_file in sorted(set(files_for_current_period)):
            print("      ...Reading", c_file)
            try:
                xarray = nc_reader(c_file=c_file, c_variables=IncF.c_ModelVars)
            except:
                continue
                #breakpoint()
                nc_reader(c_file=c_file, c_variables=IncF.c_ModelVars)

            if xarray is None:
                #breakpoint()
                #nc_reader(c_file=c_file, c_variables=IncF.c_ModelVars)
                continue

            # Models as DataFrames
            # NOTE(xZevalx): this should be done only if Dataframes are
            # better/required for the subsequent tasks. I "think"
            # xarray.Dataset is more situable for managing data, as is keeps
            # the spatial form of the data
            model_df, model_mesh = model_to_pandas(f_model=xarray)
            try:
                model_mesh=model_mesh_old if not model_mesh_old.empty else model_mesh
                model_mesh_old=model_mesh.copy()
            except:
                model_mesh_old=model_mesh.copy()


            # Cross sections (only for 3D data)
            if not model_df.empty:
                if '3D' in c_file:
                    models_x_sections += cross_section.main(f_model=xarray)
                    model_dfs.append(model_df)
                else:
                    model_dfs.append(model_df)
        ###########
        ## Merge data from different files
        try:
            ## Data was read for the model
            try:
                units=xarray.units
            ## Specific netcdf file wasnt read
            except AttributeError:
                BoolUnits=False
                if len(model_dfs)!=0:
                    for df in model_dfs:
                        try:
                            units=df.meta.units 
                            BoolUnits=True
                            break
                        except:
                            pass

                if not BoolUnits:
                    units={}

            model_dfs = pd.concat(model_dfs).sort_index()
            model_dfs.meta.units=units
            data[model] = model_dfs
            model_info[model] = gpd.GeoDataFrame(model_mesh, geometry=gpd.points_from_xy(model_mesh.lon, model_mesh.lat), crs="epsg:4326")

            nonexistent_vars= list(set(IncF.c_ModelVars).difference(model_dfs.columns) )
            if len(nonexistent_vars)!=0:
                print("   ...The variables", " ".join(nonexistent_vars), "  were omitted. They are not in the files")
            for index, var in enumerate(IncF.c_ModelVars):
                if var in nonexistent_vars:
                    IncF.c_ModelVars[index]=None

        ## No days were read. Either didnt find files or date specification is wrong. It raises that xarray doesnt exist
        except NameError:
            ## Nothing was read for the model
            IncF.c_ModelNames[index]=None
            print("   ...No files were found. Check dates specified")
        ## Files were found but no requested vars were inside so we received empty datasets.
        except ValueError:
            ## Nothing was read for the model
            IncF.c_ModelNames[index]=None
            print("   ...Requested variables are not in the files. Check model name and variables")
        except Exception as E:
            print("!!!!!!Exception not caught! Check this out!!!!!!")
            print("!!!!!!!!!Exception:   ", E)


    if len(data)==0:
        return {}, {}, {}, {}
    ##############
    ## Filter Models to store less information
    ## TODO: Could be done one step higher to save on stored memory
    print("...Filtering Models")
    t_filters, model_info, data = select_models.filter_data(model_info, data)
    model_info={IncF.c_ModelAlias[IncF.c_ModelNames.index(model)]:model_info[model] for model in model_info}
    t_filters={IncF.c_ModelAlias[IncF.c_ModelNames.index(model)]:t_filters[model] for model in t_filters}


    ####################################################
    ####################################################
    ### Part 2: Create final output structure

    ###############
    ## Check what frequencies were asked for
    vals_freqs = { aux.return_value_freq(f)[1]: aux.return_value_freq(f)[0] \
                   for f in IncF.c_TimeMeanRes }
    wanted_freqs = sorted(list(vals_freqs.keys()))


    ###############
    ## Initialize into dics
    keys = vals_freqs.values()
    model_data = dict.fromkeys(keys)
    for time_res in keys:
        model_data[time_res] = dict()



    ################
    ## Check if frequency available is the lowest one asked or even lower
    print("")    
    print("...Creating required time frequencies ", ", ".join(IncF.c_TimeMeanRes))
    for model in data:
        alias_model=IncF.c_ModelAlias[IncF.c_ModelNames.index(model)]
        print("   ...For model   ", model)
        ###########
        ## Find frequency of original data
        time_array = data[model].sort_index().index.get_level_values(0).unique()
        freq = aux.find_freq(time_array)

        ###########
        ## Figure out index of first frequency that can be created or exists and
        # save its index
        index = None
        for count, asked_freq in enumerate(wanted_freqs):
            if freq <= asked_freq:
                index = count
                break

        ###########
        ## Save data to dictionaries
        ## Use minimum index to rearrange the data
        offset=data[model].index[0][0].replace(hour=0)-min(IncF.d_Start_Date).replace(hour=0)
        offset=pd.Timedelta(days=0) 
        min_freq = vals_freqs[wanted_freqs[index]]
        df_min=aux.change_freq(data[model], aux.return_freq_value(freq), aux.return_freq_value(wanted_freqs[index]), offset=offset.floor("D") )
        df_min.meta.units=units
        model_data[min_freq][alias_model] = df_min

        ## Get other frequencies
        for re_freq in wanted_freqs[index+1:]:
            current_freq = vals_freqs[re_freq]
            ## Here there needs to be something procedural that will go H->D->M->Y
            print("      ...Resampling. Might take a while", end= " ", flush=True)
            df_new_freq=aux.change_freq(data[model], freq, re_freq)
            df_new_freq.meta.units=units
            model_data[current_freq][alias_model]=df_new_freq
            print("...Done", flush=True)

    # Concatenate cross sections. This should be done earlier if this data is
    # incorporated into de resampling process
    models_x_sections = cross_section.concat_cross_sections(models_x_sections)
    IncF.c_ModelNames=IncF.c_ModelAlias
    return t_filters, model_info, model_data, models_x_sections


#####################################
def _overlap(first_inter, second_inter):
    for f,s in ((first_inter, second_inter), (second_inter,first_inter)):
        # will check both ways
        both = True
        for time in (f["Start"], f["End"]):
            if s["Start"] < time < s["End"]:
                both = both and True
            else:
                both = both and False
        if both:
            return True
    else:
        return False
