# -*- coding: utf-8 -*-
import os
import re

from datetime import datetime
from typing import List

import numpy as np
import pandas as pd
import unicodedata
import yaml

import src.IncludeFile as IncF

def clean_network(network):
    return "_".join(network.split("_")[:-1])


def clean_dic_lists(dic):
    return {k:v for k, v in dic.items() if len(v)!=0}

###############
### Trial of flatten
def flat_trial(d, out=[]):
    """ Trial function to flatten dics of lists and dics of arrays of lists """
    if isinstance(out, list):
        flat(d, out=out)
    elif isinstance(out, np.ndarray):
        flat_array(d, out=out)

###############
### flatten dics of lists
def flat(d, out=[]):
    """ Flattens dictionary of lists. Returns list with all values of dic"""
    if not d:
        return []
    for val in d.values():
        if isinstance(val, dict):
            out=flat(val, out=out)
        else:
            out+=val
    return out

###############
### Flatten dics of arrays of lists
def flat_array(d, out=np.array([])):
    """ Flattens dictionary of arrays of lists. Returns array with all values of arrays"""
    for val in d.values():
        if isinstance(val, dict):
            out=flat(val, out=out)
        else:
            out=np.concatenate((out, np.concatenate(val)), axis=0)
    return out

###############
## Translate names for plots
def find_equiv_var(var, obs, models, databases={"O":[], "M":[]}, variables={"O":[], "M":[]}):
    indexes=var.split("_")[1:]
    ###############
    ## Find variables and networks for observations
    if "O" in var:
        ## Names
        database=rename_groups(IncF.c_ObsNetwork)
        if indexes[0]!="*":
            index=int(indexes[0])-1
            database=[database[index]]

        ## Removes databases that didnt import anything 
        database=[c for c in database if c]
        ## Vars
        if indexes[1]=="*":
            var=[IncF.c_ObsVars[get_index_network(d)] for d in database] 
        else:
            index=int(indexes[1])-1
            var=[[IncF.c_ObsVars[get_index_network(d)][index]] for d in database if database] 
   
        ## Add data so that all requested vars can be saved
        databases["O"]=databases["O"]+database
        variables["O"]=variables["O"]+var

    ###############
    ## Find variables and name of models for models
    elif "M" in var:
        ## Names
        if indexes[0]=="*":
            database=IncF.c_ModelNames
        else:
            index=int(indexes[0])-1
            database=[IncF.c_ModelNames[index]]
        database=[c for c in database if c]

        ## Vars
        if indexes[1]=="*":
            var=[IncF.c_ModelVars for d in database] 
        else:
            index=int(indexes[1])-1
            var=[[IncF.c_ModelVars[index]] for d in database] 

        var=[[d for d in v if d] for v in var]

        ## Add data so that all requested vars can be saved
        databases["M"]=databases["M"]+database
        variables["M"]=variables["M"]+var

 
    return databases, variables


def find_equiv_vars(list_vars):
    """ Get the name of networks, models and vars asked for plotting """
    models = IncF.c_ModelNames
    obs = IncF.c_ObsNetwork
    def merge_vars(databases, variables):
        if len(variables)==0:
            return {}
        sums = {}
        for key, value in zip(databases, variables):
            try:
                sums[key] += value
            except KeyError:
                sums[key] = value
        return sums

    databases={"O":[], "M":[]}
    variables={"O":[], "M":[]}

    for var in list_vars:
        databases, variables=find_equiv_var(var, obs, models, databases=databases, variables=variables)

    #########
    ## Create final structure
    data_O=merge_vars(databases["O"], variables["O"])#dict(zip(databases["O"], variables["O"]))
    data_M=merge_vars(databases["M"], variables["M"])#dict(zip(databases["M"], variables["M"]))
    ## Creates dic
    data={"O":data_O, "M":data_M}
    ## Removes stations and models that werent read
    data={"O":{k:v for k,v in data["O"].items()}, "M":{k:v for k,v in data["M"].items()}} 
    #data={"O":{k:v for k,v in data["O"].items() if k in obs}, "M":{k:v for k,v in data["M"].items() if k in models}}
    ## Removes if there is an empty dic
    data={k:v for k,v in data.items() if v}

    ## TODO
    ## Make it so that repeated variables of a network are erased so that it is there only one
    return data

######################################################
######################################################
##################
def rename_groups(c_ObsNetwork):
    ################
    ## Rename repeated observations networks
    ## Adds _1 or _2 to names
    #breakpoint()
    Net_set=list(set(c_ObsNetwork))
    New_names=[]
    count=[0 for x in range(len(Net_set))]
    for network in c_ObsNetwork:
        index=Net_set.index(network)
        counter=count[index]+1
        count[index]=counter
        if network is None:
            New_names.append(None)
        else:
            New_names.append(network+"_"+str(counter))

    return New_names

##################
## Returns index in c_ObsNetwork of renamed variable group
def get_index_network(group):
    ## Find index of group
    index=[count for count, g in enumerate(IncF.c_ObsNetwork) if g=="_".join(group.split("_")[:-1])]
    #if len(index)==0:
    #    return None
    #else:
    index=index[int(group.split("_")[-1])-1]
    return index

##############################################
def remove_accents(word):
    if not isinstance(word, (str, bytes)):#not isinstance(word, (str, unicode)):
        return word

    word=word.strip().replace("\n", " ")
    try:
        word=unicodedata.normalize(u'NFKD', word).encode('ascii', errors='ignore')
    # TODO: add expection class
    except:
        try:
            word = unicode(word, "utf-8")
            word = unicodedata.normalize(u'NFKD', word).encode('ascii', errors='ignore')
        # TODO: add expection class
        except:
            pass
    return word.decode("utf-8") #str(word,"utf-8")#word.decode("utf-8")

######################################################
######################################################
####### Retrieving data from dataframes
##############

# TODO: we should consider extending Pandas.DataFrame or
# GeoPandas.GeoDataFrame to include all these functions as methods
#
# https://pandas.pydata.org/docs/development/extending.html

def get_stations(df, stations):
    return df.loc[:,df.columns.get_level_values(1).isin(stations)]
##############
def get_info_group(mi_info, group):
    return mi_info.loc[group]

##############
def get_group_data(dic, group):
    return dic[group]

##############
def get_station_data(df_group, station):
    idx = pd.IndexSlice
    try:
        df=df_group.loc[idx[:], idx[:, station]]
        df.columns = df.columns.droplevel(1)
    except KeyError:
        df=pd.DataFrame(columns=list(set(list(df_group.columns.get_level_values(0)))), index=df_group.index)

    return df

##############
def get_var_station_data(df_group, station, var):
    idx = pd.IndexSlice
    df=df_group.loc[idx[:], idx[:, station]]
    df.columns = df.columns.droplevel(1)
    return df[[var]]

##############
def get_var_group_data(df_group, var):
    return df_group[var]


################
def read_units(group):
    directories_db=read_yaml("data")
    return directories_db["UNITS"]

#################
def read_yaml(yaml_dir):
    if yaml_dir=="data":
        with open("config/data_directories.yaml") as f:
            directories_db = yaml.safe_load(f)
    elif yaml_dir=="shp":
        with open("config/shp_directories.yaml") as f:
            directories_db = yaml.safe_load(f)
    elif yaml_dir=="vars":
        with open("config/equiv_vars.yaml") as f:
            directories_db = yaml.safe_load(f)

    return directories_db
######################################################
######################################################
##################
## Adds column of variable name with the same value to all subcolumns in dataframe data
def add_columns(data, values, name):
    ## Stations read
    stations=list(data.columns.levels[1])

    ## Create data frame that will hold the time of day
    new=pd.DataFrame().reindex_like(data).drop(columns=data.columns)
    for station in stations:
        new[name, station]=values
    new.columns = pd.MultiIndex.from_tuples(new.columns)

    ## Concatenate original data and the time of day one
    data=pd.concat([new, data], axis=1)
    return data

##################
## Get duration of events
def get_durations(All, var, group):
    def find_events(All, var, group):
        group_index=get_index_network(group)
        index=IncF.c_vars_events[group_index].index(var)
        e_type=IncF.c_event_type[group_index][index]
        e_value=IncF.f_event_values[group_index][index]
        ##
        if e_type=="std":
            mean=All[var].mean()
            std=All[var].std()
            All=All[All[var]>=mean+e_value*std]
            ##
            name_event="%s >= mean + %s std"%(var, e_value)
        ##
        elif e_type=="value":
            All=All[All[var]>=e_value]
            ##
            name_event="%s >= %s"%(var, e_value)
        ##
        else:
            print("...No appropriate definition of events has been given. Returning original DataFrame")
            return All
        ## Remove values that do not exist for the variable
        All=All.dropna(subset=[var])
        return All, name_event

    All, name_event=find_events(All, var, group)
    ###########
    #### Define group events
    ## Define an hour difference
    hour = pd.Timedelta('1H')
    All["T"]=All["Start"]
    dt = All["T"]
    ## Group consecutive hours
    groups=(dt.diff() != hour).cumsum()
    All["Events"]=groups

    ########
    ## Start, End and duration for each event
    Info=All.groupby(['Events', groups])['T'].agg(['min','max']).reset_index().rename(columns={'min':'Start','max':'End'}).drop(['Events', "T"], axis=1)
    ## Get duration of events
    Info["Duration"]=(Info["End"]-Info["Start"]).map(lambda a: a.total_seconds()/3600.0 +1)

    ## Df that will hold duration per hour
    Info=Info.reindex(Info.index.repeat(Info.Duration)).reset_index(drop=True)
    Info.index=dt.index

    return pd.concat([All.drop(columns=["Start", "Duration", "T"]), Info], axis=1), name_event

######################################################
######################################################
##################
## Find frequency
def find_freq(time_array):
    """
    Gets time frequency as hour units of given time array.

    Args:
        # TODO: check if it also works for numpy arrays with datetimes
        time_array: Pandas.DatetimeIndex list of timestamps.

    Returns:
        Frequency as hour units.
    """
    time_delta = np.diff(time_array).min()
    return time_delta/np.timedelta64(1, 's')/3600

##################
## Input string and outputs a numerical value to frequency
def return_value_freq(freq):
    timefreq={"T":1.0/60.0 ,"H":1.0,  "D":24.0, "M":720.0, "Y":12423.0}
    chars=[char for char in freq]
    numbers=[]
    string=[]
    for char in chars:
        try:
            n=int(char)
            numbers.append(char)
        except:
            string.append(str(char))
    string= list(set(string))
    if len(numbers)!=0:
        number=int("".join(numbers))
        return ["%s%s"%(number, "".join(string)),  number*timefreq[string[0]]]

    else:
        return [string[0], timefreq[string[0]]]

##################
## Input numerical value outputs string that pandas understands of frequencies
def return_freq_value(value):
    # NOTE(xZevalx): I think we should use Pandas.to_timedelta instead
    if isinstance(value, (str, bytes)):
        return value
    timefreq={"T":1.0/60.0 ,"H":1.0,  "D":24.0, "M":720.0, "Y":12423.0}
    swaped={value:key for key, value in timefreq.items()}
    for value_freq in sorted(swaped.keys()):
        ## Find value that is greater
        ## Previous value is the basis of the thing
        if value<value_freq:
            break
        ## Save old value
        old_val=value_freq

    ## Generate string of the frequency
    if (value/old_val).is_integer():
        times=int(value/old_val)
        string=swaped[old_val]

        return "%s%s"%(times, string)
    else:
        print("Frequency value ", value ,"does not match integer values of T(minute), H(hour), D(day), M(month), or Y(year)")

# TODO: this should be in 'manip'
##################
## Will change between frequencies
def change_freq(dff, freq1n, freq2n, percentage=0.75, offset=0):
    ## Same frequency, nothing needs to be done.
    freq1=return_freq_value(freq1n)
    freq2=return_freq_value(freq2n)
    if freq1 == freq2 and offset==pd.Timedelta(days=0) :
        return dff
    ## TODO Probably not good
    if freq1 ==freq2:
        try:
            dff.index=dff.index-offset
        except:
            dff.index=dff.index-pd.Timedelta(days=offset)
        return dff

    ## Figure out how many freq1 are in freq2
    dummy = pd.DataFrame(columns=['a'], index=pd.date_range('2018-01-01', periods=2, freq=freq2))
    dummy = dummy.asfreq(freq1)
    max_dates=len(dummy.index)-1
    if max_dates==0:
        print("   ...Asking to change from a smaller frequency to a bigger one. Can't do.\nData is being returned as is.")
        return dff

    ############
    ## Resample function
    def resample_1D(dff, freq1, freq2):
        #########
        ## Resample
        df=dff.copy()
        df=df.asfreq(freq1)
        ## Separate by the desired frequency and get the mean
        df_mean=df.groupby(pd.Grouper(freq=freq2, offset=-offset)).mean()

        ## Check how many non values there were
        df_count=df.resample(freq2).count()
        ## Check how many values there were. It's treated separetly since there
        # arent the same amount of days in a month for example
        df_len=df.bfill().resample(freq2).count()
        ## Get percentage of information
        df_count=df_count/df_len
        ## Boolean dataframe: True for values that should be nan
        variables=list(df.columns)
        vars_aod=[c for c in variables if c[0] in ["AOD", "AE"]]
        vars_rest=[c for c in variables if c[0] not in ["AOD", "AE"]]
        df_count[vars_rest]=df_count[vars_rest]<=percentage
        df_count[vars_aod]=df_count[vars_aod]<=0
        ## Change True values to np.nan
        df_mean=df_mean.mask(df_count)
        return df_mean
    ## Observations
    # NOTE(xZevalx) so this happens only when observations are read?
    if len(dff.index.names)==1:
        ## Resample
        df_mean=resample_1D(dff, freq1, freq2)
    else:
        # NOTE(xZevalx): esto sería muuuucho más rápido (no sé si más sencillo)
        # con xarray. La API de xarray es muy parecida a la de Pandas:
        # http://xarray.pydata.org/en/stable/time-series.html#resampling-and-grouped-operations
        # En este código tbn falta dar soporte a la data 3D, específicamente a
        # 'height'. La data 3D usa la misma grilla (lat, lon) a lo largo del
        # tiempo pero la altura cambia!. Es decir:
        #   (model_data['height'].loc[t, i, j] == model_data['height'].loc[t+1, i, j]) -> False
        #
        # Por lo tanto, 'height' también debe promediarse, es una variable más.
        # Con la lógica actual esto se lograría haciendo groupby ['i_height', 'i_lat', 'i_lon']
        # Ojo, este tipo de mensajes los dejo pues estoy enfocado en otras
        # partes del sistema y no puedo priorizar hacer refactoring a esta parte
        # o donde sea que haya un NOTE. Así que cualquier refactoring lo hago
        # para que funcione la lógica ya implementada.

        #print("   ...Starting slow part")
        dims = dff.index.names
        t_dim = dff.index.get_level_values('time')

        if 'height' in dff.index.names:
            spatial_indexes_names = ['i_height', 'i_lat', 'i_lon']
        else:
            spatial_indexes_names = ['i_lat', 'i_lon']


        df_groupbyresample=dff.groupby(spatial_indexes_names).resample(freq2, level="time")
        df=df_groupbyresample.mean()
        #df_count=df_groupbyresample.count()
        #mask=df_count/(freq2n/freq1n) >=percentage
        #df_mean=df[mask]i
        #df_mean=df_mean.reorder_levels(dims)
        df_mean=df
        """
        df = dff.droplevel('time').reset_index().set_index(t_dim)
        if 'height' in df.columns:
            spatial_indexes_names = ['i_height', 'i_lat', 'i_lon']
        else:
            spatial_indexes_names = ['i_lat', 'i_lon']
        grouped=df.groupby(spatial_indexes_names)
        df_mean = grouped.apply(resample_1D, freq1, freq2).drop(columns=spatial_indexes_names)
        df_mean = df_mean.reorder_levels(dims).sort_index()
        # NOTE(xZevalx): 'time' is already in the index, is it necessary to add it as a column?
        df_mean['time'] = df_mean.index.get_level_values(0)
        #print("   ...End of slow part")
        """
    return df_mean.dropna(how="all")


#################################################
###### Arrays operations
def argmin(arr, value):
    """
    Finds index or indexes (row, col) of the closest float in the array.

    Args:
        arr: List or ndarray of shape 1 or 2.
        value: float. Value to find closest index

    Returns:
        Index: int or tuple(int, int) if shape is 1 or 2 respectively

    """
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr)
    ndims = len(arr.shape)
    i = np.abs(arr - value).argmin()
    if  ndims == 1:
        return i
    elif ndims == 2:
        col_size = arr.shape[1]
        return int(i/col_size), i%col_size


##########################################################
# Files operations

def find_files_with_pattern(parent_dir='.', pattern=r'.*'):
    """
    Gets all files under 'parent_dir' and subdirs with the given pattern.

    Args:
        parent_dir: str Directory to start the search
        pattern: str Regular expression in python format

    Returns:
        List of files with complete path that matches given pattern
    """
    walker = os.walk(parent_dir)
    # crea lista con todos los archivos y sus paths completos
    matches = []
    for cur_path, _, files in walker:
        matches += [cur_path + '/' + f for f in files if re.match(pattern, f)]

    return matches

###########################################
# Dates and timestamps

def get_datetime_from_within(line, datetime_format="%Y%m%d%H"):
    """
    Gets a date contained in 'line' with the format 'datetime_format'.

    Args:
        line: str String containing a date.
        datetime_format: str Format string supported by datetime.datetime.strptime

    Returns:
        Date contained in the string as datetime object

    See also:
        https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior
    """
    # Notices these regex don't follow strictly the equivalences declared in the
    # dictionary below
    regex_hour = r'(((0|1)[0-9])|(2[0-3]))' # from 00 to 23
    regex_month = r'((0*[1-9])|(1[0-2]))'  # 01 to 12 or 1 to 12
    regex_day = r'((0*[1-9])|([1-2][0-9])|(3[0-1]))'  # 01 to 31 or 1 to 31
    regex_year = r'((19|20)([0-9][0-9]))'  # only years from s xx and s xxi

    format_equivalences = {
        '%Y': regex_year,
        '%m': regex_month,
        '%d': regex_day,
        '%H': regex_hour
    }

    # replace the % expressions with their regex equivalences in datetime_format
    # brute force replacement
    original_pattern = datetime_format
    for key in format_equivalences:
        datetime_format = datetime_format.replace(key, format_equivalences[key])

    match = re.search(pattern=datetime_format, string=line)

    if match:
        string_date = match.group()
        if IncF.c_ModelTimeFreq=="DD":
            string_date = line[:-3][-10:]
        
        return datetime.strptime(string_date, original_pattern)
    else:
        raise RuntimeError(f"Can't find date with format {original_pattern} in {line}")


equiv_vars={"DUST":["PMC", "DUST"], "optdaero":["AOD", "optdaero"], "PMC":["PMC", "DUST"], "AOD":["AOD", "optdaero"]}
def equiv_vars(lista):
    def get_equiv(var):
        try:
            return equiv_vars[var]
        except:
            return [var]

    variables=list(set(lista))
    new_vars=[]
    for var in variables:
        new_vars+=get_equiv(var)

    return list(set(new_vars))




#######################################################
########### Model helpers

def get_model_bbox(f_model) -> List[float]:
    # lon-east, lon-west, lat-south, lat-north
    return list(map(lambda val: float(val), [f_model.lon.max(), f_model.lon.min(), f_model.lat.min(), f_model.lat.max()]))
