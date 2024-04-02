# -*- coding: utf-8 -*-
from datetime import datetime
from os import listdir
from os.path import isfile, join

import geopandas as gpd
import numpy as np
import pandas as pd

from shapely.geometry import Point

import src.IncludeFile as IncF
import src.common as aux
from src.manip import select_stations as ss
from src.plot.ticker_formatters import LatDMS, LonDMS
pd.options.mode.chained_assignment = None 


np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 
#import warnings
#warnings.simplefilter("error", np.VisibleDeprecationWarning)


######################################################
def get_freqs_info(ID, ggroup, directories_db=None):
    ##################
    ## Obtain minimum frquency that was asked for
    ## Get network without the number
    group="_".join(ggroup.split("_")[:-1])

    ## Path to database
    datapath=directories_db["OBSERVATIONS"][group] 

    ## Numeric value of frequencies they want
    val_freq=[aux.return_value_freq(freq)[1] for freq in IncF.c_TimeMeanRes]
    ## Sorted list of frequencies they want
    freq_sort=[IncF.c_TimeMeanRes[val_freq.index(val)] for val in sorted(val_freq)]
    ## Minimum of frequency they want
    freq_min=freq_sort[0]


    ############
    ## Files in directory that have the ID of station
    onlyfiles = [f for f in sorted(listdir(datapath)) if isfile(join(datapath, f)) and "ID-%s"%(ID) in f ]
    freqs_file=list(set([file.split("_")[-1].split(".")[0] for file in onlyfiles]))

    return freq_min, len(freqs_file)


######################################################
### SINCA  file name
def pname(ID, ggroup, types, freq_min, directories_db=None, freqs_used=[]):
    ## Get network without the number
    group="_".join(ggroup.split("_")[:-1])

    ## Path to database
    datapath=directories_db["OBSERVATIONS"][group] 

    ## Files in directory that have the ID of station
    onlyfiles = [f for f in sorted(listdir(datapath)) if isfile(join(datapath, f)) and "ID-%s--%s"%(ID, types) in f ]
    onlyfiles = [f for f in onlyfiles if aux.return_value_freq(f.split(".")[-2].split("_")[-1])[0] not in freqs_used]

    freqs=[file.split("_")[-1].split(".")[0] for file in onlyfiles]
    freqs_file=sorted(set(freqs), key=freqs.index) 
    freqs_file=[aux.return_value_freq(f)[0] for f in freqs_file]
    freqs_file=[f for f in freqs_file if f not in freqs_used]    


    if freq_min in freqs_file:
        file=onlyfiles[freqs_file.index(freq_min)]
        freqs_used.append(freq_min)
        #print("1 File: ", join(datapath, file))
        return join(datapath, file)        


    if len(onlyfiles)==1:
        #print("1 File: ", join(datapath, onlyfiles[0]))
        freqs_used.append(freqs_file[0])
        return join(datapath, onlyfiles[0])
    elif len(onlyfiles)==0:
        #print("4 File: ", "No file")
        return None
    ## If desired frequency is in files upload it
    ## NEed to add that
    else:
        ## Minimum of the freqs available
        val_freq=[aux.return_value_freq(freq)[1] for freq in freqs_file]
        freq_min=freqs_file[val_freq.index(min(val_freq)) ]
        ##
        file=onlyfiles[freqs_file.index(freq_min)]
        freqs_used.append(freq_min)
        return join(datapath, file)

        
######################################################
### Imports individual CSV with all of its pre-processsing
def import_csv(name, timeframe):
    if name is None:
        return pd.DataFrame()
    ## Read csv
    try:
        df=pd.read_csv(name, dtype={'Fecha': 'str'})
        ## Converts to date time object
        df['Fecha']=pd.to_datetime(df['Fecha'])
    except:
        df=pd.read_csv(name, dtype={IncF.i_date_column: 'str'})
        ## Converts to date time object
        df['Fecha']=pd.to_datetime(df[IncF.i_date_column])        

    ## Extracts given timefrme
    df=df.loc[(df["Fecha"]>=timeframe[0]) & (df["Fecha"]<=timeframe[1])]
    ## Removes empty columns
    df=df.dropna(axis=1, how='all')

    ## Convert everything but date to float
    ## TODO: make it so that we can check range
    cols=[i for i in df.columns if i not in ["Fecha", "WINDU", "WINDV", "WDIR"]]
    for col in cols:
        df[col]=pd.to_numeric(df[col], errors="coerce",  downcast='float')
        df.loc[df[col]<0, col]=np.nan


    ## Remove days that dont have any measurement
    #df=df.dropna(how="all", axis=1).dropna(thresh=2)
    df.index=df["Fecha"]
    df.index.name="time"
    df.drop(columns=["Fecha"], inplace=True)

    return df.sort_index()


######################################################  
def read_data(ID, group, freq_min, variables, timeframe, directories_db=None, freqs_used=[]):
    ## Import individual CSVs
    try:
        #breakpoint()
        Cal=import_csv(pname(ID, group, "Cal", freq_min, directories_db=directories_db, freqs_used=freqs_used), timeframe)    
        BoolCal=True
        if Cal.empty:
            BoolCal=False        
    except:
        BoolCal=False

    try:
        Met=import_csv(pname(ID, group, "Met", freq_min, directories_db=directories_db, freqs_used=freqs_used), timeframe)
        BoolMet=True
        if Met.empty:
            BoolMet=False
    except:
        BoolMet=False


    ## DataFrame with all the Cal and Met info
    if BoolMet and BoolCal:
        All=pd.merge(Cal, Met, left_index=True, right_index=True, how="outer")
    elif BoolCal:
        All=Cal
    elif BoolMet:
        All=Met
    else:
        return None, None

    ########
    ### Keep only columns of interest
    ## Find the names of the variables 
    dic_vars=ss.find_equiv_variables(All.columns, variables)
    dic_vars_rev={v: k for k, v in dic_vars.items()}
    cols=list(dic_vars.values())


    units={c.split("_")[0].split("--")[0]:c.split("_")[-1].split("--")[0] for c in All.columns}
    #breakpoint()
    units={dic_vars_rev[k]:v for k, v in units.items() if k in cols}

    ## Only keep variables of interest  
    All=All[cols]

    #if All.shape[1]!=len(variables):
    #    return None, None


    #breakpoint()
    if len(units)==0:
        #units=read_units(group)
        ncols=[dic_vars_rev[c] for c in cols]
        units={k:v for k,v in aux.read_units("_".join(group.split("_")[:-1])).items() if k in ncols} 
    ## If there are no columns
    if not cols:
        return None, None

    return All, units
    ## Return None when there are missing variables
    try:
        All.columns=[dic_vars_rev[c] for c in cols]
        All.meta.units=units
        if "*" in variables:
            return All, units
    except:
        ## Deals with:
        ## When there is no longer a height for the station within the range
        ## And if a variable completely dissapeared in the time range
        #print("      ...Missing variable for station", ID, " for timeframe:  ",timeframe)
        return None, None

    return All, units


######################################################
### Import all data
def import_all(ID, group, variables, vars_all,  timeframe, directories_db=None):
    freq_min, nfreqs = get_freqs_info(ID, group, directories_db=directories_db)
    freqs_used=[]
    c=0
    All=None
    while c!=nfreqs:
        c+=1
        All, units=read_data(ID, group, freq_min, vars_all, timeframe, directories_db=directories_db, freqs_used=freqs_used)

        if All is not None:
            break

    if All is None:
        return None

    ## Remove height from name
    All.columns=[col.split("--H")[0] for col in All.columns]
    #All.columns=[dic_vars_rev[c] for c in cols] 

    ######
    ## U:Zonal, V:Meridional
    ## Create Zonal and meridial direction
    if "D-U" in variables or "D-V" in variables:
        ## Make the math
        All["D-U"], All["D-V"] = np.sin(All["WDIR"]*np.pi/180), np.cos(All["WDIR"]*np.pi/180)
        units["D-U"]=units["WDIR"]
        units["D-V"]=units["WDIR"]
    ## Create Zonal and meridial speed
    if "WINDU" in variables or "WINDV" in variables:
        ## Make the math
        All["WINDU"]=-np.sin(np.deg2rad(All["WDIR"]).astype(np.float64))*All["WSPD"]
        All["WINDU"]=pd.to_numeric(All["WINDU"], errors="coerce",  downcast='float')
        All["WINDV"]=-np.cos(np.deg2rad(All["WDIR"]).astype(np.float64))*All["WSPD"]
        All["WINDV"]=pd.to_numeric(All["WINDV"], errors="coerce",  downcast='float')
        units["WINDU"]=units["WSPD"]
        units["WINDV"]=units["WSPD"]

    ## Create coarse variable
    if "PMC" in variables:
        pm10=pd.to_numeric(All["PM10"], errors="coerce")
        pm25=pd.to_numeric(All["PM25"], errors="coerce")
        pm10=pm10.mask(pm10<=0)
        pm25=pm25.mask(pm25<=0)
        All["PMC"]=pm10-pm25
        All["PMC"]=All["PMC"].mask(All['PMC']<=0)
        All["PMC"]=pd.to_numeric(All["PMC"], errors="coerce",  downcast='float')
        units["PMC"]=units["PM10"]

    ## Keep requested vars
    All=All[[v for v in variables if v in All.columns]]
    All.meta.units=units

    ## Return Dataframe
    return All


######################################################
### Imports to wanted format
def read_to_format(c_ObsNetwork_local, t_Models_Info):
    ## Obtains filters to apply
    list_filters=IncF.i_StationFilters
    info_dfs=[]



    #########################
    ## Go through each observations network
    n_groups=len(IncF.c_ObsNetwork)
    data={}
    networks_used=[]
    t_filters={}
    for index in range(n_groups):
        units={}
        group=c_ObsNetwork_local[index]
        varia=IncF.c_ObsVars[index]

        ## Directories
        directories_db = aux.read_yaml("data")

        print("...Filtering stations of network: ", group)
        t_filters_group, df = ss.filter_stations(t_Models_Info, group, varia, directories_db=directories_db)
        t_filters = {**t_filters, **t_filters_group}

        
        print("   ...Variables:  ", " , ".join(varia))
        print("   ...Number of stations found:  ", len(list(df["Nombre"].values)))
        ###########################
        ###########################
        

        
        ## If no stations then skip the import part
        if df.empty:
            print("   ...No stations to import")
            IncF.c_ObsNetwork[index]=None
            continue

        print("   ...Importing data of: ", group)
        ################
        if len(IncF.d_Start_Date)==len(IncF.d_Last_Date):
            data_red=[]
            units={}
            variables=ss.get_all_vars_requested(varia)
            for ID in df["ID"]:
                if "f" == ID[0] and "-" == ID[2]:
                    continue                

                name_station=df[df["ID"]==ID]["Nombre"].values[0]
                print(f"      ...Importing  {ID}: {name_station}")# , end="\n         ")
                ID_data=[]
                ## Retrieve data for all timeframes in IncludeFile
                for index_timeframe in range(len(IncF.d_Start_Date)):
                    date_ini = IncF.d_Start_Date[index_timeframe]
                    date_fin = IncF.d_Last_Date[index_timeframe]

                    ## Import data
                    info=import_all(ID, group, varia, variables, [date_ini, date_fin], directories_db=directories_db)
                    ## Checks if data was found
                    if info is None or info.isnull().values.all() or info.empty:
                        continue
                    ## Append retrieved data
                    else:
                        ID_data.append(info)
                        units.update(info.meta.units)

                ## Check if there is no data for all timeframes asked in IncludeFile
                if len(ID_data)==0:

                    ## Removes station from list if no data was found
                    idx = df[df['ID']==ID]
                    df.drop(idx.index, axis=0, inplace=True)
                    print("         ...No data for timeframe or variables requested")
                    continue
                else:
                    data_red.append(pd.concat(ID_data))
                    #print("         ...Done")

            ## Check if no station for network was read
            if len(data_red)==0:
                print("   ...No stations were imported for network  ", group)
                IncF.c_ObsNetwork[index]=None
                continue
            else:
                df_imported=df[~df["ID"].str.contains('^f[0-9]-+')]
                print("      ...Stations Imported (%s): "%len(df_imported["ID"].values), " , ".join(df_imported["ID"].values))
                networks_used.append(group)
                mi_data=pd.concat(data_red, axis=1, keys=df["ID"]).swaplevel(axis=1).sort_index(axis=1)
                data[group]=mi_data
                mi_data.meta.units=units
                info_dfs.append(df.reset_index(drop=True))
                ## Remove stations without data from filters
                t_filters[group]={f:{ k:[s for s in v if s in df["ID"].values] for k, v in t_filters[group][f].items() } for f in t_filters[group]}


    ###########
    ## Multi-Indexed Dataframe with all the info of the groups
    #local_c_ObsNetwork=["SINCA", "SINCA_2"]
    ## Condition if repeated add sufix
    if len(info_dfs)!=0:
        mi_info=pd.concat(info_dfs, axis=0, keys=networks_used)
        geometry=[Point(xy) for xy in zip(mi_info.Longitud, mi_info.Latitud)]
        mi_info=gpd.GeoDataFrame(mi_info, crs='epsg:4326', geometry=geometry)
        mi_info.crs='epsg:4326'
        mi_info["location"]=mi_info.apply(lambda  df: "%s %s"%(LatDMS(df.Latitud).to_string(resolution='minutes'), LonDMS(df.Longitud).to_string(resolution='minutes')), axis=1)
        return [t_filters, mi_info, data]
    else:
        ## TODO: Add missing columns
        return {}, pd.DataFrame(columns=["ID", "lat", "lon"]), {}


######################################################
def main(t_Models_Info):
    if not IncF.c_ObsNetwork:
        #t_Stations_Filters, t_Stations_Info, t_ObsStationData
        print("...No observations were asked to be read. Skipping this step.")
        return {}, pd.DataFrame(), {f:pd.DataFrame() for f in IncF.c_TimeMeanRes}

    ###############
    ## Rename repeated observations networks
    ## Adds _1 or _2 to names
    c_ObsNetwork_local = aux.rename_groups(IncF.c_ObsNetwork)

    ##########
    ## Get observations data
    t_filters, stations_info, data = read_to_format(c_ObsNetwork_local, t_Models_Info)



    ###############
    ## Frequencies that were asked for
    vals_freqs = {aux.return_value_freq(f)[1]: aux.return_value_freq(f)[0] for f in IncF.c_TimeMeanRes}
    wanted_freqs = list(vals_freqs.keys())

    ###############
    ## Initialize into dics
    keys=list(vals_freqs.values())
    obs_data=dict.fromkeys(keys)
    for time_res in keys:
        obs_data[time_res]=dict({})



    ################
    ## Check if frequency available is the lowest one asked or even lower
    print("")
    print("...Creating required time frequencies", ", ".join(IncF.c_TimeMeanRes))
    for red in data:
        print("  ...For network   ", red)
        ################
        ## Data for red
        df = data.get(red, None)
        if df is None:
            continue

        ## Get units of Network
        group="_".join(red.split("_")[:-1])

        if len({k:v for k,v in df.meta.units.items() if k!=v})!=0:
            units=df.meta.units
        else:
            units=aux.read_units(group)


        ## Find frequency of original data
        freq=aux.find_freq(df.dropna(thresh=1).index)



        ##########
        ## Figure out index of first frequency that can be created or exists and save its index
        for count, asked_freq in enumerate(wanted_freqs):
            if freq<=asked_freq:
                index=count
                break

        ###########
        ## Save data to dictionaries
        ## Use minimum index to rearrange the data 
        min_freq=vals_freqs[wanted_freqs[index]]
        offset=df.index[0].replace(hour=0)-min(IncF.d_Start_Date).replace(hour=0)
        offset=pd.Timedelta(days=0) 

        obs_data[min_freq][red] = aux.change_freq(df, aux.return_freq_value(freq), aux.return_freq_value(wanted_freqs[index]), offset=offset.floor("D") )
        obs_data[min_freq][red].meta.units = units
        obs_data[min_freq][red] = obs_data[min_freq][red].reindex(sorted(obs_data[min_freq][red].columns), axis=1)

        ## Add info of filters
        for f in t_filters[red]:
            groups = t_filters[red][f]
            for al, group in groups.items():
                ## Stations data
                df = obs_data[min_freq][red]
                ## Filters data
                data_group = obs_data[min_freq][red].iloc[:, obs_data[min_freq][red].columns.get_level_values(1).isin(group)].mean(level=0, axis=1)
                ## Add filters data to stations data
                df[[(c, f"{f}-{al}") for c in data_group.columns]] = data_group
                obs_data[min_freq][red] = df

        obs_data[min_freq][red].meta.units=units

        ## Get other frequencies
        for re_freq in wanted_freqs[index+1:]:
            current_freq=vals_freqs[re_freq]
            ## Here there needs to be something procedural that will go H->D->M->Y
            #breakpoint()
            obs_data[current_freq][red]=aux.change_freq(df, aux.return_freq_value(freq), aux.return_freq_value(re_freq))
            obs_data[current_freq][red].meta.units=units

    return t_filters, stations_info, obs_data
