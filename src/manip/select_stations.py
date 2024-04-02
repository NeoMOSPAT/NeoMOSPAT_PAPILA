# -*- coding: utf-8 -*-
from os import listdir
from os.path import isfile, join

import numpy as np
import pandas as pd
import geopandas as gpd
import src.IncludeFile as IncF
import src.common as aux


def get_all_vars_requested(variables):
    if "*" in variables:
        return  dict(zip(columns, columns))
                
    ## Transform Zonal variables into its parts
    if "WINDU" in variables or "WINDV" in variables:
        variables = [var for var in variables if var not in ["WINDU", "WINDV"]]+["WSPD", "WDIR"]
    if "D-U" in variables or "D-V" in variables:
        variables = [var for var in variables if var not in ["D-U", "D-V"]]+["WDIR"]
    if "PMC" in variables:
        variables = [var for var in variables if var not in ["PMC"]] + ["PM10", "PM25"]
    return list(set(variables))

############################
### Get what is the lowest frequency the station has
def get_freq_station(group, ID, types):
    #######  Files
    ## Path to database
    datapath=IncF.database[group]["path"]
    ## Files in directory that have the ID of station
    onlyfiles = [f for f in listdir(datapath) if isfile(join(datapath, f)) and "ID-%s--%s"%(ID, types) in f ]
    if len(onlyfiles)==0:
        return None
    ## Frequencies available in files
    freqs=[file.split("_")[-1].split(".")[0] for file in onlyfiles]
    vals=[aux.return_value_freq(freq)[1] for freq in freqs]


    #######   What was asked
    ## Numeric value of frequencies they want
    val_freq=[aux.return_value_freq(freq)[1] for freq in IncF.c_TimeMeanRes]
    ## Sorted list of frequencies they want
    freq_sort=[IncF.c_TimeMeanRes[val_freq.index(val)] for val in sorted(val_freq)]
    ## Minimum of frequency they want
    freq_min=freq_sort[0]


    #######
    ## Go through order of frequencies that are wanted.
    ## If same or lower frequency exists than minimum it will first return that
    ## return [frequency that exists, frequency wanted]
    for freq in freq_sort:
        ## If the frequency asked is available use this 
        if freq in freqs:
            return [freq, freq]
        ## Elif there is a lower frequency in files use that
        elif min(vals)<aux.return_value_freq(freq_min)[1]:
            return [freqs[vals.index(min(vals))], freq]

    ## If nothing was found then no frequency available for what was asked
    return None


############################
### SINCA  file name
def pname(ID, group, types, directories_db=None):
    ## Path to database
    datapath=directories_db["OBSERVATIONS"][group] 

    ## Files in directory that have the ID of station
    onlyfiles = [f for f in listdir(datapath) if isfile(join(datapath, f)) and "ID-%s--%s"%(ID, types) in f ]

    #print("\n", types, onlyfiles )
    ## Numeric value of frequencies they want
    val_freq=[aux.return_value_freq(freq)[1] for freq in IncF.c_TimeMeanRes]
    ## Sorted list of frequencies they want
    freq_sort=[IncF.c_TimeMeanRes[val_freq.index(val)] for val in sorted(val_freq)]
    ## Minimum of frequency they want
    freq_min=freq_sort[0]
    #print(freq_min)

    
    ## Frequencies available in files
    freqs=[file.split("_")[-1].split(".")[0] for file in onlyfiles]
    #print(freqs)

    if freq_min in freqs:
        file=onlyfiles[freqs.index(freq_min)]
        #print("1 File: ", join(datapath, file))
        return join(datapath, file)     
    if len(onlyfiles)==1:
        #print("1 File: ", join(datapath, onlyfiles[0]))
        return join(datapath, onlyfiles[0])
    elif len(onlyfiles)==0:
        #print("4 File: ", "No file")
        return "No File"
    ## If desired frequency is in files upload it
    ## NEed to add that
    else:
        ## Frequencies desired
        if freq_min in freqs:
            #print("2 File: ", join(datapath, "ID-%s--%s_%s.csv"%( ID, types, freq_min)))
            return  join(datapath, "ID-%s--%s_%s.csv"%( ID, types, freq_min))
        else:
            ## Minimum of the freqs available
            val_freq=[aux.return_value_freq(freq)[1] for freq in freqs]
            freq_min=freqs[val_freq.index(min(val_freq)) ]
            #print("3 File: ", join(datapath, "ID-%s--%s_%s.csv"%( ID, types, freq_min)))
            return join(datapath,"ID-%s--%s_%s.csv"%(ID, types, freq_min) )
            #common= list(set(freq).intersection(list2))
    
    
############################
## Find name of columns for variables
def find_equiv_variables(columns, variables):
    if "*" in variables:
        return  dict(zip(columns, columns))
     
    ## Transform Zonal variables into its parts
    if "WINDU" in variables or "WINDV" in variables:
        variables=[var for var in variables if var not in ["WINDU", "WINDV"]]+["WSPD", "WDIR"]
    if "D-U" in variables or "D-V" in variables:
        variables=[var for var in variables if var not in ["D-U", "D-V"]]+["WDIR"]
    if "PMC" in variables:
        variables=[var for var in variables if var not in ["PMC"]]+["PM10", "PM25"]

    variables=list(set(variables))
    cat=aux.read_yaml("vars")

    columns_met=[c for c in columns if c.split("--H")[0] in cat["Met"].split()] 
    columns_cal=[c for c in columns if c not in columns_met and (c not in ["Fecha", "date", "date_local", "date_utc"] and "CAT" not in c)]
    columns_cal=[c for c in columns_cal if "RadiaciA3n" not in c]
    ## Fix so that if value is not written to get error message
    columns_cal_right=[[var for var in cat["Equiv"] if c.split("_")[0] in cat["Equiv"][var]][0]  for c in columns_cal]
    #[var for var in cat["Equiv"] if 'MP2.5' in cat["Equiv"][var]][0]  
    #### Find the correspondoing name for the variables asked

    ## Variables whose name is without height specified
    dic_vars_cal={var:columns_cal[columns_cal_right.index(var)] for var in variables if var in columns_cal_right}
    
    ## Met variables equivalent  or those with height specified
    dic_vars_met={x:x.split("--H")[-1] for x in columns_met}
    ## Removes variables that dont have height
    dic_vars_met={k:(float(v)if v.isdigit() else 0) for (k,v) in dic_vars_met.items()}

    ## Removes variables whose heights are zero or are not within range
    if IncF.b_include_zero:
        dic_vars_met={k:v for (k,v) in dic_vars_met.items() if  v==0 or IncF.i_height-IncF.i_height_me<=v<=IncF.i_height+IncF.i_height_me}
    else:
        dic_vars_met={k:v for (k,v) in dic_vars_met.items() if IncF.i_height-IncF.i_height_me<=v<=IncF.i_height+IncF.i_height_me}
    

    dic_vars_met={var:list({k:v for (k,v) in dic_vars_met.items() if var in k}.values()) for var in variables if var not in dic_vars_cal}
    ## Not all wanted metereological variables have information in the wanted heights
    dic_vars_met={k:k+"--H"+str(int(min(v, key=lambda x:abs(x-IncF.i_height)))) for (k,v) in dic_vars_met.items() if v}

    dic_vars_met={k:(v if v in columns_met else k) for k,v in dic_vars_met.items()}

    ## All Variables equivalent
    dic_vars={**dic_vars_met, **dic_vars_cal}

    return dic_vars


############################
## Apply all filters from IncludeFile
def applyfilters(df_group, group):
    ## List of integers     
    list_filters = IncF.i_StationFilters
    ## Find index of group
    index = [count for count, g in enumerate(IncF.c_ObsNetwork) if g=="_".join(group.split("_")[:-1])]
    index = index[int(group.split("_")[-1])-1]

    t_filters={group:{}}
    df_group=gpd.GeoDataFrame(df_group, geometry=gpd.points_from_xy(df_group.Longitud, df_group.Latitud))

    #########
    ## No filter: read all
    if len(list_filters)==0:
        t_filters[group]["f0"]=[list(df_group["ID-Stored"].values)]
        return t_filters, df_group

    df_filters=pd.DataFrame()
    
    #########
    ## Apply corresponding filters
    if 1 in list_filters:
        t_filters[group]["f1"]={}
        ## Mask filter
        masks=IncF.c_RegionMask # _Mod[index]
        alias = IncF.c_RegionMaskAlias 
        print("   ...Filter  1:   Region Mask ", masks)
        for al, area in zip(alias, masks):
            if isinstance(area, str):
                ## Read geojson
                gdf = gpd.read_file(area, crs="epsg:4326")
                df_filters=df_filters.append({"ID":f"f1-{al}", "Latitud":gdf.centroid.y.values[0], "Longitud":gdf.centroid.x.values[0]}, ignore_index=True)
                new_mesh=gpd.overlay(df_group, gdf, how='intersection')
                if not new_mesh.empty:
                    t_filters[group]["f1"][al]=list(new_mesh["ID-Stored"].values)


    #########
    if 2 in list_filters:
        t_filters[group]["f2"]={}
        ## Geographical area
        squares=IncF.i_SqrRegions #_Obs[index]
        print("   ...Filter  2:   Square geographical area ", squares)
        alias=IncF.c_SqrRegionsAlias #_Obs[index]
        l_squares=[]
        for al, square in zip(alias, squares):
            df_filters=df_filters.append({"ID":f"f2-{al}", "Latitud":(square[1]+ square[0])/2, "Longitud":(square[2]+ square[3])/2}, ignore_index=True)
            df_group_f=df_group.loc[df_group['Latitud'].between(square[1], square[0])]
            df_group_f=df_group_f.loc[df_group_f['Longitud'].between(square[2], square[3])]
            if not df_group.empty:
                t_filters[group]["f2"][al]=list(df_group["ID-Stored"].values) 

    #########
    if 3 in list_filters:
        t_filters[group]["f3"]={}
        stations=IncF.c_Stations_Obs[index]
        alias=IncF.c_Stations_Obs_Aliases[index] 
        if len(stations)!=0:
            ## Stations list
            print("   ...Filter  3:   List of stations ", stations)
            l_stations=[]
            for al, station in zip(alias, stations):
                if len(station)!=0:
                    t_filters[group]["f3"][al]=[str(s) for s in station if str(s) in df_group["ID-Stored"].values]

                    cent=df_group[df_group["ID-Stored"].isin(station)].dissolve().centroid
                    df_filters=df_filters.append({"ID":f"f3-{al}", "Latitud":cent.y.values[0], "Longitud":cent.x.values[0]}, ignore_index=True)
    #########
    if 4 in list_filters:
        t_filters[group]["f4"]={}
        regions=IncF.c_Regions #_Obs[index]
        alias=IncF.c_RegionsAlias
        if len(regions)!=0:
            ## Chilean region
            print("   ...Filter  4:   By Chilean region ", regions)
            l_regions=[]
            for al, region in zip(alias, regions):
                region=[int(c) for c in region]
                stations=list(df_group.loc[df_group['Region'].isin(region)]["ID-Stored"].values)
                if len(stations)!=0:
                    cent=df_group[df_group["ID-Stored"].isin(stations)].dissolve().centroid
                    df_filters=df_filters.append({"ID":f"f4-{al}", "Latitud":cent.y.values[0], "Longitud":cent.x.values[0]}, ignore_index=True)
                    t_filters[group]["f4"][al]=stations
    #########
    if 5 in list_filters:
        t_filters[group]["f5"]={}
        comunas=IncF.c_Comunas #_Obs[index]
        if len(comunas)!=0:
            ## Chilean  Comuna
            print("   ...Filter  5:   By Chilean comuna ", comunas)
            for index, comuna in enumerate(comunas):
                stations=list(df_group.loc[df_group['cod_comuna'].isin(comuna)]["ID-Stored"].values)
                if len(stations)!=0:
                    cent=df_group[df_group["ID-Stored"].isin(stations)].dissolve().centroid
                    df_filters=df_filters.append({"ID":f"f5-{al}", "Latitud":cent.y.values[0], "Longitud":cent.x.values[0]}, ignore_index=True)
                    t_filters[group]["f5"][al]=stations

    #########
    if 6 in list_filters:
        t_filters[group]["f6"]={}
        elevation=IncF.i_Elevation #_Obs[index]
        if elevation!=0:
            ## Elevation
            print("   ...Filter  6:   By elevation ", elevation)
            l_elev=[]
            for el in elevation:
                al=str(el)
                df_mini=df_group.loc[pd.to_numeric(df_group['Elevacion'], errors='coerce').between(min(el), max(el))]
                cent=df_mini.dissolve().centroid
                df_filters=df_filters.append({"ID":f"f6-{al}", "Latitud":cent.y.values[0], "Longitud":cent.x.values[0]}, ignore_index=True)
                t_filters[group]["f6"][al]=list(df_mini["ID"].values)

    
    t_filters[group]={k:v for k, v in t_filters[group].items() if len(v)!=0}
    #########
    vals=[*t_filters[group].values()] 
    vals=[list(v.values()) for v in vals]
    vals=sum(sum(vals, []) , [])
    list_stations=list(set(vals))
    df_group=df_group[df_group["ID-Stored"].isin(list_stations)]
    return t_filters, df_group, df_filters


############################
## Find stations that have all variables asked
def filter_variables(df_group, group, variables, directories_db=None):
    return df_group
    ###########################################
    ############ Check variables outside of Summary CSV
    ## List of stations' ID
    IDs=[]
    ## List with number of variables each station has out of the ones in unused
    varia=[]    
    for ID in df_group["ID-Stored"]:
        ###########
        ## Read column names
        columns=[]
        try:
            Cal=pd.read_csv(pname(ID, group, "Cal", directories_db=directories_db), nrows=1)
            Cal.columns=[c.split("_")[0] for c in Cal.columns]
            columns=columns+list(Cal.columns)
        except Exception as e:
            pass
        try:
            Met=pd.read_csv(pname(ID, group, "Met", directories_db=directories_db), nrows=1)
            Met.columns=[c.split("_")[0] for c in Met.columns]
            columns=columns+list(Met.columns)
        except Exception as e:
            pass
        

        ###########
        #### Find the correspondoing name for the variables asked
        dic_vars=find_equiv_variables(columns, variables)
        

        ###########
        ## Checking that everything exists
        vars_have=list(dic_vars.values())
        if not len([var for var in vars_have if var in columns])==len(vars_have):
            ## There are missing variables so do not append station
            continue


        ## Append ID variable to list
        IDs.append(ID) 


    ##########
    ## Return summary information of all the stations that have all of the asked parameters
    df_group=df_group.loc[df_group["ID-Stored"].isin(IDs)]


    ## This dataframe has Lat, Lon, region, and all of the information needed to get the data into the MOSPAT layout.
    df_group.drop(list(df_group.filter(regex="Bool")), axis=1, inplace=True)
    return df_group


def filter_area_model(df_group, t_Models_Info):
    ## If models were read and no filter is applied for stations then consider only the stations that are within the area of the model
    if len(t_Models_Info)!=0 and len(IncF.i_StationFilters)==0:
        polygons = []
        dfs=[]
        for model in t_Models_Info:
            df=df_group.copy()
            Polygon=t_Models_Info[model].unary_union.convex_hull
            geometry=gpd.GeoSeries([Polygon])
            minx, miny, maxx, maxy = geometry.total_bounds
            df=df[df['Latitud'].between(miny, maxy)]
            df=df[df['Longitud'].between(minx, maxx)]
            dfs.append(df)
        new_df_group=pd.concat(dfs).drop_duplicates(keep="first")
        return new_df_group

    else:
        return df_group
######################################################
### Applies filters:
###       Only take stations that have the variables of interest
###  TODO:  Apply geographical filter. Look at mospat. Should work as is, the geographical one is currently being done later, so it is still fine.
def filter_stations(t_Models_Info, ggroup, variables, directories_db=None):
    ###########################################
    ## Get network without the number
    group="_".join(ggroup.split("_")[:-1])

    path_database=directories_db["OBSERVATIONS"][group] 
    df=pd.read_csv(path_database+"Chile-Stations.csv", dtype={"ID-Stored": str})
    


    ###########################################
    ## Filter by the general filters
    df_group = df.copy()
    df_group = filter_area_model(df_group, t_Models_Info)
    t_filters, df_group, df_filters = applyfilters(df_group, ggroup)


    ###########################################
    ## Keep only the stations that have all variables asked
    if "*" not in variables:
        df_group = filter_variables(df_group, group, variables, directories_db=directories_db)

    ######
    df_group = df_group[["ID-Stored", "Nombre", "Latitud", "Longitud", "Elevacion", "Region", "Comuna"]]
    df_group.columns = ["ID", "Nombre", "Latitud", "Longitud", "Elevacion", "Region", "Comuna"]
    df_group = df_group.append(df_filters, ignore_index=True)

    ###########
    ## Removes stations that dont have the variables from the filters
    list_stations=df_group["ID"].values
    for f in t_filters[ggroup].keys():
        ini_stations=t_filters[ggroup][f]
        #t_filters[ggroup][f]=[[s for s in group_f if s in list_stations] for group_f in ini_stations]
        t_filters[ggroup][f]={k:[s for s in v if s in list_stations] for k, v  in t_filters[ggroup][f].items() }


    
    return t_filters, df_group.reset_index(drop=True)

