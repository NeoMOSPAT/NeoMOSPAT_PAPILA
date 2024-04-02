# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from functools import reduce

import src.IncludeFile as IncF
#
import src.common as aux
#
from os.path import join
from os import makedirs
import itertools

IncF.i_Stats=[int(c) for c in str(IncF.i_Stats)]
#####################################
def get_average_stations(data_or, group):
    ## Retrieve data of stations of group
    data=aux.get_stations(data_or, group)
    ## Create average for each variable
    data=data.mean(axis=1, level=0)
    return data


#####################################
## Make Statistics
def make_stats(data_models, data_obs, column_obs, freq="D"):
    data_models=data_models.resample(freq).mean()
    data_obs=data_obs.resample(freq).mean()

    ## Substract models with obs
    models_subs = data_models.sub(data_obs[column_obs], axis=0)
    ## Add models with obs
    models_add = data_models.add(data_obs[column_obs], axis=0)
    ## Obs average
    obs_avg = data_obs.mean()
    ## Create model Average
    models_avg = data_models.mean()
    div_avg = models_avg/obs_avg.values[0]
    ## Sigmas
    models_sigma = np.sqrt( ( (data_models-models_avg)**2 ).mean() )
    obs_sigma = np.sqrt( ( (data_obs-obs_avg.values[0])**2 ).mean() )
    div_sigma = models_sigma/obs_sigma.values[0]

    ################
    data_count = 100*data_models.count()/(data_models.shape[0]-1)
    obs_count = 100*data_obs.count()/data_obs.shape[0]

    ##
    data_bias = models_subs.mean()
    data_rmse = np.sqrt( (models_subs**2).mean())
    data_nmbias = (2*models_subs/models_add).mean()
    data_fge = abs(2*models_subs/models_add).mean()
    data_r = ((data_models-models_avg).mul((data_obs-obs_avg.values[0])[column_obs], axis=0)).mean()/(models_sigma.mul(obs_sigma[column_obs], axis=0))
    data_coefvar_obs = obs_sigma/obs_avg
    data_coefvar_mod = models_sigma/models_avg

    df=pd.concat([models_avg, models_sigma, data_bias, data_nmbias, data_rmse, data_fge, data_r, div_sigma, div_avg, data_coefvar_mod, data_count ], keys=["Mean", "STD", "Bias", "NMBias", "RMSE", "FGE", "R", "razon_STD", "razon_AVG", "CV", "Count"], axis=1)
    df.index=df.index
    df.loc["Obs", "Mean"]=obs_avg.values[0]
    df.loc["Obs", "STD"]=obs_sigma.values[0]
    df.loc["Obs", "CV"]=data_coefvar_obs.values[0]
    df.loc["Obs", "Count"]=obs_count.values[0]
    print(df)
    #breakpoint()
    return df.T



def make_statistics(data, variables, freq="D", direc=""):
    merge={}
    for var in variables:
        data_var=data[[c for c in data.columns if var in c]]
        data_var.columns=[c.split(" ")[1] for c in data_var.columns]

        ## If there is network and model variable
        if data_var.shape[1]!=1:
            data_obs=data_var[[data_var.columns[0]]]
            data_models=data_var[[c for c in data_var.columns if c not in data_obs.columns]]
            if data_obs.dropna().empty:
                print("         ...Observations dont have the variable  ", var)
                continue
            df = make_stats(data_models, data_obs, data_var.columns[0], freq=freq)          
            merge[var]=df
        else:
            print("         ...Models dont have the variable  ", var)

    if bool(merge):
        return pd.concat(merge, axis=1, keys=variables).dropna(how="all")
    else:
        return pd.DataFrame()


#####################################
## Make previous step before making statistics
def make_csv(ID, t_ObsData, t_ModData, freq="D", direc=""):
    network=list(t_ObsData.keys())[0] 
    clean_network="_".join(network.split("_")[:-1]) 
    net_data=t_ObsData[network]
    variables=list(t_ObsData[network].keys())
    variables=[v for v in variables if v!="groupby"]
    net_data.index.name="Time"
    net_data.columns=["%s  %s %s"%(var, clean_network, ID) for var in net_data.columns]

    data=[net_data]

    for model in t_ModData:
        mod_data=t_ModData[model]
        #variables=variables+list(t_ModData[model].keys())
        mod_data.index.name="Time"
        mod_data.columns=["%s %s  %s %s"%(var, model, clean_network, ID) for var in mod_data.columns]
        data.append(mod_data)

    #print(net_data)
    #print(mod_data)
    data=reduce(lambda df1, df2: pd.merge(df1, df2, left_index=True, right_index=True, how="outer"), data)

    dic = make_statistics(data, list(set(variables)), direc=direc, freq=freq)
    """
    ##################
    dic={}
    for timeframe in zip(IncF.c_Start_Date, IncF.c_Last_Date):
        date1="".join(timeframe[0].split("-")[::-1])
        date2="".join(timeframe[1].split("-")[::-1])
        #print(data)
        dic[f"{date1}-{date2}"]=make_statistics(data.loc[date1:date2], list(set(variables)), direc=direc, freq=freq)
    """
    return dic


#####################################
def main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_Model_Filters, t_Stations_Filters):
    if 0 in IncF.i_Stats:
        print("...Statistics were not requested.")
        return None
    if t_Stations_Info.empty:
        print("...No stations were read. Skipping this step.")
        return None 
        
    ## Check models and obs that ended up being read
    models=list(t_ModelStationData[list(t_ModelStationData.keys())[0]].keys())
    obs=list(t_ObsStationData[list(t_ObsStationData.keys())[0]].keys())   

    path=join(IncF.c_FigureDir, "Statistics")
    makedirs(path, exist_ok = True) 
    for plot_freq in IncF.i_Stats_freq:
        freq, plot_freq = plot_freq.split(",")
        for network in obs:
            print("...Creating stats for network:  ", "_".join(network.split("_")[:-1]))
            info=aux.get_info_group(t_Stations_Info, network)
            clean_network="_".join(network.split("_")[:-1]) 
             

            data_observations = t_ObsStationData[freq][network]
            data_model = t_ModelStationData[freq]

            if plot_freq=="Whole":
                c=0
                for timeframe in zip(IncF.c_Start_Date, IncF.c_Last_Date):
                    c+=1
                    date1="".join(timeframe[0].split("-")[::-1])
                    date2="".join(timeframe[1].split("-")[::-1])
                    data_observations.loc[date1:date2, "groupby"] = c 
                grouped_obs=data_observations.groupby("groupby")
            else:
                grouped_obs=data_osservations.groupby(pd.Grouper(freq=plot_freq))

            list_dates=list(grouped_obs.groups.keys())
            data={}
            for group in list_dates:
                try:
                    data_obs=grouped_obs.get_group(group)
                except:
                    data_obs=pd.DataFrame()

                if data_obs.empty:
                    continue  

                ## timeframe
                date1 = data_obs.index[0].strftime("%Y-%m-%d")
                date2 = data_obs.index[-1].strftime("%Y-%m-%d") 
                timeframe =  date1 + "-" + date2 
                data[timeframe] = {}


                ##############
                ## Individual plots  
                print("   ...Creating stats for timeframe  ", timeframe)
                for station in  info["ID"]:
                    available_sta_mod={}
                    vailable_sta_obs=pd.Series()
                    if station[0]=="f" and station[2]=="-":
                        stations_group = t_Stations_Filters[network][station.split("-")[0]][station.split("-")[1]]
                        available_sta_mod = {mod:data_model[mod][network].loc[date1:date2, (slice(None), stations_group)].dropna(how="all", axis=1).groupby(level=0, axis=1).size() for mod in models } 
                        

                        available_sta_obs = data_obs.loc[:, (slice(None), stations_group)].dropna(how="all", axis=1).groupby(level=0, axis=1).size()
                        available_sta_obs = pd.DataFrame(available_sta_obs)
                        available_sta_obs["level"] = "Obs"
                        available_sta_obs = available_sta_obs.set_index("level", append=True)
                        available_sta_obs.index.names = [None, None]


                        if 2 not in IncF.i_Stats:
                            continue
                    else:
                        if 1 not in IncF.i_Stats:
                            continue

                    print("      ...Creating stats for ", station)
                    ## Retrieve data of stations of group
                    data_net={network:aux.get_station_data(data_obs, station)}
                    ##
                    data_mod={mod:aux.get_station_data(data_model[mod][network], station).loc[date1:date2] for mod in models}

                    ############
                    df=make_csv(station, data_net, data_mod, direc=path, freq=freq)
                    
                    if bool(available_sta_mod):
                        available_sta_mod = pd.concat(available_sta_mod).swaplevel()
                        merge = pd.concat([available_sta_obs, pd.DataFrame(available_sta_mod)]).T
                        merge.index = ["# Stations"]
                        df = pd.concat([df, merge])

                    data[timeframe][station]=df


                mini=info.copy()
                mini=mini[~mini["ID"].str.contains('^f[0-9]-+')]
                minif=info[info["ID"].str.contains('^f[0-9]-+')]
                if 1 in IncF.i_Stats:
                    print("   ...Creating stats for all stations")
                    ## Retrieve data of stations of group
                    data_net=[aux.get_station_data(data_obs, station) for station in mini["ID"]]
                    data_net={clean_network:pd.concat(data_net)}
                    ##
                    data_mod={mod:[aux.get_station_data(data_model[mod][network], station) for station in mini["ID"]] for mod in models }   
                    data_mod={mod:pd.concat(data_mod[mod]) for mod in models}   
                    data[timeframe]["All-Stations"]=make_csv(clean_network, data_net, data_mod, direc=path, freq=freq)

            
            def function(x):
                ## Obs have no info
                if x.shape[0]==1:
                    return None
                nstations = x.loc[[c for c in x.index if c[1]=="# Stations"]][[c for c in x.columns if c[1]=="Obs"]].values[0][0]
                meanvalue = x.loc[[c for c in x.index if c[1]=="Mean"]][[c for c in x.columns if c[1]=="Obs"]].values[0][0]
                availabil = x.loc[[c for c in x.index if c[1]=="Count"]][[c for c in x.columns if c[1]=="Obs"]].values[0][0]
                
                ## Obs have no info but models do
                if np.isnan(meanvalue):
                    return None
                ## No stations
                elif np.isnan(nstations) or nstations==0:
                    return 
                ## No data
                elif np.isnan(availabil) or availabil==0:
                    return None
                else:
                    x = x.T
                    ## For individual stations
                    x = x.fillna({c:1 for c in x.columns if c[1]=="# Stations"}) 
                    x = x.T
                    for col in x.columns:
                        for city in x.index.get_level_values(0):
                            avail = x.loc[(city,"Count")][col]
                            if avail==0:
                                x.loc[(city, slice(None)), col] = np.nan
                    return x.reset_index(drop=True, level=0)

            #####################################
            for timeframe, df in data.items():
                ## Add Stations=1 for statistics of individual stations
                df1 = pd.concat(df)
                df1.index.names = ["Station", "Stats"]
                #df1 = df1.reset_index()
                df1 = df1.groupby("Station").apply(lambda x: function(x))
                #df1 = df1.fillna({c:1 for c in df1.columns if c[1]=="# Stations"}) 
                #df = df1.T

                for var in set(df1.columns.get_level_values(0)):
                    makedirs(f"{path}/{freq},{plot_freq}/{clean_network}", exist_ok = True)
                    df_var=df1[var]
                    #print(df_var)
                    #print(df_var.shape)
                    name=f"{path}/{freq},{plot_freq}/{clean_network}/{var}__{clean_network}__{timeframe}.csv"
                    df_var.to_csv(name, float_format='%.4f')
                    print("   ...Saving to:  ", name)


