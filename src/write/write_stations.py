# -*- coding: utf-8 -*-
from os import makedirs
from os.path import join
import src.IncludeFile as IncF
import pandas as pd

################################################
def save_station_like(data, info, freq, path=IncF.c_FigureDir):
    makedirs(path, exist_ok=True)
    ## Extract individual stations data
    stations=set(data.columns.get_level_values(1))
    for ID in stations:
        sta_data = data.xs(ID, axis=1, level=1).dropna(how="all")
        ## Save
        name=join(path, f"ID-{ID}_{freq}.csv")
        sta_data.to_csv(name)
        print("   ...Writing data to ", name)

    ## Save info stations
    info.to_csv(join(path, "Stations.csv"), index=False)



################################################
def main(t_ObsStationData, t_StationsInfo, t_Stations_Filters, model=None):
    ## Observations
    if model is None:
        for freq in t_ObsStationData:
            networks=list(set(["_".join(net.split("_")[:-1]) for net in list(t_ObsStationData[freq])]))
            for network in networks:
                ## Gather all pollutants that were read from one network
                data_network=[]
                for net in t_ObsStationData[freq]:
                    if network in net:
                        data_network.append(t_ObsStationData[freq][net])
                data_network=pd.concat(data_network, axis=1)

                ####
                stations=data_network.columns.get_level_values(1).values
                if 1 in IncF.i_Write_net:                         
                    selected_stations=[c for c in stations if (c[0]!="f") and (c[2]!="-")]
                    data_net=data_network.iloc[:, data_network.columns.get_level_values(1).isin(selected_stations)]
                    info=t_StationsInfo.loc[net]
                    info=info[info["ID"].isin(selected_stations)]
                    ## Save the concatenation of all pollutants
                    save_station_like(data_network, t_StationsInfo.loc[net], freq,  path=join(IncF.c_FigureDir, "Data", "Observations", network))
                if 2 in IncF.i_Write_net:
                    selected_stations=[c for c in stations if (c[0]=="f") and (c[2]=="-")]
                    data_net=data_network.iloc[:, data_network.columns.get_level_values(1).isin(selected_stations)]
                    info=t_StationsInfo.loc[net]
                    info=info[info["ID"].isin(selected_stations)]
                    ## Save the concatenation of all pollutants
                    save_station_like(data_net, info, freq,  path=join(IncF.c_FigureDir, "Data", "Observations", network+"_filters"))

    ## Models
    else:
        for freq in t_ObsStationData:
            for model in t_ObsStationData[freq]:
                networks=list(set(["_".join(net.split("_")[:-1]) for net in list(t_ObsStationData[freq][model])]))
                for network in networks:
                    c=0
                    for net in t_ObsStationData[freq][model]:
                        if c!=0:
                            continue
                        if network in net:
                            save_station_like(t_ObsStationData[freq][model][net], t_StationsInfo[model].loc[net], freq,  path=join(IncF.c_FigureDir, "Data", model, network) )
                            c=1
    