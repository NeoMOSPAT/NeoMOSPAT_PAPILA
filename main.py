# -*- coding: utf-8 -*-
import src.includefile_checker as checker
checker.run()
#
import src.common as aux
# Adds customs extensions to Pandas
import src.mospat_pandas_extensions
#
## read
from src.read import read_observations as read_obs
from src.read import read_models as read_model
#from src.read import read_aod
#
## manip
from src.manip import modelstationdata
from src.manip import bias_correction_SCM_Factor as bias_SCM_Factor
from src.manip import re_mapping as remap
#
## plots
#from src.plot import PCA
#from src.plot import Characteristics as char
from src.plot import Initial_Map as map_filters
#from src.plot import cross_section_plotter
from src.plot import Scatter as sct
from src.plot import Time_Series as ts
from src.plot import TS_ens as ts_ens
from src.plot import Time_Series_Summary as ts_sum
from src.plot import Time_Series_Statistics as ts_stat
from src.plot import TS_mosaic as ts_mosaic
from src.plot import Map_mosaic as maps_mosaic
from src.plot import Map_mosaic_models as maps_mosaic_mod
from src.plot import Surface_Maps as maps
from src.plot import Surface_Maps_Summary as maps_sum
from src.plot import Surface_Maps_Statistics as maps_stats
#
## write
from src.write import statistics as stats
from src.write import write_data as write_data
#
## Debug
from os import environ
environ["PYTHONBREAKPOINT"] = "pudb.set_trace"

import traceback
import warnings
import sys

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

warnings.showwarning = warn_with_traceback


"""
import sys
if not sys.warnoptions:
    import os, warnings
    warnings.simplefilter("error")
#"""



###################################################
###################################################
def main():
    print("................................................")
    print("...................Neo-MOSPAT...................")



    #############################################
    ### Model Data
    print("")
    print("")
    print(".................Reading Models..................")
    t_Model_Filters, t_Models_Info, t_ModelData, t_Cross_Sections = read_model.main()


    #############################################
    ### Observation Data
    print("")
    print("")
    print("..............Reading Observations..............")
    t_Stations_Filters, t_Stations_Info, t_ObsStationData = read_obs.main(t_Models_Info)



    #############################################
    ### Make summary map of run ares of interest
    print("")
    print("")
    print("..................Filters Map..................")
    map_filters.main(t_Stations_Info, t_Stations_Filters)



    #############################################
    ### Model Station Data
    print("")
    print("")
    print("...............Model-Station Data...............")
    t_ModelStation_Info, t_ModelStationData = modelstationdata.main(
        t_Stations_Info,
        t_ObsStationData,
        t_Models_Info,
        t_ModelData, 
        t_Model_Filters
    )


    print("")
    print("")
    print("...................Write Data...................")
    write_data.main(t_ObsStationData, t_ModelStationData, t_ModelData, t_Stations_Info, t_ModelStation_Info, t_Models_Info, t_Stations_Filters, t_Model_Filters)


    #############################################
    ### Plot general behaviour of data
    #print(" ")
    #print(" ")
    #print(".........Plot:  General Characteristics.........")
    #char.main(t_Stations_Info, t_ObsStationData)



    #print("")
    #print("")
    #print("..................Cross Section.................")
    #cross_section_plotter.main(t_Cross_Sections=t_Cross_Sections, t_ModelFilters=t_Model_Filters)



    print("")
    print("")
    print("...................Surface Maps.................") 
    maps.mainn(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info)
    maps_mosaic_mod.main(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info)
    maps_mosaic.main(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info)
    maps_sum.main(t_ObsStationData, t_ModelData, t_Stations_Info, t_Model_Filters, t_Models_Info)


    #maps_stats.main(t_ObsStationData, t_Stations_Info, t_ModelStationData,  t_ModelStation_Info,  t_ModelData, t_Model_Filters, t_Models_Info)



    print("")
    print("")
    print("...................Time Series..................")
    ts.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData)
    ts_sum.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData)
    ts_stat.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData)
    #ts_ens.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData)
    #ts_mosaic.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData)
    #ts_ens.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters, t_ModelData)   




    print("")
    print("")
    print(".....................Scatter....................")   
    sct.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_ModelStation_Info, t_Model_Filters, t_Stations_Filters)



    print("")
    print("")
    print("...................Statistics...................")
    stats.main(t_ObsStationData, t_ModelStationData, t_Stations_Info, t_Model_Filters, t_Stations_Filters)




    ############################################################
    #### Things that are so niche that we dont want them popping up all the time
    #print("...................Re Mapping...................")
    remap.main(t_ObsStationData,  t_ModelData, t_Stations_Info, t_Models_Info)
    
    #print(".................Bias Correction................")
    t_BiasCorrModelData,  t_BiasCorrModels_Info = bias_SCM_Factor.main(t_ObsStationData, t_ModelStationData, t_ModelData, t_Stations_Info, t_ModelStation_Info, t_Models_Info)



############################################################
if __name__== "__main__":
    main()
