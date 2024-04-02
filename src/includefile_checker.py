# -*- coding: utf-8 -*-
# In this file are checked all restrictions for IncludeFile options. This helps 
# to eliminate double checking along the code, simplifying it and to catch
# errors earlier in the execution. Also, makes some transformation to data when
# necessary.

import sys
from datetime import datetime



######################################################
old_out = sys.stdout
class St_ampe_dOut:
    """Stamped stdout."""
    nl = True
    def flush(self):
        pass
    def write(self, x):
        """Write function overloaded."""
        if x == '\n':
            old_out.write(x)
            self.nl = True
        elif self.nl:
            old_out.write('[%s]    %s' % (str(datetime.now().strftime("%d/%m/%Y %H:%M:%S")), x))
            self.nl = False
        else:
            old_out.write(x)

sys.stdout = St_ampe_dOut()



__all__ = ['run']
## Converts nested lists of objects to nested lists of strings
def list_to_string(x):
    if isinstance(x, list):
        return list(map(list_to_string, x))
    else:
        return str(x)

## Converts nested lists of objects to nested lists of strings
def list_to_int(x):
    if isinstance(x, list):
        return list(map(list_to_int, x))
    else:
        return int(x)



def run(ignore_exceptions=False):
    """
    Perform sanity checks of IncludeFile params

    Adds the following variables to IncludeFile:

    - d_Start_Date and d_Last_Date: dates ranges as lists of datetime.datetime
    objects.
    - c_Filters: list of characters made from i_Filters

    Args:
        ignore_exceptions: bool default False. If True all exceptions raised
        by no 'assert' statements are ignored, but are still reported through
        standard output. Only for debugging purposes.
    """
    from datetime import datetime
    from numbers import Number
    import numpy as np
    import pandas as pd
    import yaml
    import os
    from  shutil import copy2
    
    #IF=__import__(IncludeFile.replace(".py", ""))


    """
    TODO: When doing logging add color
    def d(*v): 
        return '\x1B['+';'.join(map(str, v))+'m' 

    print(' '.join([d(k,i)+str(i%10)+d(0) for i in list(range(30,38))+list(range(40,48)) for k in range(2)]))
    print(' '.join([d(k,i)+"k:%s  i:%s"%(k, i)+d(0) for i in list(range(30,38))+list(range(40,48)) for k in range(2)]))
    CSI = "\x1B["
    print(CSI+"31;40m" + "Colored Text" + CSI + "0m")
    print(CSI+"31;40m" + u"\u2588" + CSI + "0m")
    #https://stackoverflow.com/a/287944/5513902
    """

    print("................................................")
    print("..............Checking IncludeFile..............")

    if ignore_exceptions:
        print(f"...WARNING!! YOU'RE IGNORING EXCEPTIONS AT {run}")

    ##############################################
    #####             IncludeFile            #####
    ##############################################  
    ### Copy Includefile to target so that more than one run can be made
    try:
        IncludeFile = sys.argv[1]
    except:
        IncludeFile = "IncludeFile.py"

    if os.path.isfile(IncludeFile):
        copy2(IncludeFile, "src/IncludeFile.py")
        os.system("cat config/IncludeFileBase.py %s > src/IncludeFile.py"%IncludeFile)
    else:
        print(f"...Used IncludeFile  '{IncludeFile}'  doesn't exist.")
        print("...Check the path of your file and how you ran NeoMOSPAT.")
        print("   ...Correct ways running the script are:\n                                  '$ python main.py YourIncludeFile.py'  \n                                  '$ python main.py'      (searches for IncludeFile.py by default) ")
        sys.exit()
        assert False
        #, f"...IncludeFile ile  {IncludeFile}  doesn't exist. Check your file.\nCorrect usage:    python main.py YourIncludeFile.py  or  python main.py (searches for `IncludeFile.py` by default) " 



    ##############################################
    #####               GENERAL              #####
    ##############################################
    import src.IncludeFile as IF
    print("...Checking general settings")

    ############
    ## Check the variable type of the input in IncludeFile
    assert isinstance(IF.c_TimeMeanRes, list), "   ...c_TimeMeansRes must be a list"
    assert len(IF.c_TimeMeanRes)!=0, "   ...c_TimeMeansRes has no specified frequency of data to be read"
    freqs_avail = ["B", "C", "D", "3D", "W", "M", "SM", "BM", "CBM", "MS", "SMS", "BMS", "CBMS", "Q", "BQ", "QS", "BQS", "A", "Y", "BA", "BY", "AS", "YS", "BAS", "BYS", "BH", "H", "T", "min", "S", "L", "ms", "U", "us", "N"]
    freq_asked_without_numbers=set([f.lstrip('0123456789') for f in IF.c_TimeMeanRes])
    assert freq_asked_without_numbers.intersection(freqs_avail) == freq_asked_without_numbers, "   ...Unrecognized values in c_TimeMeanRes"

    assert has_nested_type(IF.c_Start_Date, list, str), "   ...c_Start_Date must be a \
        list of strings dates written in dd-mm-yyyy format"
    assert has_nested_type(IF.c_Last_Date, list, str), "   ...c_Last_Date must be a \
        list of strings dates written in dd-mm-yyyy format"
    assert len(IF.c_Start_Date) == len(IF.c_Last_Date), "   ...c_Start_Date and \
        c_Last_Date must have the same length"

    ############
    ## date to Datetime objects
    IF.d_Start_Date = [datetime.strptime(date, "%d-%m-%Y") for date in IF.c_Start_Date]
    IF.d_Last_Date = [datetime.strptime(date, "%d-%m-%Y") for date in IF.c_Last_Date]


    ############
    ## Changes directory to save to include subfolder if requested to
    if IF.b_CreateNameSubFolder:
        IF.c_FigureDir=os.path.join(IF.c_FigureDir,IF.c_RunName.replace(" ", "_"))
    ## Creates figure directory
    if IF.c_FigureDir and not os.path.exists(IF.c_FigureDir):
        try:
            os.makedirs(IF.c_FigureDir, exist_ok = True)
        except PermissionError as e:
            print(f"   ...{e}. Can't create c_FigureDir")
            if not ignore_exceptions:
                raise e



    ##############################################
    #####               MODELS               #####
    ##############################################
    print("...Checking Models")
    ############
    ## Check the variable type of the input in IncludeFile
    assert has_nested_type(IF.c_ModelNames, list, str), "   ...c_ModelNames must be a \
        list of strings"

    assert has_nested_type(IF.c_ModelAlias, list, str), "   ...c_ModelAlias must be a \
        list of strings"

    assert isinstance(IF.i_TimeZone, int), "   ...i_TimeZone must be an integer"

    assert has_nested_type(IF.c_ModelVars, list, str), "   ...c_ModelVars must be a \
        list of strings"

    assert isinstance(IF.b_OnlyLand, bool), "   ...b_OnlyLand should be True or False"

    c_model_time_freqs = ('PP', 'DD')
    assert IF.c_ModelTimeFreq in c_model_time_freqs, \
           f"   ...Unrecognized c_ModelTimeFreq {IF.c_ModelTimeFreq}. Should be one of {c_model_time_freqs}"


    ############
    ## Check if model directory exists
    assert os.path.exists(IF.c_ModelDir), "   ...Model directory does not exist in \
        your filesystem"
    ## Check all models are different
    assert not any(IF.c_ModelNames.count(element) > 1 for element in IF.c_ModelNames),  "   ...There is a duplicate model in existence in  c_ModelNames. Each element of that list must be unique"
    ## Check for the existence of the model data
    folders_in_dir=[x[0].split("/")[-1] for x in os.walk(IF.c_ModelDir)]
    #assert len([v for v in IF.c_ModelNames if  v not in folders_in_dir]) == 0, f"   ...A model requested does not exist. Models available for dir  {IF.c_ModelDir}  are  {folders_in_dir}"
    for mod in IF.c_ModelNames:
        if mod not in folders_in_dir:
            assert len([v for v in IF.c_ModelNames if  v not in folders_in_dir]) == 0, f"   ...A model requested does not exist. Models available for dir  {IF.c_ModelDir}  are  {folders_in_dir}"



    ############
    ## Check if Alisases exist correctly
    if len(IF.c_ModelNames) != len(IF.c_ModelAlias):
        # TODO: this should use logging.warning
        # TODO: rellenar los model aliases
        print("   ...WARNING: Number of model aliases and names does not match => \
            aliases have been replaced by model names")
        IF.c_ModelAlias=IF.c_ModelNames
    ## Check if there are duplicate alias
    assert not any(IF.c_ModelAlias.count(element) > 1 for element in IF.c_ModelAlias), "   ...There is a duplicate alias in existence in  c_ModelAlias. Each element of that list must be unique"
    ## Check that the amoung of characters is below a certain threshold
    if not np.all([len(alias)==15 for alias in IF.c_ModelAlias]): 
        print("   ...WARNING: Aliases are longer than 15 characters. They have been truncated to the smallest possible that will keep them unique.")
    ## Try to truncate to 15 characters. If that creates repeated aliases then keep increasing one character until they are all different aliases.
    b_allsame=True
    c=-1
    while(b_allsame):
        c+=1
        new_model_alias=[alias[:15+c].replace(" ", "_") for alias in IF.c_ModelAlias] 
        b_allsame=any(new_model_alias.count(element) > 1 for element in new_model_alias)       
    IF.c_ModelAlias=new_model_alias

    


    ##############################################
    #####           OBSERVATIONSS            #####
    ##############################################
    print("...Checking Observations")

    ############
    ## Check the variable type of the input in IncludeFile
    assert has_nested_type(IF.c_ObsNetwork, list, str), "   ...c_ObsNetwork must be a list of strings"

    assert has_nested_type(IF.c_ObsVars, list, list, str), "   ...c_ObsVars must be a list of lists of strings"

    assert isinstance(IF.i_height, int), "   ...Height for MET data should be an integer"

    assert isinstance(IF.i_height_me, int), "   ...Height margin error should be an integer"

    assert isinstance(IF.b_include_zero, bool), "   ...b_include_zero should be True or False"

    assert isinstance(IF.i_date_column, int), "   ...i_date_column should be an integer"
    
    if IF.i_date_column==1:
        IF.i_date_column = "date_utc"
    elif IF.i_date_column==2:
        IF.i_date_column = "date_local"
    else:
        IF.i_date_column = "date"

    ############
    ## Check that the variable requested exists
    data_dirs=yaml.load(open("config/data_directories.yaml"), Loader=yaml.FullLoader)
    unused_vars_net={}
    for count, network in enumerate(IF.c_ObsNetwork):
        vars_asked=IF.c_ObsVars[count]
        vars_network=list(data_dirs["UNITS"].keys())
        unused_vars=[]
        for v in vars_asked:
            if v not in vars_network:
                unused_vars.append(v)

        if len(unused_vars)!=0:
            print(f"   ...Variables  {', '.join(unused_vars)}  for network  {network}  don't exist")
            unused_vars_net[network]=unused_vars

        #assert len([v for v in vars_asked if  v not in vars_network]) == 0, f"   ...At least one requested variable for network {network} doesnt exist. Allowed variables are: {vars_network}"

    if len(unused_vars_net)!=0:
        assert False, f"   ...At least one requested variable for network doesnt exist. Check these variables:  {unused_vars_net}"


    for count, network in enumerate(IF.c_ObsNetwork):
        vars_asked=IF.c_ObsVars[count]
        with open("config/equiv_vars.yaml") as f:
            cat = yaml.safe_load(f)

        columns_cal=[c for c in vars_asked if c in cat["Cal"].split()]
        columns_met=[c for c in vars_asked if c in cat["Met"].split()]

        if len(vars_asked)!=len(columns_cal)+len(columns_met):
            print("   ...Not all variables requested are categorized.") 
            print("   ...Please go to  config/equiv_vars.yaml and define all variables available in data base according to if they are metereological or air quality variables")
            assert False, "   ...Not all variables requested are categorized"


    ##############################################
    #####               FILTERS              #####
    ##############################################
    print("...Checking Filters")

    ############
    ## Model Filters
    if not IF.i_ModelFilters == 0:
        assert isinstance(IF.i_ModelFilters, int), "   ...i_ModelFilters must be an integer"
        c_Filters = list(str(IF.i_ModelFilters))
        assert len(set(c_Filters).difference(['1','2','3','4','5','6'])) == 0, "   ...You have set unrecognized filters codes. Only 1-6 are allowed"
        IF.c_ModelFilters = c_Filters
        IF.i_ModelFilters = list(list_to_int(c_Filters))
    else:
        IF.c_ModelFilters = ['0']
        IF.i_ModelFilters = []

    ############
    ## Station Filters
    if not IF.i_StationFilters == 0:
        assert isinstance(IF.i_StationFilters, int), "   ...i_ModelFilters must be an integer"
        c_Filters = list(str(IF.i_StationFilters))
        assert len(set(c_Filters).difference(['1','2','3','4','5','6'])) == 0, "   ...You have set unrecognized filters codes. Only 1-6 are allowed"
        IF.c_StationFilters = c_Filters
        IF.i_StationFilters = list(list_to_int(c_Filters))

        if 3 in IF.i_StationFilters:
            assert len(IF.c_ObsNetwork) <= len(IF.c_Stations_Obs), "You have more stations filters than networks"

            ## TODO: Document when implemented and remove comment
            #assert has_nested_type(IF.c_RegionMask, list, str), "c_RegionMask must \
            #    be a list f strings"

            assert has_nested_type(IF.c_Stations_Obs, list, list, list), "   ...c_stations must be a list of lists of lists."
            IF.c_Stations_Obs=list(list_to_string(IF.c_Stations_Obs))
    else:
        IF.c_StationFilters = ['0']
        IF.i_StationFilters = []


    ############
    ## Shared Filters
    if not IF.i_ModelFilters==[0] or IF.i_StationFilters==[0]:
        ########
        ## F1: Shapefile and or netcdf. Undefines yet
        #assert not (len(IF.c_ObsNetwork) < len(IF.i_RegFilterLimits)), "You have defined more region filter lists than observation networks"


        ########
        ## F2: Square regions
        assert has_nested_type(IF.i_SqrRegions, list, list, Number), "   ...You can define many square regions for each network, make sure they are correctly setup in i_SqrRegions"

        # Region limits checks.
        # Latitudes from -90 to 90 and longitudes from -180 to 180
        i_region_limits = np.array(IF.i_SqrRegions, dtype=np.int32).flatten()
        assert len(i_region_limits) % 4 == 0, "   ...Some of yours regions at i_RegFilterLimits lack of coordinates"

        i_regs = int(len(i_region_limits) / 4)
        lats_selector = [True, True, False, False] * i_regs
        lons_selector = [False, False, True, True] * i_regs

        assert all(np.abs(i_region_limits[lats_selector]) <= 90), "   ...Latitudes must be in -90, 90 range"

        assert all(np.abs(i_region_limits[lons_selector]) <= 180), "   ...Longitudes must be in -180, 180 range"
        
        assert has_nested_type(IF.c_SqrRegionsAlias, list, str), "   ...c_RegFilterAlias must be a list of lists of strings"


        ########
        ## F3: Chilean's Region
        assert has_nested_type(IF.c_Regions, list, list), "   ...c_Regions must be a list of lists with the regions of interest."
        IF.c_Regions=list(list_to_string(IF.c_Regions))


        ########
        ## F5: Chilean's Comunas
        assert has_nested_type(IF.c_Comunas, list, list), "   ...c_comunas must be a list of lists."
        IF.c_Comunas=list(list_to_string(IF.c_Comunas))
        comunas_codes=pd.read_csv("config/comunas_codes.csv")["cod_comuna"].astype(str).values
        assert not any([(c not in comunas_codes) for g in IF.c_Comunas for c in g]), "   ...Not all comunas in   c_Comunas  exist. Check the file 'config/comunas_codes.csv' to know which ones are available"



        ########
        ## F6: Elevation 
        assert has_nested_type(IF.i_Elevation, list, list, int), "   ...i_elevation must be a list of lists of integers"



    ##############################################
    #####           CROSS SECTIONS           #####
    ##############################################
    print("...Checking Plots")
    print("   ...Checking Cross-Sections")
    # in this case it only checks that dict keys are strings

    assert has_nested_type(IF.t_cross_sections_config, list, dict, str), \
    "      ...t_cross_sections_config must be a list of dicts with string keys"

    for cross_section_dict in IF.t_cross_sections_config:
        # check that all cross sections have 'path' field and at two points
        assert 'path' in cross_section_dict, "      ...Missing field 'path' in cross section config"
        assert has_nested_type(cross_section_dict['path'], list, tuple, Number)



    ##############################################
    #####               PLOTS                #####
    ##############################################
    ##############################
    ## Opens up representation of plots, opens up * to numbers
    def open_plots(plots_IF):
        try:
            plots=list(set(list(np.hstack(plots_IF).flatten())))
        except:
            return []
        ########
        ## Open first value
        new_plots=[]
        for var in set(plots):
            split_version=var.split("_")
            if split_version[1]=="*":
                if "O" in var:
                    len_obs=len(IF.c_ObsNetwork)
                    for i in range(len_obs):
                        new_plots.append("O_%s_%s"%(i+1, split_version[2]))
                elif "M" in var:
                    len_mod=len(IF.c_ModelNames)
                    for i in range(len_mod):
                        new_plots.append("M_%s_%s"%(i+1, split_version[2]))
            else:
                indexes = split_version[1].split(",")
                for index in indexes:
                    split_version[1]=index
                    new_plots.append("_".join(split_version))
        #######
        ## Open second value
        plots=[]
        for var in set(new_plots):
            split_version=var.split("_")
            if split_version[-1]=="*":
                if "O" in var:
                    len_vars=len(IF.c_ObsVars[int(split_version[1])-1])
                    for i in range(len_vars):
                        plots.append("O_%s_%s"%(split_version[1], i+1))
                elif "M" in var:
                    len_vars=len(IF.c_ModelVars)
                    for i in range(len_vars):
                        plots.append("M_%s_%s"%(split_version[1], i+1))
            else:
                indexes = split_version[2].split(",")
                for index in indexes:
                    split_version[2]=index
                    plots.append("_".join(split_version))

        return list(set(plots))

    ## Makes sure the variables defined within plots were asked to be read
    def check_plots_vars_exist(variables):
        for var in variables:
            split=var.split("_")
            split_1=int(split[1])
            split_2=int(split[2])

            if "O" in var:
                try:
                    IF.c_ObsNetwork[split_1-1]
                except:
                    raise AssertionError(f"      ...Network O_{split_1} doesn't exist. Check  IncludeFile  parameters" )
                try:
                    IF.c_ObsVars[split_1-1][split_2-1]
                except:
                    raise AssertionError(f"      ...Network var O_{split_1}_{split_2} doesn't exist. Check  IncludeFile  parameters" )            


            if "M" in var:
                try:
                    IF.c_ModelNames[split_1-1]
                except:
                    raise AssertionError(f"      ...Model M_{split_1} doesn't exist. Check  IncludeFile  parameters" )
                try:
                    IF.c_ModelVars[split_2-1]
                except:
                    raise AssertionError(f"      ...Model var M_{split_1}_{split_2} doesn't exist. Check  IncludeFile  parameters" )


    #####################
    print("   ...Checking Time Series")
    ## If Time series requested
    if IF.i_TS!=0:
        ## Check frequencies were asked to be read
        freqs=list(set([f.split(",")[0] for f in IF.c_TS_freq]))
        assert set(freqs).intersection(freqs_avail) == set(freqs), "   ...Unrecognized plot frequency in c_TS_freq. Check frequencies specified before the comma"
        for freq in freqs:
            if freq not in IF.c_TimeMeanRes:
                print(f"      ...WARNING: requested freq  {freq}  doesn't exist in  c_TimeMeanRes. It will be added to variable")
                IF.c_TimeMeanRes=IF.c_TimeMeanRes+[freq]

        ## Check frequencies of plots that were asked
        freqs=set([f.split(",")[1] for f in IF.c_TS_freq])
        assert freqs.intersection(freqs_avail+["Whole"]) == freqs, "   ...Unrecognized plot frequency in c_TS_freq. Check frequencies specified after the comma"

        ## Check that plots requested exist
        plots=open_plots(IF.c_TS_Plots)
        #check_plots_vars_exist(plots)

        ## Check requested plot types exist
        if len(set([s for s in str(IF.i_TS)]).difference(['0', '1'])) != 0:
            print("      ...WARNING: You have set an unrecognized value in  i_TS. Only 0-3 are allowed")


    #####################
    print("   ...Checking Scatter")
    ## If Time series requested
    if IF.i_SCT!=0:

        freqs=list(set([f.split(",")[0] for f in IF.c_SCT_freq]))
        assert set(freqs).intersection(freqs_avail) == set(freqs), "   ...Unrecognized plot frequency in c_TS_freq. Check frequencies specified before the comma"
        for freq in freqs:
            if freq not in IF.c_TimeMeanRes:
                print(f"      ...WARNING: requested freq  {freq}  doesn't exist in  c_TimeMeanRes. It will be added to variable")
                IF.c_TimeMeanRes=IF.c_TimeMeanRes+[freq]

        ## Check frequencies of plots that were asked
        freqs=set([f.split(",")[1] for f in IF.c_SCT_freq])
        assert freqs.intersection(freqs_avail+["Whole"]) == freqs, "   ...Unrecognized plot frequency in c_TS_freq. Check frequencies specified after the comma"
        ## Check that plots requested exist
        plots=open_plots(IF.c_SCT_Plots)
        check_plots_vars_exist(plots)

        ## Check requested plot types exist
        if len(set([s for s in str(IF.i_SCT)]).difference(['0', '1','2','3', "4"])) != 0:
            print("      ...WARNING: You have set an unrecognized value in  i_SCT. Only 0-4 are allowed")
 
        ## Check percentiles make sense
        for perc in IF.i_SCT_perc:
            assert perc<=100,  "      ...Percentile  {perc}  requested is greater than 100. Please change to an appropriate value"
            assert perc>=0,  "      ...Percentile  {perc}  requested is lower than 0. Please change to an appropriate value"




    #####################
    print("   ...Checking Surface Maps")
    ## If Time series requested
    if IF.i_Map!=0:
        freqs=list(set([f.split(",")[0] for f in IF.c_Map_freq]))
        assert set(freqs).intersection(freqs_avail) == set(freqs), "   ...Unrecognized plot frequency in c_TS_freq. Check frequencies specified before the comma"
        for freq in freqs:
            if freq not in IF.c_TimeMeanRes:
                print(f"      ...WARNING: requested freq  {freq}  doesn't exist in  c_TimeMeanRes. It will be added to variable")
                IF.c_TimeMeanRes=IF.c_TimeMeanRes+[freq]
        ## Check frequencies of plots that were asked
        freqs=set([f.split(",")[1] for f in IF.c_Map_freq])
        assert freqs.intersection(freqs_avail+["Whole"]) == freqs, "   ...Unrecognized plot frequency in c_TS_freq. Check frequencies specified after the comma"
        ## Check that plots requested exist
        plots=open_plots(IF.c_Map_Plots)
        check_plots_vars_exist(plots)

        ## Check requested plot types exist
        if len(set([s for s in str(IF.i_Map)]).difference(['0', '1','2','3'])) != 0:
            print("      ...WARNING: You have set an unrecognized value in  i_Map. Only 0-3 are allowed")

        ## Check requested shapefiles 
        if len(set([s for s in str(IF.i_Map_Shapefiles)]).difference(['0', '1','2','3', '4', '5'])) != 0:
            print("      ...WARNING: You have set an unrecognized value for shapefile in   i_Map_Shapefiles. Only 0-5 are allowed")



    IF.c_RMap_savetype = [int(c) for c in str(IF.c_RMap_savetype)]
    IF.i_Write_net = [int(c) for c in str(IF.i_Write_net)]




    ##############################################
    #####          FINAL MESSAGE             #####
    ##############################################
    if ignore_exceptions:

        print("")
        print("...Congratulations IncludeFile without problems, but you have ignored exceptions [`-´]")
    else:
        print("")
        print("...Congratulations IncludeFile without problems [^.^]")
    
    print("")
    print("")
    print("")





##########
## Check that it has nested types
def has_nested_type(object, *types):
    """
    Returns True or False whether 'object' is a nested object with the types
    specified at 'types'.

    Important: this checks all elements in the nested object, so it should be
    used carefully to not impact performance.

    Example:

    > o1 = [[1, 2], [3, 4]]
    > has_nested_type(o1, list, list, int)
    True

    > o2 = [['1', '2'], ['3', '4'], '5']
    > has_nested_type(o1, list, list, str)
    False

    > o3 = [['1', '2'], ['3', '4'], []]
    > has_nested_type(o3, list, list, str)
    True
    """
    # Checks first container
    if not isinstance(object, types[0]):
        return False

    # Base case, only one nested level
    if len(types) == 2:
        for elem in object:
            if not isinstance(elem, types[1]):
                return False
        return True

    # More than two levels of nested objects
    elif len(types) > 2:
        return all([has_nested_type(elem, *types[1:]) for elem in object])


if __name__ == '__main__':
    # for testing purposes
    from sys import path
    path.append('src/..')
    run(ignore_exceptions=True)

    # nested object tests
    o1 = [[1, 2], [3, 4]]
    print(has_nested_type(o1, list, list, int))  # True

    o2 = [['1', '2'], ['3', '4'], '5']
    print(has_nested_type(o2, list, list, str))  # False

    o3 = [['1', '2'], ['3', '4'], []]  # True
    print(has_nested_type(o3, list, list, str))
