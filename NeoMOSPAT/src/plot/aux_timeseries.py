# -*- coding: utf-8 -*-
import src.IncludeFile as IncF
#
import numpy as np
import pandas as pd
#
import matplotlib.pyplot as plt
#
import src.common as aux
import src.plot.aux_plots as aux_plots
#
from os.path import join
from os import makedirs
#
import geopy.distance
from shapely.ops import nearest_points



TS_plots=[int(p) for p in str(IncF.i_TS)]
bool_stations=False
######################################
settings=aux_plots.TSConfig
##
color1, color2=settings.colors
figsize = settings.imagesize
ms = settings.markersize
ls = settings.linestyle
legend_location = settings.legend_location
nele_legend=settings.nele_legend
linewidth_highlight=settings.linewidth_highlight
linewidth_all=settings.linewidth_all
color_highlight=settings.color_highlight
color_highlight_obs=settings.color_highlight_obs
if IncF.b_TS_hlmodel:
    ls=["-o"]+ls
    color2=color_highlight[:len(IncF.c_TS_hlmodel)]+color2
if IncF.b_TS_hlobs:
    ls=["-o"]+ls
    color1=color_highlight_obs[:len(IncF.c_TS_hlobs)]+color1

titlesize = settings.titlesize
ticksize = settings.ticksize
labelsize = settings.labelsize
overlap_percentage = settings.overlap_percentage
overlap_twovar = 0.4
legendsize = settings.legendsize
b_closest = settings.closest
min_perc_window = settings.percwindow
bool_center = settings.center
b_showshared = False #settings.b_showshared
subtitlesize= settings.subtitlesize
perc_top=settings.perc_top
perc_add=(100-perc_top)/100

#####################################
def prepare_metadata(fig_meta, ax_meta, locations, clean_var_x, clean_var_y):
    ## Figure Metadata
    fig_meta=fig_meta or plt.subplots(1, figsize=figsize)
    ax_meta=ax_meta or plt.gca()

    ax_meta.set_axis_off()
    
    fig_meta.suptitle('%s vs %s'%(clean_var_x, clean_var_y), fontsize=titlesize, y=1)
    ## Titles and subtitles
    subtitle='\n'.join(map(str, locations[1:]))
    ax_meta.set_title(subtitle, fontsize=subtitlesize, y=0.975, loc="left", fontdict={'verticalalignment':"top"})


#####################################
def do_legend(fig, ax1, handles, new_legend, legend_title, plot_type="ts"):
    cols=len(handles)
    #cols=int(len(handles)/2)+1
    cols= nele_legend if cols>nele_legend else cols
    ax1.legend(handles, new_legend, loc=legend_location, framealpha=0, ncol=cols, fontsize=int(0.9*legendsize), handletextpad=0.4)#, title=legend_title, title_fontsize=legendsize
    if plot_type!="ts":
        if plot_type=="DC":
            myFmt = mdates.DateFormatter('%H')
            ax1.xaxis.set_major_formatter(myFmt)
        elif plot_type=="AC":
            myFmt = mdates.DateFormatter('%m-%d')
            ax1.xaxis.set_major_formatter(myFmt)

#####################################
def get_color_columns(columns):
    def pick_color(c, counter_obs=[0], counter_model=[0]):
        if "Obs" in c:
            counter_obs[0]+=1
            return color1[counter_obs[0]]
        elif "Model" in c:
            counter_model[0]+=1
            return color2[counter_model[0]]
    counter_obs=[-1]
    counter_model=[-1]
    color_dict = {c: pick_color(c, counter_obs=counter_obs, counter_model=counter_model) for c in columns}
    return color_dict


######################################################
## Finds the closest grid points of stations in df to the models
def find_closest_points(info1, info2):
    #########
    def near(point, pts=info2.geometry.unary_union):
        nearest = info2.geometry == nearest_points(point, pts)[1]
        return info2[nearest].ID.values[0], info2[nearest].geometry.values[0]

    ## First nearest neighbor informations
    info1[['ID_Nearest', "Nearest"]] = pd.DataFrame(info1.apply(lambda row: near(row.geometry), axis=1).values.tolist(), index=info1.index)
    ## Distance
    info1["Distance"]=info1.apply(lambda a: geopy.distance.geodesic((a.geometry.y, a.geometry.x), (a.Nearest.y, a.Nearest.x)), axis=1)

    return info1


#####################################
def prepare_labels(columns, columns_1, columns_2, var1, var2, units):
    ##
    model=list(set([c.split()[1] for c in columns if "Model" in c]))
    network=list(set([c.split()[1] for c in columns if "Obs" in c]))
    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    evar1=[k for  k, v in variables_equivalences.items() if var1 in v] 
    evar1=evar1[0] if evar1 else var1
    evar2=[k for  k, v in variables_equivalences.items() if var2 in v] 
    evar2=evar2[0] if evar2 else var2
    
    ## Dont repeat station because it will always be in reference to it
    if len(network)==1:
        if len(model)!=0:
            ##
            new_legend=[c.replace("Obs ", "").replace("Model ", "") for c in columns_1]+ [c.replace("Obs ", "").replace("Model ", "") for c in columns_2]
            ## Removes var from legend
            new_legend=[" ".join(c.split(" ")[:-1]) for c in new_legend]
            new_legend=[c.replace("PAPILA-", "") for c in new_legend]
            ## Remove _1 _2
            #new_legend=["_".join(c.split("_")[:-1]) for c in new_legend]
            new_legend=["_".join(c.split("_")[:-1]) if "_" in c else c for c in new_legend]
        else:
            ## If no models then only show ID of station
            new_legend=[c.split()[-1] for c in columns_1]+ [c.split()[-1] for c in columns_2]

        
        y1label="%s [%s]"%(evar1, aux_plots.get_pretty_units(units[var1]))
        y2label="%s [%s]"%(evar2, aux_plots.get_pretty_units(units[var2]))

        return y1label, y2label, new_legend

    elif len(network)==2:
        ##
        new_legend=[c.replace("Obs ", "").replace("Model ", "") for c in columns_1]+ [c.replace("Obs ", "").replace("Model ", "") for c in columns_2]
        ## Removes var from legend
        new_legend=[" ".join(c.split(" ")[:-1]) for c in new_legend]

        units_1=units[var1]#aux.read_units(columns_1[0].split()[1])
        units_2=units[var2]#aux.read_units(columns_2[0].split()[1])
        y1label="%s [%s]"%(evar1, aux_plots.get_pretty_units(units[var1]))
        y2label="%s [%s]"%(evar2, aux_plots.get_pretty_units(units[var2]))
        return y1label, y2label, new_legend 

    elif len(network)==0:
        new_legend=[c.replace("Obs ", "").replace("Model ", "") for c in columns_1]+ [c.replace("Obs ", "").replace("Model ", "") for c in columns_2]
        new_legend=[c.replace(var1, "").replace(var2, "") for c in new_legend]

        #units=aux.read_units("SINCA")
        y1label="%s [%s]"%(evar1, aux_plots.get_pretty_units(units[var1]))
        y2label="%s [%s]"%(evar2, aux_plots.get_pretty_units(units[var2]))
        return y1label.replace("[]", "").replace("[ad]", ""), y2label.replace("[]", "").replace("[ad]", ""), new_legend     

#####################################
def highlight_first_line(ax, labels, names, case=""):
    if case=="Model":
        colorhl=color_highlight
    else:
        colorhl=color_highlight_obs

    for index in names:
        col=labels[index]
        color=colorhl[index]
        for line in ax.get_lines():
            if line.get_label() == col:
                line.set_linewidth(linewidth_highlight) 
                line.set_linestyle("-")
                line.set_color(color)
                line.set_zorder(10000)

#####################################
## Add extra padding to xaxis so that the plot looks prettier
def align_xaxis(ax1):
    xmin, xmax = ax1.get_xlim()
    rangex = xmax - xmin
    #print(xmin, xmax)
    ax1.set_xlim(xmin-0.025*rangex, xmax+0.025*rangex)

#####################################
## Align y axis
def align_yaxis(ax1, ax2, BoolTwoVars=True, var1="", var2="", nyticks=7):

    if ax2 is None:
        #ax2=ax1
        ax2 = ax1.twinx()  

    
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.patch.set_visible(False)
    ax2.patch.set_visible(True)
       
    ############
    ax1.xaxis.set_tick_params(which='both', labelsize=ticksize)
    ax1.yaxis.set_tick_params(which='both', labelsize=ticksize)
    ax2.xaxis.set_tick_params(which='both', labelsize=ticksize)
    ax2.yaxis.set_tick_params(which='both', labelsize=ticksize)

    ax1.xaxis.label.set_size(labelsize)
    ax1.yaxis.label.set_size(labelsize)
    ax2.xaxis.label.set_size(labelsize)
    ax2.yaxis.label.set_size(labelsize)


    ax1.axes.get_xaxis().get_label().set_visible(False)
    ax2.axes.get_xaxis().get_label().set_visible(False)

    #############
    lines1=ax1.lines
    lines2=ax2.lines
    try:
        ## Works if ax1 and ax2 have data,  or if ax1 has data and ax2 doesnt
        y1=np.concatenate([l.get_ydata() for l in lines1], axis=0)
        y2=np.concatenate([l.get_ydata() for l in lines2], axis=0) if len(lines2)!=0 else y1
    
    except ValueError:
        ## Works If ax2 has data and ax1 doesnt
        y2=np.concatenate([l.get_ydata() for l in lines2], axis=0)
        y1=np.concatenate([l.get_ydata() for l in lines1], axis=0) if len(lines1)!=0 else y2

    ##
    per_max=100 if IncF.i_TS_ShowAll==1 else 99.5
    per_min=0 if IncF.i_TS_ShowAll==1 else 0.5
    min1, max1=np.nanpercentile(y1, per_min), np.nanpercentile(y1, per_max)
    min2, max2=np.nanpercentile(y2, per_min), np.nanpercentile(y2, per_max)




    ##
    overmax = min(max1, max2) 
    overmin = max(min1, min2)
    
    overlap = 0 if overmax-overmin<0 else overmax-overmin
    maxg = max(max1, max2)
    ming = min(min1, min2)
    range_data = maxg-ming 
    range_data = range_data if range_data!=0 else 0.1



    variables_equivalences=aux.read_yaml("vars")["Equiv"]
    equiv1=[k for  k, v in variables_equivalences.items() if var1 in v][0]
    equiv2=[k for  k, v in variables_equivalences.items() if var1 in v][0]
    
    #"""
    if equiv1==equiv2:
        BoolTwoVars=False
    else:
        BoolTwoVars=True
    #"""


    if IncF.i_TS_SharedYRange and BoolTwoVars:
        print("One vars:  1")
        if not np.isnan(ming) and not np.isnan(maxg):
            ax1.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax2.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax1.set_yticks(np.linspace(ming, maxg+perc_add*range_data, nyticks))
            ax2.set_yticks(np.linspace(ming, maxg+perc_add*range_data, nyticks))
            ax1.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            ax2.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
        return True      


    if IncF.i_TS_SharedYRange and not BoolTwoVars:
        print("One vars:  2")
        if not np.isnan(ming) and not np.isnan(maxg):
            ax1.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax2.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax1.set_yticks(np.linspace(ming, maxg+0.05*range_data, nyticks))
            ax1.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            #ax2.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            ax2.set_yticks([])
        ax1.yaxis.label.set_color('black')
        ax1.tick_params(axis='y', colors='black')
        ax2.set_ylabel("")
        return True     


    ## One axis for one var if there is overlap of at least 40% of data 
    if  not BoolTwoVars and (overlap!=0 and overlap/range_data>=overlap_percentage):
        print("One vars:  3")
        if not np.isnan(ming) and not np.isnan(maxg):
            ax1.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax2.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax1.set_yticks(np.linspace(ming, maxg+0.05*range_data, nyticks))
            ax1.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            #ax2.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            ax2.set_yticks([])
        ax1.yaxis.label.set_color('black')
        ax1.tick_params(axis='y', colors='black')
        ax2.set_ylabel("")
        return True

    ## One axis ticks for two vars if there is overlap of at least 70% of data
    if  BoolTwoVars and (overlap!=0 and overlap/range_data>=overlap_twovar):
        print("One vars:  4")
        if not np.isnan(ming) and not np.isnan(maxg):
            ax1.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax2.set_ylim(ming-0.05*range_data, maxg+perc_add*range_data)
            ax1.set_yticks(np.linspace(ming, maxg+0.05*range_data, nyticks))
            ax2.set_yticks(np.linspace(ming, maxg+0.05*range_data, nyticks))
            ax1.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            ax2.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
        return True

    ## Two different axis plot for two vars or one var if there is no overlap
    else:
        print("Two vars")
        ##
        range1=max1-min1 
        range2=max2-min2
        if not np.isnan(range1):
            ax1.set_ylim(min1-0.05*range1, max1+perc_add*range1)
            ## Looks prettier
            #ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
            ## Ticks aligned in both axis
            ax1.set_yticks(np.linspace(min1, max1+0.05*range1, nyticks))
            ax1.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            #ax2.yaxis.set_major_locator(plt.MaxNLocator(nyticks))

        if not np.isnan(range2):
            ax2.set_ylim(min2-0.05*range2, max2+perc_add*range2)
            ## Looks prettier
            #ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
            ## Ticks aligned in both axis
            ax2.set_yticks(np.linspace(min2, max2+0.05*range2, nyticks))
            #ax1.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
            ax2.yaxis.set_major_locator(plt.MaxNLocator(nyticks))
        return False


#####################################
def get_new_data(vars_to_plot, t_ObsData, t_ModData, plot_freq):
    BoolObs=False
    BoolMod=False

    ## Merge obs data 
    try:
        obs_data={f"Obs {group} {var}":t_ObsData[group][var] for group in vars_to_plot["O"] for var in vars_to_plot["O"][group]}
        #obs_units={aux.clean_network(group):{var:t_ObsData[group].meta.units[var] for var in vars_to_plot["O"][group]} for group in vars_to_plot["O"] }
        obs_units={var:t_ObsData[group].meta.units[var] for group in vars_to_plot["O"] for var in vars_to_plot["O"][group]}
        obs_data=pd.concat(list(obs_data.values()), keys=list(obs_data.keys()), axis=1 )
        obs_data.meta.units=obs_units
        #obs_data.dropna(how="all", axis=1, inplace=True)
        BoolObs=True

    except:
        obs_data=pd.DataFrame()
        obs_units={}
    
    ## Merge model data 
    try:
        mod_data={f"Model {model} {var}":t_ModData[model][group][var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model] for group in vars_to_plot["O"]}
        mod_units={var:t_ModData[model][group].meta.units[var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model] for group in vars_to_plot["O"]}
        mod_data=pd.concat(list(mod_data.values()), keys=list(mod_data.keys()), axis=1)
        mod_data.meta.units=mod_units
        #mod_data.dropna(how="all", axis=1, inplace=True)
        BoolMod=True
    except:
    
        try:
            group=list(t_ModData[list(t_ModData.keys())[0]].keys())[0]
            mod_data={f"Model {model} {var}":t_ModData[model][group][var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model]}
            mod_data={k:v.loc[:, v.columns.str.contains('^f[0-9]-+')] for k,v in mod_data.items() }
            mod_data=pd.concat(list(mod_data.values()), keys=list(mod_data.keys()), axis=1)
            mod_units={var:t_ModData[model][group].meta.units[var] for model in vars_to_plot["M"] for var in vars_to_plot["M"][model]}
            BoolMod=True
        except:
            mod_data=pd.DataFrame()
            mod_units={}


    if BoolObs and BoolMod:
        if IncF.i_TS_Sum_SharedXRange==1:
            idx_common = obs_data.index.intersection(mod_data.index)
            mod_data.loc[mod_data.index.difference(idx_common)]=np.nan


        ## Units
        mod_data.meta.units=mod_units
        obs_data.meta.units=obs_units

    ## Delete
    #obs_data=obs_data.applymap(lambda l: l if not np.isnan(l) else np.random.choice([1, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5])) 2
    data=pd.concat([obs_data, mod_data], axis=1)#.dropna()

    data.meta.units={**obs_units, **mod_units}#{"O":obs_units, "M":mod_units}
    return data, obs_data, mod_data
    

#####################################
def prepare_data_one_var(data, color_dict):
    ## Separate data
    data_1=data[[v for v in data.columns if "Obs" in v]]
    data_2=data[[v for v in data.columns if "Model" in v]]
    
    if data_1.empty and not data_2.empty:
        case="models"
        color1_0="black"
        color2_0=color_dict[data_2.columns[0]]
    elif not data_1.empty and not data_2.empty:
        case= "model-obs"
        color1_0=color_dict[data_1.columns[0]]
        color2_0=color_dict[data_2.columns[0]]
        if data_2.dropna(thresh=1).empty or data_1.dropna(thresh=1).empty:
            return None, None, None, None, None
    elif data_2.empty and not data_1.empty: 
        networks=sorted(list(set([c.split()[1] for c in data.columns])))
        try:
            case= "obs-obs"
            data_1=data[[c for c in data.columns if networks[0] in c]]
            data_2=data[[c for c in data.columns if networks[1] in c]]
            color1_0=color_dict[data_1.columns[0]]
            color2_0=color_dict[data_2.columns[0]]
            if data_2.dropna(thresh=1).empty or data_2.dropna(thresh=1).empty:
                return None, None, None, None, None
        except:
            case= "obs"
            data_1=data
            data_2=pd.DataFrame()

            color1_0=color_dict[data_1.columns[0]]
            color2_0="black"
    else:
        return None, None, None, None, None

    return case, color1_0, color2_0,  data_1, data_2


#####################################
def get_model(df, freq, group):
    data = df.reset_index(level=["time"])
    data=data.loc[group].groupby("time").mean()
    var=data.columns[0] 
    return data.resample(freq).mean()[var]
