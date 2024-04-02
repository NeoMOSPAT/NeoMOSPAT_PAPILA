# Basic plots for cross sections. To make it simpler it will plot all
# calculated cross sections when they exists
import src.IncludeFile as IncF
from os import path, makedirs
from typing import List, Literal, Union

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from pandas import DataFrame  # only for typing at this point

from src.plot.aux_plots import XSSConfig
from src.plot import maps_utils

FORMAT = Literal['png', 'svg']

def get_savepath(model: str, figure_file_name: str='') -> str:
    return f"{IncF.c_FigureDir}/CROSS_SECTIONS/{model}/{figure_file_name}"


class CrossSectionPlotter:
    """
    Class that encapsulates all cross sections plots.

    Parameters:
        bbox: model bbox. Default to make maps
        data: Cross section data is a dictionary with the following keys:
            - 'lat'
            - 'lon'
            - 'height'
            - 'time'
            - model variables keys
        figsize: size in inches of plot's width
        model: Name or label of the model from where data was extracted.
        name: Name or label of the cross section.
        shape: ('time', 'z', 'length') where 'z' is the number of horizontal
               levels and 'length' how many points conform the cross sections.
        tight_bbox: tight bounding box to create maps. It has a margin of 5% for
                    each side respect to a tight bounding box.
    """

    def __init__(self, cross_section_data: dict):
        self.data = dict(cross_section_data)  # always inits a copy
        self.name = self.data.pop('c_name')
        self.model = self.data.pop('c_model')
        self.shape = self.data['height'].shape
        # compute bounding box
        self.tight_bbox = CrossSectionPlotter._set_maps_margin(self.data['lat'],
                                                         self.data['lon'],
                                                         percentage=25)
        # TODO: add margin and also consider tight_bbox if it goes out model bbox
        self.bbox = cross_section_data['f_bbox']
        self.figsize = XSSConfig.figheight  # default for now


    def plot_path(self, f_bbox: List[float]=[], f_surface_model_data: Union[dict, DataFrame]=None, save: bool=True) -> None:
        """
        Plots cross section path on a map using self.bbox as limits.
        Plot will be saved in the specified figure directory at IncludeFile with
        the path:

        [figure_dir]/CROSS_SECTIONS/[model_name]/[cross_section_name]_MAP.[format]

        Args:
            f_bbox: If given, overrides ignores self.bbox. Useful when plotting
            filtered regions.
            f_surface_model_data: Dict or Pandas.DataFrame that contains
            two-dimensional data array to plot a pcolormesh or fill-contour
            (Not decided yet). If dict:

            {
                'lat': numpy.ndarray,
                'lot': numpy.ndarray,
                'time': Union[List[datetime], numpy.ndarray],
                'model_var_1': numpy.ndarray,
                'model_var_2': numpy.ndarray,
                ...
                'model_var_3': numpy.ndarray
            }

            if Pandas.DataFrame, must have the same columns as the keys in the
            dict, where 'time' can be the column. TODO: contrast with t_ModelData.

        Returns:
            None
        """
        fig = plt.figure()
        width = self.figsize
        bbox = self.bbox
        fig.set_size_inches(w=width, h=width*abs((bbox[3]-bbox[2])/(bbox[1]-bbox[0])))
        # TODO: make larger side equals to self.figsize
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

        ax.set_title(f"{self.name} from model {self.model}")

        # base map
        maps_utils.create_map(bbox, width, ax=ax, gridlines=True)
        # maps_utils.plot_gridlines(lats, lons, ax=ax)
        # TODO: add pcolormesh of integrated data
        ax.plot(self.data['lon'], self.data['lat'], marker='*', color='red',
                markersize=15, transform=ccrs.PlateCarree())

        plt.tight_layout()

        # NOTE: save == False is used in dev process, please keep it
        if save:
            fname = self.get_savepath(figure_file_name=f"{self.name}_MAP")
            self.save(fname, XSSConfig.figformat)
        else:
            plt.show()

        plt.close(fig=fig)

        print(f"   ...Files {fname} has been saved") 


    def plot(self, representation: List[str]=XSSConfig.representation, save: bool=True) -> None:
        """
        Plots model variables values as scatter plots with colorbar.
        Horizontal axis shows coordinates of the sampled points and vertical
        axis the height.

        Args:
            representation: Which kind of drawing to represent the data.
            Possible values: ['pcolormesh'], ['pcolormesh', 'scatter']. Default: `XSSConfig.representation`.

        Plot will be saved in the specified figure directory at IncludeFile

        [figure_dir]/[model_name]/[cross_section_name]_[timestamp]_LEVELS.[format]

        """

        assert 'pcolormesh' in representation, "'pcolormesh' is mandatory as \
                                                cross section representation"
        print(f'...Plotting "LEVELS" figures for "{self.name}" cross section from model "{self.model}"')
        print(f'   ...{len(self.data["time"])} figures will be rendered for each variable')


        # xaxis is the same for all figures
        nrows = self.data['height'].shape[1]
        ncols = self.data['height'].shape[2]
        xaxis = np.array([list(range(ncols))]*nrows)

        for model_var in self.data:
            if model_var in ('f_bbox', 'lat', 'lon', 'time', 'height'):
                continue
            print(f"   ...Creating {model_var} figures")

            # global min and max value of the current variable
            var_min = np.nanmin(self.data[model_var])
            var_max = np.nanmax(self.data[model_var])

            # Color range
            norm = getattr(mpl.colors, XSSConfig.norm)(vmin=var_min, vmax=var_max)

            # Color map
            cmap = XSSConfig.cmap

            for time_index, time in enumerate(self.data['time']):
                # TODO: get data freq at instance creation to change param unit
                # to use those described in
                # https://numpy.org/doc/stable/reference/arrays.datetime.html#datetime-units
                time = np.datetime_as_string(time, unit='h')

                fig = plt.figure()
                # each colum of data should occupy a portion of the size
                fig.set_size_inches(w=ncols+1, h=self.figsize)
                ax = fig.add_subplot(1, 1, 1)
                # Separates what is drawn from the edges of the figure
                ax.use_sticky_edges = False

                yaxis = self.data['height'][time_index]
                colors = self.data[model_var][time_index]  # actual data

                xaxis, yaxis, colors = self.crop_nan_columns(xaxis, yaxis, colors)

                plot = ax.pcolormesh(xaxis,
                                     yaxis,
                                     colors,
                                     alpha=1,
                                     shading='gouraud',
                                     cmap=cmap,
                                     norm=norm)
                if 'scatter' in representation:
                    ax.scatter(xaxis,
                               yaxis,
                               # TODO: find a way to use inverted color map instead of inverted data
                               c=colors**-1,
                               edgecolor='white',
                               s=150,
                               cmap=cmap,
                               alpha=1,
                               norm=norm)
                ax.grid(True)  # Very important for correct visualization
                # TODO: this should be in grades
                xlabels = [f"({_lat:.1f}, {_lon:.1f})" \
                           for _lat, _lon in zip(self.data['lat'], self.data['lon'])]
                ax.set(
                    # TODO: add contaminant units
                    title=f"{self.model} {model_var}\n{self.name} {time}",
                    xlabel=XSSConfig.xlabel,
                    ylabel=XSSConfig.ylabel,
                    xticks=range(len(xlabels))
                )
                ax.set_xticklabels(xlabels, rotation=45)

                plt.colorbar(
                    plot, # plot associated with the colorbar
                    format=mpl.ticker.ScalarFormatter(), # el formato de los ticks de la barra, números reales en este caso
                    ticks=mpl.ticker.IndexLocator(3, 0)) # cada cuántos ticks poner una label y a partir de que valor

                # plt.tight_layout() -> eats some gridline's ticks
                # NOTE: save == False is used in dev process, please keep it
                if save:
                    fname = self.get_savepath(figure_file_name=f"{self.name}_{model_var}_{time}_PCOLORMESH")
                    self.save(fname, XSSConfig.figformat)
                else:
                    plt.show()

                plt.close(fig=fig)


    def set_maps_margin(self, percentage: int=5):
        """
        Change margin for maps.

        Args:
            percentage: (int) value between 1 and 100
        """
        assert 1 <= percentage <= 100, "Margin should be between 1% to 100%"
        self.bbox = CrossSectionPlotter._set_maps_margin(self.data['lat'],
                                                         self.data['lon'],
                                                         percentage=percentage)

    @classmethod
    def _set_maps_margin(cls, lats: List[float], lons: List[float], percentage: int=5) -> List[float]:
        min_lat, min_lon, max_lat, max_lon = np.min(lats), \
                                             np.min(lons), \
                                             np.max(lats), \
                                             np.max(lons)
        margin_ratio = percentage/100
        horizontal_margin = abs(max_lon - min_lon) * margin_ratio
        vertical_margin = abs(max_lat - min_lat) * margin_ratio
        # lon-east, lon-west, lat-south, lat-north
        return [
            max_lon + horizontal_margin,
            min_lon - horizontal_margin,
            max_lat + vertical_margin,
            min_lat - vertical_margin
        ]

    def get_savepath(self, figure_file_name: str='') -> str:
        return get_savepath(model=self.model, figure_file_name=figure_file_name)

    def save(self, figure_path: str, format: FORMAT):
        figure_path += '.' + format
        plt.savefig(figure_path, format=format)
        print(f"File {figure_path} has been saved")

    @staticmethod
    def crop_nan_columns(xaxis: np.ndarray, yaxis: np.ndarray, colors: np.ndarray) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Checks if data contains NaN columns. This method only cover the case
        when such a columns are the sides.

        Args:
            xaxis: `i` coordinates for `pcolormesh`
            yaxis: `j` coordinates for `pcolormesh`. Can contain NaNs.
            colors: data for `pcolormesh`. Can contain NaNs.
        """
        # Check columns form left to right and the from right to left and
        # records how many columns it should ignore at each side
        # It's enough to only check `yaxis`, a.k.a `heights` matrix.
        skip_at_left : int = 0
        skip_at_right: int = 0
        ncols = yaxis.shape[1]

        # from left to right
        for ncol in range(ncols):
            if all(np.isnan(yaxis[:, ncol])):
                skip_at_left += 1
            else:
                break

        if ncol != ncols - 1:  # edge case
            # from right to left
            for ncol in range(1, ncols + 1):
                if all(np.isnan(yaxis[:, -ncol])):
                    skip_at_right += 1
                else:
                    break

        slicer = (slice(None), slice(skip_at_left, ncols - skip_at_right))
        if xaxis.shape[1] < ncols: # Avoids `xaxis` being cropped multiple times. Not very clean
            return xaxis, yaxis[slicer], colors[slicer]
        else:
            return xaxis[slicer], yaxis[slicer], colors[slicer]

def serialize_xss_data(t_Cross_Sections):
    """ Only for debug purposes """
    import datetime as dt
    import json
    from src import serializers as ser
    xss = ser.preserialize(t_Cross_Sections)
    xss_file = open(f'dumped_cross_sections_{dt.datetime.now().strftime("%d%m%Y")}.json', 'w')
    json.dump(xss, xss_file)
    xss_file.close()

# TODO: Add correct typing to `t_ModelFilters` param
def main(t_Cross_Sections: List[dict], t_ModelFilters: dict) -> None:
    if not IncF.t_cross_sections_config:
        print("...Plotting was not requested.")
        return None
    elif not t_Cross_Sections:
        print("...No information exists to be plotted.")
        return None

    plotters : List[CrossSectionPlotter] = [CrossSectionPlotter(xss) for xss in t_Cross_Sections]

    # Make dirs for every model and cross section name. Figures dir should
    # already exist thanks to IncludeFile checking
    for plotter in plotters:
        savepath = plotter.get_savepath()
        if not path.exists(savepath):
            makedirs(savepath)

    # Make maps
    for plotter in plotters:
        plotter.plot_path()

    # Make pcolormesh full domain
    for plotter in plotters:
        plotter.plot()

    # Make pcolormesh filters

def test():
    pass

if __name__=='__main__':
    test()

