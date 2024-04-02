# Common functions to aid map plots

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# DEFAULT_GRIDLINER_RESOLUTION = 0.5

def create_map(bbox: list, width_inches: int, ax=None, gridlines=False):
    """
    Creates base map using Cartopy with PlateCarree projection.

    Args:
        bbox: decimal coordinates defining the map's bounding box. Its format is
        [lon-east, lon-west, lat-south, lat-north]
        size_inches: Width size of the map. Height will be adjusted to
        keep relation aspect ratio with 'bbox'.
        ax: Matplotlib Axis instance. Default 'None'. If given the map will
        be drawn on it.

    Returns:
        fig, ax: Figure and Axes instance
    """
    if ax is None:
        fig = plt.figure()
        lon_size = bbox[1]-bbox[0]
        lat_size = bbox[3]-bbox[2]
        size = (width_inches, width_inches*abs(lat_size/lon_size))
        fig.set_size_inches(*size)
        # print(size)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    else:
        fig = plt.gcf()  # get current figure

    ax.set_extent(bbox, crs=ccrs.PlateCarree())  # lon-east, lon-west, lat-south, lat-north

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, alpha=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.STATES)

    if gridlines:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=1.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

    return fig, ax

def plot_gridlines(lats: list, lons: list, bbox=None, ax=None):
    """
    Generates all horizontal and verticales lines of given mesh on the map

    Args:
        lats: Two dimesional Numpy.ndarray with latitudes mesh.
        lons: Two dimesional Numpy.ndarray with longitudes mesh.
        bbox: *NOT IMPLEMENTED* Default 'None'. If given it will optimize
        drawings to that bounding box
        ax: Matplotlib Axes instance. Default 'None'. If given gridlines
        will be plotted on it. Make sure it has a map drawn.

    Returns:
        ax: Axes instance where gridlines are plotted.


    """
    # TODO: if bbox solo dibujar la lineas entre medio
    horizontals = []
    for i in range(lats.shape[0]):
        # TODO: use named tuples
        horizontals.append(dict(lat=None, lon=None))
        horizontals[i]['lat'] = lats[i, :]
        horizontals[i]['lon'] = lons[i, :]
    for h in horizontals:
        plot_kwargs = dict(x=h['lon'], y=h['lat'], c='grey', linestyle='--', transform=ccrs.PlateCarree())
        if ax is None:
            plt.plot(**plot_kwargs)
        else:
            ax.plot(**plot_kwargs)

    verticals = []
    for i in range(lats.shape[1]):
        # TODO: use named tuples
        verticals.append(dict(lat=None, lon=None))
        verticals[i]['lat'] = lats[:, i]
        verticals[i]['lon'] = lons[:, i]
    for v in verticals:
        plot_kwargs = dict(x=v['lon'], y=v['lat'], c='grey', linestyle='--', transform=ccrs.PlateCarree())
        if ax is None:
            plt.plot(**plot_kwargs)
        else:
            ax.plot(**plot_kwargs)

    return ax

def test():
    bbox = [-49.5145919, -73.6160657, -37.8251535, -28.9861358]
    fig, ax = create_map(bbox, 8, gridlines=True)
    plt.savefig('map_test.png')
    plt.show()

if __name__=='__main__':
    test()