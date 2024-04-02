# Functions to compute horizontals (not yet) and verticals planes over 3D model
# data. Only supports grids where latitude and longitude are 2D and not single
# array. Although last case is easier, we prioritize harder one xD.
# Check out examples in "planes_on_model_data" folder to know more.
# The abbreviation 'xss' is used a lot to refer to cross section(s)


import itertools
from typing import List, Tuple, Sequence

import numpy as np
import pyproj

from scipy import interpolate

from src.common import argmin, get_model_bbox

from src.IncludeFile import t_cross_sections_config, c_ModelNames
__all__ = [
    "argmin",
    "closest_idx_in_grid",
    "cross_section",
    "geodesic",
    "main"
]
# TODO: add typing

# global var
WGS84 = pyproj.Geod(ellps="WGS84")  # crs='EPSG:4326'

def sort_xss_configs_by_model(c_model_names: List[str], t_xss_config: dict) -> dict:
    """
    Split cross sections by models. That way when a model is received by the
    main routine redundant computations are avoided. A dict with model names as
    keys and a list of cross sections configs indexes as values is returned.

    Args:
        c_model_names: List[String] Names of the models beeing read by NeoMOSPAT
        t_xss_config: List[Dict[String]] Configuration for cross sections.

    Returns:
        config_by_models: dictionary as described above and in the example.

    Example:
    Consider this configuration in the IncludeFile

    c_ModelNames = ['model1', 'model2', 'model3']

    t_cross_sections_config = [
        {
            "path": [...]
        },
        {
            "path": [...],
            "models": [0]
        },
        {
            "path": [...],
            "models": [0, 2]
        },
    ]

    Then is transformed to this:

    config_by_models = {
        "model0": [0, 1, 2]
        "model1": [0]
        "model2": [0, 2]
    }
    """

    out = dict.fromkeys(c_model_names)
    # init dict values, since dict.fromkeys(c_model_names, []) don't work
    # because all keys will be pointing the same list instance
    for key in out:
        out[key] = []

    for idx_config, config in enumerate(t_xss_config):
        models = config.get('models', None)
        if models:
            # cross sections must be calculated only for these models
            _models_names = [c_model_names[idx] for idx in models]
            for name in _models_names:
                out[name].append(idx_config)
        else:
            # cross sections must be calculated for all models
            for key in out:
                out[key].append(idx_config)

    return out

def relative_position(point: Tuple[float, float], reference_point: Tuple[float, float]) -> Tuple[str, str]:
    """
    Indicates wether 'point' is at NW, NE, SW or SE from 'reference_point'.

    If 'reference_point' is equal to 'point' an empty string is returned.

    Args:
        point: lat lon sequence
        reference_point: lat lon sequence

    Returns:
        direction_SN, direction_WE: Returns 'S', 'N' or ''
        (south, north or empty string) for 'direction_SN' and 'W', 'E' or ''
        (west, east or empty string) for 'direction_WE'.

    Example:

    > arica_pos = -18.4746, -70.29792
    > rio_pos = -22.90642, -43.18223
    > relative_position(arica_pos, rio_pos)
    ('N', 'W')  # 'Arica' it's at north west of 'Rio de Janeiro'
    """
    direction_SN = ''
    if point[0] < reference_point[0]:
        direction_SN = 'S'
    elif point[0] > reference_point[0]:
        direction_SN = 'N'

    direction_WE = ''
    if point[1] < reference_point[1]:
        direction_WE = 'W'
    elif point[1] > reference_point[1]:
        direction_WE = 'E'

    return direction_SN, direction_WE


# TODO: make tests and move to common.py
def object_id_cache(function):
    """
    Similar to lru_cache but works hashing the args' object ids.

    Check https://docs.python.org/3/library/functools.html#functools.lru_cache

    It expects that decorated function only receives objects.
    """
    # TODO: add all functions with this decorator a to global dictionary to
    # perform garbage (i.e. del calls_cache) collection at some point in
    # NeoMOSPAT
    if not hasattr(function, 'calls_cache'):
        function.calls_cache = {}
    # only for *args. **kwargs not supported yet
    def wrapper(*args):
        args_ids = tuple(list(map(id, args)))
        dictkey = args_ids
        if dictkey not in function.calls_cache:
            ret = function(*args)
            function.calls_cache[dictkey] = ret
        else:
            # print('Returning cached value')
            ret = function.calls_cache[dictkey]
        return ret

    return wrapper


@object_id_cache
def calc_deltas(lats_arr, lons_arr):
    """
    Returns maximum variation (resolution) at directions south-north and
    west-east.

    Args:
        lats_arr: float two-dimensional array of latitudes
        lons_arr: float two-dimensional array of longitudes

    Returns:
        (vertical_resolution, horizontal_resolution) Tuple(float, float)

    """
    return np.max(np.diff(lats_arr, axis=0)), np.max(np.diff(lons_arr))

def distance(p1, p2):
    """
    Distance in meters between two points in WGS84 coordinate system using
    PyProj.

    Args:
        p1: Seq(float, float), first point in lat lon order
        p2: Seq(float, float), second point in lat lon order

    """
    lat1, lon1 = p1
    lat2, lon2 = p2
    _,_, d = WGS84.inv(lon1, lat1, lon2, lat2)  # distance in meters
    return d


# @object_id_cache  --> for some reason this generates a bug in rgua-fcfm
# example
def closest_idx_in_grid(lats_arr, lons_arr, point):
    """
    Calculates closest index of the coordinate for a given grid.

    lats_arr and lons_arr are two meshes that define a latlon grid (WGS84). This
    performs the calculation following these steps:

    1. Gets maximum resolution for each array

    2. Using that looks for a vicinity in order to minimize searching
    space.

    Lets say 'o' is the real point we want to approximate in the grid, where
    each 'x' represents a point in the grid.

    In the index space the algorithm will select:

    x    x    x    x    x    x                 x    x    x   |x    x|   x
                                                             |      |
                                              ---------------|------|-----
    x    x    x    x    x    x                 x    x    x   |x    x|   x
                                                             |      |
                     o                                       |  o   |
    x    x    x    x    x    x     ======>     x    x    x   |x    x|   x
                                              ---------------|------|-----
                                                             |      |
    x    x    x    x    x    x                 x    x    x   |x    x|   x
                                                             |      |
                                                             |      |
    x    x    x    x    x    x                 x    x    x   |x    x|   x
                                                             |      |
                                                             |      |
    x    x    x    x    x    x                 x    x    x   |x    x|   x

    3. Only one of those four points[1] (according to the drawing) is the
    closest point to 'o'. So, using the appropiate distance criteria, it
    calculates each distance 'o' - 'x'.

    4. Select point with less distance and return it.

    Args:

        lats_arr: ndarray. Latitudes mesh
        lons_arr: ndarray. Longitudes mesh
        point: Tuple(float, float): lat, lon point to find closest.

    Returns:
        indexes: Tuple(int, int). Row and column of closest point

    Example:

    > lat = np.array([[-37.38575745, -37.38787079, -37.3889389 , -37.3889389 ],
                      [-37.11677933, -37.11888123, -37.11994171, -37.11994171],
                      [-36.84789276, -36.84999847, -36.85104752, -36.85104752],
                      [-36.57912064, -36.58119965, -36.5822525 , -36.5822525 ]])

    > lon = np.array([[-71.52651978, -71.18792725, -70.8493042 , -70.51068115],
                      [-71.52319336, -71.18591309, -70.84863281, -70.51135254],
                      [-71.51989746, -71.18395996, -70.84799194, -70.51202393],
                      [-71.51663208, -71.18200684, -70.84732056, -70.51266479]])

    > i, j = closest_idx_in_grid(lat, lon, (-36.9, -70.7))
    > (i, j)
    (2, 2)
    > print(lat[i, j], lon[i, j])
    -36.85104752 -70.84799194

    -------------

    [1]: Four points is the easiest case to visualize, but depending the
    requested point, the 'windows' sizes and grid size can vary.
    So, it's better to think about a neighborhood of points.

    """
    mg = closest_coords_idxs_candidates(lats_arr, lons_arr, point)
    closest_coords_idxsv2 = list(zip(mg[0].flatten(), mg[1].flatten()))

    if len(closest_coords_idxsv2) == 1:
        return closest_coords_idxsv2[0]

    elif len(closest_coords_idxsv2) > 1:
        closest_coords = list(zip(lats_arr[mg].flatten(), lons_arr[mg].flatten()))
        distances = list(map(lambda p: distance((point[0], point[1]), p),
                             closest_coords))
        return closest_coords_idxsv2[np.argmin(distances)]


def closest_coords_idxs_candidates(lats_arr, lons_arr, point):
    """
    Returns all neighbours' indexes of 'point'. These neighbours are those close
    enough that get into the window around the point.

    Args:
        lats_arr: float two-dimensional array of latitudes
        lons_arr: float two-dimensional array of longitudes
        point: Tuple(float, float), lat lon point

    Returns:
        Tuple of rows and columns indexes.

    Example:

    > lat = np.array([[-37.38575745, -37.38787079, -37.3889389 , -37.3889389 ],
                      [-37.11677933, -37.11888123, -37.11994171, -37.11994171],
                      [-36.84789276, -36.84999847, -36.85104752, -36.85104752],
                      [-36.57912064, -36.58119965, -36.5822525 , -36.5822525 ]])

    > lon = np.array([[-71.52651978, -71.18792725, -70.8493042 , -70.51068115],
                      [-71.52319336, -71.18591309, -70.84863281, -70.51135254],
                      [-71.51989746, -71.18395996, -70.84799194, -70.51202393],
                      [-71.51663208, -71.18200684, -70.84732056, -70.51266479]])

    > meshgrid = closest_coords_idxs_candidates(lat, lon, (-36.9, -70.7))
    > meshgrid
    (array([[2]]), array([[2]]))
    > print(lat[meshgrid], lat[meshgrid])
    -36.85104752 -70.84799194

    """
    plat, plon = point
    delta_lat, delta_lon = calc_deltas(lats_arr, lons_arr)

    # windows
    half_delta_lat = delta_lat/2
    half_delta_lon = delta_lon/2
    winlat = (plat - half_delta_lat,  plat + half_delta_lat)
    winlon = (plon - half_delta_lon,  plon + half_delta_lon)

    # what to do if any these are empty? To raise an exception seems just right
    idxs_lats = np.where(np.logical_and(lats_arr > winlat[0],
                                        lats_arr < winlat[1]))
    idxs_lons = np.where(np.logical_and(lons_arr > winlon[0],
                                        lons_arr < winlon[1]))

    closest_coords_idxsv1 = (list(set(idxs_lats[0]).intersection(idxs_lons[0])),
                             list(set(idxs_lats[1]).intersection(idxs_lons[1])))

    mg = tuple(np.meshgrid(*closest_coords_idxsv1))

    return mg

def geodesic(point1, point2, i_points_in_between, order='latlon', sort=True):
    """
    Calculates geodesic path between two lat lon coordinates. It generates a
    geodesic of lenght 2 + i_points_in_between. Geodesic is calcuted using
    WGS84.

    Args:
        point1: Tuple(float, float), start lat lon point
        point2: Tuple(float, float), end lat lon point
        i_points_in_between: int, number of points generated besides start and
        end points.
        order: String, 'latlon' or 'lotlan' order for coordinates. Default
        'latlon'
        sort: Boolean, default True. Sorts generated path from west to east if
        the path angle is between 0° and 45° or from south to north if angle is
        between >45° and 90°

    Returns:
        geodesic: List[Tuple(float, float)]. Generated geodesic path.
    """
    y1, x1 = point1
    y2, x2 = point2

    if sort:
        inclination_angle = np.abs(np.arctan((y2-y1)/(x2-x1))*(180/np.pi))

        if 0 <= inclination_angle <= 45:
            if x1 < x2:
                lonlat_in_between_geodesic = WGS84.npts(x1, y1, x2, y2,
                                                        npts=i_points_in_between)
                pstart = point1
                pend = point2
            else:
                lonlat_in_between_geodesic = WGS84.npts(x2, y2, x1, y1,
                                                        npts=i_points_in_between)
                pstart = point2
                pend = point1
        elif 45 < inclination_angle <= 90:
            if y1 < y2:
                lonlat_in_between_geodesic = WGS84.npts(x1, y1, x2, y2,
                                                        npts=i_points_in_between)
                pstart = point1
                pend = point2
            else:
                lonlat_in_between_geodesic = WGS84.npts(x2, y2, x1, y1,
                                                        npts=i_points_in_between)
                pstart = point2
                pend = point1
    else:
        lonlat_in_between_geodesic = WGS84.npts(x1, y1, x2, y2,
                                                npts=i_points_in_between)
        pstart = point1
        pend = point2

    if order=='latlon':
        latlon_in_between_geodesic = [(y, x) for x, y in lonlat_in_between_geodesic]
        return [pstart] + latlon_in_between_geodesic + [pend]
    elif order=='lonlat':
        return [(pstart[1], pstart[0])] + lonlat_in_between_geodesic + \
               [(pend[1], pend[0])]

def i_cross_section(lats_arr, lons_arr, point1, point2):
    """
    Returns two lists of indexes, one for rows and another for cols that define
    the cross section using "closest points method".

    It calculates the minimal quantity points the geodesic needs to
    create a path between 'point1' and 'point2'. Then for each of these points
    the algorithm gets the closest point in the grid to it. The step of getting
    that minimal helps to decrease computation time, as if more points where
    placed along the geodesic some of them would share the closest point in the
    grid.

    Check out "vertical_planes" notebook at examples.
    """

    closest_p1, closest_p2 = (closest_idx_in_grid(lats_arr, lons_arr, point1),
                              closest_idx_in_grid(lats_arr, lons_arr, point2))
    # gets number of points for geodesic
    y1, x1 = closest_p1
    y2, x2 = closest_p2

    # NOTE(xZevalx): During the development of 'two_points_interpolation' I
    # figured out that npts can be calculated by kwowing the inclination of the
    # cross section and intersections with horizontals or verticals. Before, I
    # used pithagoras and that could lead to generate unneeded points.
    inclination_angle = np.abs(np.arctan((y1-y2)/(x1-x2))*(180/np.pi))
    # graphically npts can also be viewed as the number of spaces between each
    # consecutive point in the final geodesic
    npts = (np.abs(x2 - x1) if 0 <= inclination_angle <= 45 \
                            else np.abs(y2 - y1)) - 1
    # TODO: add warning if npts < 0
    geodesic_path = geodesic(point1=point1, point2=point2,
                             i_points_in_between=npts, order='latlon',
                             sort=True)
    # adjust geodesic_path to the grid
    closest_idxs_to_geodesic = list(
        map(lambda p: closest_idx_in_grid(lats_arr, lons_arr, p),
            geodesic_path))

    # Transform idxs_geodesic to allow easier indexing
    ii = []
    jj = []
    for y,x in closest_idxs_to_geodesic:
        ii.append(y)
        jj.append(x)

    return ii, jj

def cross_section(model_ds, point1, point2, variables, include_coordinates=True,
                  method='four_points_interpolation'):
    """
    Get cross section for named variables in model_ds.

    A cross section is a vertical view, like a wall or profiles collection
    between two coordinates. The number of generated points in between or, in
    other words, the length of the cross section varies depending the method you
    choose.
    In 'closest_point' method these quantity is equal as the number of
    'cells' that separate 'point1' and 'point2'.
    In 'two_points_interpolation' will be the number of verticals/horizontals if
    'point1' and 'point2' form a roughly horizontal/vertical line.
    And finally, for 'four_points_interpolation' the quantity is such that all
    points in the geodesic are equidistant and that distance is at most the
    model resolution, because, getting more points than that will result in a
    poor/short cross section.

    Args:
        model_ds: xarray.Dataset like the one you got from `read_models.nc_reader`
        point1: Tuple(float, float), start lat lon point
        point2: Tuple(float, float), end lat lon point
        variables: Seq(String), names of variables that get cross section
        include_coordinates: Boolean, default True. Whether or not to include
        latitude, longitude and height data.
        methods: String, one of the following 'two_points_interpolation',
        'four_points_interpolation', 'closest_point'. Default
        'four_points_interpolation'.

    Returns:
        Dictionary with the keys 'variables' and 'lat', 'lon', 'height' if
        'include_coordinates' was set True.
    """
    planes = None
    if method == 'closest_point':
        planes = closest_point_cross_section(model_ds, point1, point2,
                                             variables, include_coordinates)
    elif method == 'two_points_interpolation':
        planes = two_points_interpolation_cross_section(model_ds, point1,
                                                        point2, variables,
                                                        include_coordinates)
    elif method == 'four_points_interpolation':
        planes = four_points_interpolation_cross_section(model_ds, point1,
                                                         point2, variables,
                                                         include_coordinates)
    # TODO: add routine to represent as Pandas.DataFrame or
    # GeoPandas.GeoDataFrame
    planes['time'] = model_ds.time.values

    return planes

def closest_point_cross_section(model_ds, point1, point2, variables,
                                include_coordinates=True):
    """
    Gets cross sections using existing points in the grid.
    """
    lats_arr = model_ds.lat.values
    lons_arr = model_ds.lon.values
    ii, jj = i_cross_section(lats_arr, lons_arr, point1, point2)
    planes = {}
    is3D = len(model_ds.dims) == 4  # (time, height, lat, lon)
    plane_slicer = None
    if is3D:
        plane_slicer = (slice(None), slice(None), ii, jj)
    else:
        plane_slicer = (slice(None), ii, jj)
    for variable in variables:
        planes[variable] = model_ds[variable].values[plane_slicer]

    if include_coordinates:
        planes['lat'] = lats_arr[ii, jj]
        planes['lon'] = lons_arr[ii, jj]
        if is3D:
            planes['height'] = model_ds.height.values[plane_slicer]

    return planes

def two_points_interpolation_cross_section(model_ds, point1, point2, variables,
                                           include_coordinates=True):
    """
    Gets cross section using interpolation on the verticals or horizontals as
    corresponds.
    """
    # 1. Calculate angle between start and end point
    # TODO: make helper func to calc inclination angle
    inclination_angle = np.abs(
                        np.arctan((point1[0]-point2[0])/(point1[1]-point2[1]))) \
                        *(180/np.pi)

    lats_arr = model_ds.lat.values
    lons_arr = model_ds.lon.values

    # 2. Get neighbours to interpolate later (it's finished below)
    ii1, jj1 = get_neighbours_indexes(lats_arr, lons_arr, point1)
    ii2, jj2 = get_neighbours_indexes(lats_arr, lons_arr, point2)


    start_end_rel_pos = relative_position(point1, point2)
    # this defines an index bbox to the cross section
    left_limit = None
    right_limit = None
    bottom_limit = None
    top_limit = None
    # min and max are taken in order to avoid asking if closest points are in
    # the edges
    # graphically it gets a bounding box using existing grid points as follows:
    #
    #   x    x    x    x    x    x                 x   |x    x    x    x|   x
    #                                                  |                |
    #                                             -----|----------------|-----
    #   x    x    x    x    x    x                 x   |x    x    x    x|   x
    #                                                  |                |
    #                    P1                            |            P1  |
    #   x    x    x    x    x    x     ======>     x   |x    x    x    x|   x
    #                                                  |                |
    #                                                  |                |
    #   x    x    x    x    x    x                 x   |x    x    x    x|   x
    #                                                  |                |
    #          P2                                      |  P2            |
    #   x    x    x    x    x    x                 x   |x    x    x    x|   x
    #                                             -----|----------------|-----
    #                                                  |                |
    #   x    x    x    x    x    x                 x   |x    x    x    x|   x
    #
    # The purpose of this bbox is to get all cells that the geodesic pass
    # through. This is achieved in the geodesic generation step by creating one
    # more point than the "cell distance" between start and end point.
    # By "Pidgeon Hole principle" there'll be at least one point in each cell
    # that the geodesic traverse. Then for each point in the geodesic it gets
    # the closest south-west point in the grid and this identifies the cells

    # Once those cells are identified you can determine all intersections
    # Lets say the geodesic pass through the cells (x, y) and (x,y+1), then its
    # clear the geodesic intersects a 'horizontal' and the left and right
    # coordinates to that intersection are (x, y+1) and (x+1, y+1) (make a draw
    # if it helps)

    if start_end_rel_pos[1] == 'W':
        left_limit = np.max(np.array(jj1)) - 1
        right_limit = np.min(np.array(jj2)) + 1
    elif start_end_rel_pos[1] == 'E':
        left_limit = np.max(np.array(jj2)) - 1
        right_limit = np.min(np.array(jj1)) + 1

    if start_end_rel_pos[0] == 'S':
        bottom_limit = np.max(np.array(ii1)) - 1
        top_limit = np.min(np.array(ii2)) + 1
    elif start_end_rel_pos[0] == 'N':
        bottom_limit = np.max(np.array(ii2)) - 1
        top_limit = np.min(np.array(ii1)) + 1


    # 3. Obtener geodésica y sus respectivos intersecciones a la grilla en
    # índices y lat, lon

    i_points_in_between = int(((right_limit - left_limit)**2 + \
                              (top_limit - bottom_limit)**2)**0.5) + 1 - 2
    _geodesic = geodesic(point1, point2, i_points_in_between, sort=True)
    # again
    ii1, jj1 = get_neighbours_indexes(lats_arr, lons_arr, _geodesic[0])
    ii2, jj2 = get_neighbours_indexes(lats_arr, lons_arr, _geodesic[-1])

    # 3.1 Get cells touched by the geodesic i.e. all south-west indexes
    all_neighbours = map(lambda p: get_neighbours_indexes(lats_arr, lons_arr, p),
                         _geodesic)
    cells = [(nbs[0][0], nbs[1][0]) for nbs in all_neighbours]

    geodlats = []
    geodlons = []
    for y, x in _geodesic:
        geodlats.append(y)
        geodlons.append(x)


    # get lats or lons that intersect horizontals or verticals respectively
    # this DEPENDS HEAVELY on the sortint convention used when the geodesic is
    # generated

    # 'final path' stores final geodesic when all points are placed at
    # intersections except first and final
    final_path = [_geodesic[0]]
    neighbours_idxs = []  # neighbours of intermediate points
    if 0 <= inclination_angle <= 45: # horizontal
        interpolation_args = []
        for i in range(len(cells) - 1):
            if (cells[i] == cells[i+1]) or (cells[i][1] < cells[i+1][1]):
                interpolation_args.append(lons_arr[cells[i+1]])
                # points down and above of intersection
                neighbours_idxs.append(
                    (cells[i+1],
                    (cells[i+1][0] + 1, cells[i+1][1]))
                )
        # y = a*x + b
        interp_geodesic = interpolate.interp1d(geodlons, geodlats, kind='linear')
        for _lat, _lon in zip(interp_geodesic(interpolation_args), interpolation_args):
            final_path.append((_lat, _lon))

    elif 45 < inclination_angle <= 90: # vertical
        interpolation_args = []
        for i in range(len(cells) - 1):
            if (cells[i] == cells[i+1]) or (cells[i][0] < cells[i+1][0]):
                interpolation_args.append(lats_arr[cells[i+1]])
                # points left and right of intersection
                neighbours_idxs.append(
                    (cells[i+1],
                    (cells[i+1][0], cells[i+1][1] + 1))
                )
        # x = c*y + d
        interp_geodesic = interpolate.interp1d(geodlats, geodlons, kind='linear')
        for _lat, _lon in interpolation_args, interp_geodesic(interpolation_args):
            final_path.append((_lat, _lon))

    final_path.append(_geodesic[-1])

    # 4. Get interpolation of variables (contaminants, wind, temperature, etc)
    # with a loop
    # 4.1 Allocate memory
    planes = {}
    timesteps, height_levels = model_ds.height.shape[:2]
    if include_coordinates:
        variables += ['height']
        planes['lat'] = []
        planes['lon'] = []
        for lat, lon in final_path:
            planes['lat'].append(lat)
            planes['lon'].append(lon)
    for variable in variables:
        planes[variable] = np.full((timesteps, height_levels, len(final_path)),
                                   fill_value=np.nan,
                                   dtype=model_ds[variable].dtype)

    # 4.2 points outside the grid must be ignored
    # is_inside_array = geodesic_path_inside(lats_arr, lons_arr, final_path[1:-1])
    # NOTE: this step is no longer necessary since if a value is outside the
    # grid it will simply return NaN


    # 4.3 interpolate
    # initial and final points get 4-point interpolation. nan is are outside grid
    for t, z, variable in itertools.product(range(timesteps), \
                                            range(height_levels), variables):
        if is_inside(lats_arr, lons_arr, final_path[0]):
            planes[variable][t, z, 0] = \
                (interpolate.interp2d(x=lats_arr[ii1, jj1],
                                      y=lons_arr[ii1, jj1],
                                      z=model_ds[variable].values[t, z, ii1, jj1],
                                      kind='linear'))(*final_path[0])
        else:
            planes[variable][t, z, 0] = np.nan
        if is_inside(lats_arr, lons_arr, final_path[-1]):
            planes[variable][t, z, -1] = \
                (interpolate.interp2d(x=lats_arr[ii2, jj2],
                                      y=lons_arr[ii2, jj2],
                                      z=model_ds[variable].values[t, z, ii2, jj2],
                                      kind='linear'))(*final_path[-1])
        else:
            planes[variable][t, z, -1] = np.nan

    # remember, this array only contains the points in between
    for i, nbs in enumerate(neighbours_idxs):

        if 0 <= inclination_angle <= 45:
            x_source = lats_arr
            inperpolation_argument = final_path[i+1][0]
        elif 45 < inclination_angle <= 90:
            x_source = lons_arr
            inperpolation_argument = final_path[i+1][1]

        nb_ii = [nbs[0][0], nbs[1][0]]
        nb_jj = [nbs[0][1], nbs[1][1]]

        for t,z,variable in itertools.product(range(timesteps), \
                                              range(height_levels), variables):
            planes[variable][t, z, i+1] = \
                (interpolate.interp1d(x=x_source[nb_ii, nb_jj],
                                      y=model_ds[variable].values[t, z, nb_ii, nb_jj],
                                      kind='linear'))(inperpolation_argument)

    return planes

def four_points_interpolation_cross_section(model_ds, point1, point2, variables,
                                            include_coordinates=True,
                                            resolution='model'):
    """
    Gets cross section using equidistant points

    Args:
        Same as cross_section except:
        model: String or int. Default: 'model'. Custom distance between each
        point in the cross section. Note that real distance will vary since
        points are always placed in a equidistant way
    """

    # Gets resolution from file name's model
    if resolution == 'model':
        model_name = model_ds.attrs['model_name']
        res_location = model_name.index('KM') - 2
        model_spatial_resolution = int(model_name[res_location: res_location \
                                                  + 2]) * 1000  # in meters
    else:
        model_spatial_resolution = int(resolution)

    npts = int(distance(point1, point2) / model_spatial_resolution)
    geodesic_path = geodesic(point1, point2, npts, sort=True)

    lats_arr = model_ds.lat.values
    lons_arr = model_ds.lon.values

    planes = {}

    if include_coordinates:
        variables = list(variables + ['height'])  # avoids side effects
        planes['lat'] = []
        planes['lon'] = []
        for lat, lon in geodesic_path:
            planes['lat'].append(lat)
            planes['lon'].append(lon)

    timesteps, height_levels, _, _ = model_ds.height.shape

    # reserve memory for planes
    for variable in variables:
        planes[variable] = np.full((timesteps, height_levels,
                                    len(geodesic_path)), fill_value=np.nan,
                                   dtype=model_ds[variable].dtype)

    # get is_inside_array
    is_inside_array = geodesic_path_inside(lats_arr, lons_arr, geodesic_path)

    # loop on geodesic_path
    for i, point in enumerate(geodesic_path):
        if is_inside_array[i] == False:
            continue

        ii, jj = get_neighbours_indexes(lats_arr, lons_arr, point)

        # The coordinate concides with a point in the grid. Check how
        # 'get_neighbours_indexes' works.
        if len(ii) + len(jj) == 4:
            for variable in variables:
                planes[variable][:, :, i] = model_ds[variable] \
                                            .values[:, :, ii[0], jj[0]]
            continue

        neighbours_lat = lats_arr[ii, jj]
        neighbours_lon = lons_arr[ii, jj]

        # The coordinate is on a 'row' or a 'column', so only has two values to
        # make the interpolation
        if 2 in (len(ii), len(jj)):
            for t, z, variable in itertools.product(range(timesteps), \
                                                    range(height_levels), variables):
                x = neighbours_lat if len(ii) == 2 else neighbours_lon
                interp_func = interpolate.interp1d(x=x,
                                                   y=model_ds[variable].values[t, z, ii, jj],
                                                   kind='linear')
                planes[variable][t, z, i] = interp_func(point[0], point[1])

        # Usual case when 4 points exist to perform the interpolation
        else:
            for t, z, variable in itertools.product(range(timesteps), \
                                                    range(height_levels), variables):
                interp_func = interpolate.interp2d(x=neighbours_lat,
                                                   y=neighbours_lon,
                                                   z=model_ds[variable].values[t, z, ii, jj],
                                                   kind='linear')
                planes[variable][t, z, i] = interp_func(point[0], point[1])

    return planes

def get_neighbours_indexes(lats_arr, lons_arr, point):
    """
    Returns the indexes of at most four neighbours of 'point' in the grid.

    Examples:

    > print(lats_arr)
    [[-37.38575745 -37.38787079 -37.3889389  -37.3889389 ]
     [-37.11677933 -37.11888123 -37.11994171 -37.11994171]
     [-36.84789276 -36.84999847 -36.85104752 -36.85104752]
     [-36.57912064 -36.58119965 -36.5822525  -36.5822525 ]]

    > print(lons_arr)
    [[-71.52651978 -71.18792725 -70.8493042  -70.51068115]
     [-71.52319336 -71.18591309 -70.84863281 -70.51135254]
     [-71.51989746 -71.18395996 -70.84799194 -70.51202393]
     [-71.51663208 -71.18200684 -70.84732056 -70.51266479]]

    # A point inside
    > ii, jj = get_neighbours_indexes(lats_arr, lons_arr, (-37, -71))
    > print(ii, jj)
    [1, 1, 2, 2] [1, 2, 1, 2]

    > print(lats_arr[ii,jj])
    [-37.11888123, -37.11994171, -36.84999847, -36.85104752]

    > print(lons_arr[ii, jj])
    [-71.18591309, -70.84863281, -71.18395996, -70.84799194]

    """
    closest_idx = closest_idx_in_grid(lats_arr, lons_arr, point)
    relative_SN, relative_WE = relative_position(point, (lats_arr[closest_idx],
                                                         lons_arr[closest_idx]))
    nrows, ncols = lats_arr.shape
    if relative_SN == 'S' and closest_idx[0] != 0:
        ii = [closest_idx[0] - 1, closest_idx[0] - 1, closest_idx[0], closest_idx[0]]
    elif relative_SN == 'N' and closest_idx[0] != (nrows-1):
        ii = [closest_idx[0], closest_idx[0], closest_idx[0] + 1, closest_idx[0] + 1]
    else:
        ii = [closest_idx[0], closest_idx[0]]

    if relative_WE == 'W' and closest_idx[1] != 0:
        jj = [closest_idx[1] - 1, closest_idx[1], closest_idx[1] - 1, closest_idx[1]]
    elif relative_WE == 'E' and closest_idx[1] != (ncols-1):
        jj = [closest_idx[1], closest_idx[1] + 1, closest_idx[1], closest_idx[1] + 1]
    else:
        jj = [closest_idx[1], closest_idx[1]]

    if len(ii) == 2:
        jj = jj[:2]
    elif len(jj) == 2:
        ii = [ii[0], ii[2]]

    return ii, jj

def is_inside(lats_arr, lons_arr, point):
    """
    Tells whether 'point' is inside the grid. This doesn't means if it's a point in the grid.

    Args:
        lats_arr: Two-dimensional latitudes array
        lons_arr: Two-dimensional longitudes array
        point: Seq(float, float) lat lon point

    Returns:
        True or False.
    """
    idxs = closest_idx_in_grid(lats_arr, lons_arr, point)
    relative_SN, relative_WE = relative_position(point, (lats_arr[idxs],
                                                         lons_arr[idxs]))
    if relative_SN + relative_WE == '':
        return True
    nrows, ncolumns = lats_arr.shape
    return (not ((relative_SN == 'S' and idxs[0] == 0) or \
                 (relative_SN == 'N' and idxs[0] == (nrows - 1)))) \
           and \
           (not ((relative_WE == 'W' and idxs[1] == 0) or \
                 (relative_WE == 'E' and idxs[1] == (ncolumns - 1))))

def geodesic_path_inside(lats_arr, lons_arr, geodesic_path):
    """
    Checks if all points in the geodesic lay inside the grid.

    Equivalent to list(map(lambda p: is_inside(lats_arr, lons_arr, p),
    geodesic_path))
    but should be faster.

    It checks the path from the edges to the center, if a point is inside the
    grid, all subsequent points will be.
    """
    l = len(geodesic_path)
    ret = np.full(l, True)
    # left to right loop
    for i in range(l):
        if is_inside(lats_arr, lons_arr, geodesic_path[i]):
            break
        ret[i] = False
    # right to left loop
    for i in range(1, l+1):
        j = -i
        if is_inside(lats_arr, lons_arr, geodesic_path[j]):
            break
        ret[j] = False

    return ret

def concat_cross_sections(t_Cross_Sections):
    """
    Concatenates cross sections with same name and model.
    A model is usually splitted into several files which leads to also have
    cross sections splitted. After all cross sections are computed this function
    should be called.
    """
    # Group all xss by name and model
    out = []
    grouped_xss = {}
    for index, xss in enumerate(t_Cross_Sections):
        key = (xss['c_name'], xss['c_model'])
        if not key in grouped_xss:
            grouped_xss[key] = []
        # xss index and first date is stored to be sorted later
        # it is assumed that time arrays don't overlap since model files don't
        grouped_xss[key].append((index, xss['time'][0]))

    # I will use 'structured arrays' to simplify the sorting code
    dtype = np.dtype([('index', 'int32'), ('timestamp', 'datetime64[us]')])
    for keypair in grouped_xss:
        arr = np.array(grouped_xss[keypair], dtype=dtype)
        arr = arr[np.argsort(arr['timestamp'])]
        # Perform concatenation
        # shared data
        xss = {
            'c_name': keypair[0],
            'c_model': keypair[1],
            'f_bbox': t_Cross_Sections[arr['index'][0]]['f_bbox'],
            'lat': t_Cross_Sections[arr['index'][0]]['lat'],
            'lon': t_Cross_Sections[arr['index'][0]]['lon'],
        }
        # all other keys must be concatenated
        for xss_key in t_Cross_Sections[arr['index'][0]]:
            if xss_key in ('c_name', 'c_model', 'f_bbox', 'lat', 'lon'):
                # These keys are shared
                continue
            # Concatenates keys like time, height, PM25, PM10, etc
            xss[xss_key] = np.concatenate([
                t_Cross_Sections[index][xss_key] for index in arr['index']
            ])

        out.append(xss)

    return out

def compute_multicoordinate_cross_section(f_model, variables, xss_config, xss_name):
    """
    Computes cross sections defined with two or more coordinates.

    The algorithm computes normal cross sections by taking pairs of coordinates
    and concatenates them. This method has a drawback, data for all coordinates
    except the edges is computed twice.
    """

    # init dictionary with some keys
    xss = dict(c_name=xss_name,
               c_model=f_model.attrs['model_name'],  # TODO: use model alias
               f_bbox=get_model_bbox(f_model))  # Stores model bounding box for future plots

    partial_xss = []
    # iterate along coordinates
    for coord_idx in range(len(xss_config['path']) - 1):
        partial_xss.append(
            cross_section(model_ds=f_model,
                        point1=xss_config['path'][coord_idx],
                        point2=xss_config['path'][coord_idx + 1],
                        variables=variables)
            )
    # concatenate each partial cross section
    n_xss = len(partial_xss)
    if n_xss == 0:
        raise Exception('Unexpected error there are no cross sections. If we are here there should be')
    elif n_xss == 1:
        xss.update(partial_xss[0])
    else:  # n_xss > 1
        # Due to the cross section algorithm, each cross section can have
        # different orientation, i.e. from west to east or south to north. To
        # perform the desired concatenation we must match the edges coordinates:
        # ex:
        # Lets say we have an xss defined by the path (c1, c2, c3) and after
        # calculating cross sections taking pairs, each partial cross section
        # returns the coordinates as follows:
        #   partial_xss1_coords = [c1, ..., c2]
        #   partial_xss2_coords = [c3, ..., c2]
        # where the ellipsis (...) represents the generated coordinates.
        # In order to concatenate, the full coordinates array should be:
        # [c1, ..., c2, ..., c3] i.e., the user given order.

        # for each partial cross section we identify if the order is reversed or isn't
        variables = list(variables + ['height'])
        for i in range(len(xss_config['path']) - 1):
            # user defined coordinates
            c_i1 = np.array(xss_config['path'][i])
            # c_i2 = xss_config['path'][i+1] -> this isn't needed but is the next coord
            # `i` index also represents the xss
            pxss = partial_xss[i]
            edge1 = np.array((pxss['lat'][0], pxss['lon'][0]))
            edge2 = np.array((pxss['lat'][-1], pxss['lon'][-1]))
            # finding whether c_i1 == edge1 or c_i1 == edge2
            # the xss "is reversed" if the coordinate at the left is at the right
            # this is equivalent to calculate distance(c_i1, edge1) > distance(c_i1, edge2)
            isReversed = np.mean(np.abs(c_i1 - edge1)) > np.mean(np.abs(c_i1 - edge2))
            if isReversed:
                # reverse it, obviously
                # reverse lat, lon
                pxss['lat'] = pxss['lat'][::-1]
                pxss['lon'] = pxss['lon'][::-1]
                # reverse height and physical variables
                for var_key in variables:
                    # remember data dims/shape is (time, z, n_coords)
                    pxss[var_key] = pxss[var_key][:, :, ::-1]

        # When concatenating we slice omitting the last coordinate except for
        # the last xss
        # init this structure with the first xss
        concatenated_xss = dict(partial_xss[0])  # copy
        concatenated_xss['lat'] = concatenated_xss['lat'][:-1]
        concatenated_xss['lon'] = concatenated_xss['lon'][:-1]
        for var_key in variables:
            concatenated_xss[var_key] = concatenated_xss[var_key][:, :, :-1]
        # iterate until the penultimate xss to avoid 'if last'
        for i in range(1, len(xss_config['path']) - 2):
            pxss = partial_xss[i]
            # concatenate lat and lon
            concatenated_xss['lat'] = np.concatenate([concatenated_xss['lat'], pxss['lat'][:-1]])
            concatenated_xss['lon'] = np.concatenate([concatenated_xss['lon'], pxss['lon'][:-1]])
            # concatenate height and the model vars
            for var_key in variables:
                concatenated_xss[var_key] = np.concatenate([concatenated_xss[var_key], pxss[var_key][:, :, :-1]], axis=2)
        # concat last xss
        pxss = partial_xss[-1]
        # concatenate lat and lon
        concatenated_xss['lat'] = np.concatenate([concatenated_xss['lat'], pxss['lat']])
        concatenated_xss['lon'] = np.concatenate([concatenated_xss['lon'], pxss['lon']])
        # concatenate height and the model vars
        for var_key in variables:
            concatenated_xss[var_key] = np.concatenate([concatenated_xss[var_key], pxss[var_key]], axis=2)
        # TODO: avoid last code duplication
        xss.update(concatenated_xss)

    return xss



# TODO: create cache system with key-tuples (xss-name, model-name, xss-algorithm)
# to avoid recompute things

def main(f_model):
    """
    Computes cross sections for given model data

    This code depends on two external structures:
        - t_cross_sections_config: Configuration of cross sections from IncludeFile
        - config_by_models: computed by sort_xss_configs_by_model
    """
    try:
        config_by_models = sort_xss_configs_by_model(c_ModelNames, t_cross_sections_config)
    except ImportError:
        config_by_models = None
        print("WARNING: IncludeFile not found. I hope you're using this file in testing mode")
        
    model_xss = []

    for xss_index in config_by_models[f_model.attrs['model_name']]:
        xss_config = t_cross_sections_config[xss_index]
        xss_name = xss_config.get('name', f'CROSS_SECTION_{xss_index}')
        # Select model variables except lat, lon, time
        # This data can also be obtained from IncludeFile but this method is
        # better in order to achieve more decloupling degree
        variables = list(set(f_model.variables).difference(set(f_model.coords)))  # should be equals to model vars
        xss = compute_multicoordinate_cross_section(f_model, variables, xss_config, xss_name)
        model_xss.append(xss)

    return model_xss
