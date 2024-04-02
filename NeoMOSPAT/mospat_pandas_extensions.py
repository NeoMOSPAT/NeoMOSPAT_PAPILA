# -*- coding: utf-8 -*-
# http://xarray.pydata.org/en/stable/internals.html#extending-xarray
# https://pandas.pydata.org/docs/development/extending.html

import numpy as np
import pandas as pd

@pd.api.extensions.register_dataframe_accessor("meta")
class MOSPATModelMetadata:
    """
    Stores variable metadata
    """

    def __init__(self, pandas_obj):
        self._df = pandas_obj
        self._model_name = ''
        self._source_folder = ''
        self._units = {}

    @property
    def units(self):
        return self._units

    def get_units(self, c_model_variable: str) -> str:
        return self._units.get(c_model_variable, '')

    @units.setter
    def units(self, var_units):
        """
        usage ex.:

            df = pd.DataFrame(data);
            df.meta.units = {'PM10': '[ton]/[m^3]'}
        """
        self._units.update(var_units)

    def add_units(self, c_model_variable : str, c_units : str) -> dict:
        self._units.update({c_model_variable : c_units})
        return self._units

@pd.api.extensions.register_dataframe_accessor("filters")
class MOSPATModelFilter:
    """
    Stores spatial filters indexes and provides methods to easily selecting
    filtered data.
    """

    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        self._filters : dict = {}

    def add_filter(self, name, indexes_groups):
        # indexes groups must be lists
        _indexes = None
        # it's only enough one group to be an ndarray to assume all the groups are
        if isinstance(indexes_groups[0], np.ndarray):
            _indexes = [ndarray_to_list(indexes_group) for indexes_group in indexes_groups]
        else:
            _indexes = indexes_groups
        self._filters[name] = _indexes

    @property
    def filters(self):
        return self._filters

    @property
    def names(self):
        return list(self._filters.keys())

    # TODO: change `group_index` to a label when labeled filters became stable
    def get_data(self, filter_name, group_index: int = 0, dates=None):
        """
        Gets the data filtered by the filter

        Args:
            filter_name
            group_index, Zero by default. Any of the groups you defined in your
            IncludeFile.
        """
        # assert self._obj.index.names
        indexes = None
        try:
            indexes : List[Tuple[int]] = self._filters[filter_name][group_index]
        except KeyError as e:
            raise KeyError(f"Unknown filter {filter_name}")

        return select_data_with_lat_lon_indexes(self._obj, indexes, dates)

    def get_data_with_indexes(self, indexes, dates=None):
        """
        Get data of the dataframe using the given lat, lon indexes. There's no
        need to store the filter previously.
        """
        # Transform to list
        indexes = ndarray_to_list(indexes)
        # Perform data selection
        return select_data_with_lat_lon_indexes(self._obj, indexes, dates)


def ndarray_to_list(arr):
    """
    Transform numpy array into list
    """
    return arr.tolist() if isinstance(arr, np.ndarray) else arr

def select_data_with_lat_lon_indexes(pandas_obj: pd.DataFrame, indexes: list, dates=[]) -> pd.DataFrame:
    """Performs data selection on a Dataframe using latitude and longitude indexes.

    Args:
        pandas_obj (pandas.DataFrame): Data source
        indexes (list): List of tuples containing latitude and longitude indexes. Missing values are allowed.
        dates (list, optional): Dates to select. Defaults to [].

    Returns:
        pandas.DataFrame: selected data
    """

    if not dates:
        ret = pandas_obj
    else:
        ret = pandas_obj.loc[dates]

    selected_dates = ret.index.get_level_values('time').unique()

    partial_selection = []

    for _date in selected_dates:
        try:
            partial_df = ret.loc[_date].loc[indexes]
        except KeyError:
            # Usually because of this:
            # https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-with-list-with-missing-labels-is-deprecated
            partial_df = ret.loc[_date].reindex(indexes)

        # recreate the index because date index has been lost
        new_index = pd.MultiIndex.from_arrays(
            [[_date] * len(partial_df),
            partial_df.index.get_level_values(0),
            partial_df.index.get_level_values(1)],
            names=('times', 'i_lat', 'i_lon'))
        partial_df.index = new_index

        partial_selection.append(partial_df)

    return pd.concat(partial_selection)
