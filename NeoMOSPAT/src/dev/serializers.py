from typing import List, Any

from numpy import ndarray
from pandas import DataFrame, ExcelWriter

def ndarray_preserializer(array: ndarray) -> list:
    return array.tolist()

def dictionary_preserializer(_dict: dict) -> dict:
    for k, v in _dict.items():
        # only aware of not serializable objects
        if isinstance(v, ndarray):
            _dict[k] = ndarray_preserializer(v)
        elif isinstance(v, DataFrame):
            print("Preserializer for Pandas.DataFrame has not been implemented")
    return _dict

def preserialize(items: List[Any]) -> List[Any]:
    for i, item in enumerate(items):
        if isinstance(item, dict):
            items[i] = dictionary_preserializer(item)

    return items
