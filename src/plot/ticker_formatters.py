from matplotlib import ticker

class DMS:
    """
    Decimal, minutes, seconds

    See examples in examples/API/DMS.ipynb
    """
    def __init__(self, coord: float):
        """
        Defines self.integer, self.minutes, self.seconds and self._sign. The
        latest is useful when implementing operations like __add__ to preserve
        result sign.\n
        self.module is used in the children classes LatDMS and LonDMS to define
        effective range.
        """
        self.integer, self.minutes, self.seconds = self.transform_coord(coord)
        if self.integer >= 0:
            self._sign = 1
        else:
            self._sign = -1
        self._module = 1

    def _test_validity(self):
        if self._module != 1:
            if (abs(self.integer) > self._module) \
            or (abs(self.integer) == self._module and (self.minutes > 0 or self.seconds > 0)):
                raise Exception(f"Invalid value for coordinate {self}")

    @staticmethod
    def transform_coord(coord: float) -> tuple:
        integer = int(coord)
        decimals = coord - integer
        _minutes = int(decimals * 60)
        minutes = abs(_minutes)
        seconds = abs((decimals - _minutes/60) * 3600)
        return integer, minutes, seconds

    @staticmethod
    def _tostr(integer: int, minutes: int, seconds:int) -> str:
        # WARNING: DECIMALS OF SECONDS ARE IGNORED WHEN COVERTING TO STRING
        return f"{abs(integer)}°{minutes:02}'{int(seconds):02}\""

    def __str__(self) -> str:
        return self._tostr(self.integer, self.minutes, self.seconds) + '_'

    @classmethod
    def get_full_formatter(cls):
        @ticker.FuncFormatter
        def full_formatter(x, pos):
            return str(cls(x))
        return full_formatter

    @classmethod
    def get_minutes_formatter(cls):
        @ticker.FuncFormatter
        def minutes_formatter(x, pos):
            s = str(cls(x))
            # remove seconds
            secs_idx = s.index('\'') + 1
            # the whole string plus the letter
            return s[:secs_idx] + s[-1]
        return minutes_formatter

    @classmethod
    def get_degrees_formatter(cls):
        @ticker.FuncFormatter
        def degrees_formatter(x, pos):
            s = str(cls(x))
            # remove minutes and seconds
            degrees_idx = s.index('°') + 1
            # the whole string plus the letter
            return s[:degrees_idx] + s[-1]
        return degrees_formatter

    def to_string(self, resolution: str = 'full') -> str:
        """Transform to string with desired resolution

        Args:
            resolution (`str`, optional): Possible values `'full'`, `'minutes'`,\n
            `'degrees'`. Defaults to `'full'`.

        Returns:
            str: Formatted coordinate

        Example:
        ```
            > lat_dms = LatDMS(-33.458066)
            > print(lat_dms.to_string(resolution='full'))
            33°27'29"S

            > print(lat_dms.to_string(resolution='minutes'))
            33°27'S

            > print(lat_dms.to_string(resolution='degrees'))
            33°S
        ```
        """
        if resolution == 'minutes':
            s = str(self)
            # remove seconds
            secs_idx = s.index('\'') + 1
            # the whole string plus the letter
            return s[:secs_idx] + s[-1]
        elif resolution == 'degrees':
            s = str(self)
            # remove minutes and seconds
            degrees_idx = s.index('°') + 1
            # the whole string plus the letter
            return s[:degrees_idx] + s[-1]

        # resolution == 'full', to avoid else
        return str(self)

class LatDMS(DMS):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._module = 90
        self._test_validity()

    def __str__(self):
        direction_char = 'N' if self._sign > 0 else 'S'
        return f"{DMS._tostr(self.integer, self.minutes, self.seconds)}{direction_char}"

    def __repr__(self):
        return f"LatDMS({str(self)})"

class LonDMS(DMS):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._module = 180
        self._test_validity()

    def __str__(self):
        direction_char = 'E' if self._sign > 0 else 'W'
        return f"{DMS._tostr(self.integer, self.minutes, self.seconds)}{direction_char}"

    def __repr__(self):
        return f"LonDMS({str(self)})"


def test():
    lats = [-33.458066]
    lons = [-70.663402]

    # Print coordinates in decimal, minutes and seconds format
    print(LatDMS(lats[0]).to_string(resolution='full'))
    print(LatDMS(lats[0]).to_string(resolution='minutes'))
    print(LatDMS(lats[0]).to_string(resolution='degrees'))

    # Usage as matplotlib ticks formatters
    # https://matplotlib.org/stable/gallery/ticks_and_spines/tick-formatters.html

    # first you get the functions that will be use as formatters
    lat_degrees_formatter = LatDMS.get_degrees_formatter()
    lon_minutes_formatter = LonDMS.get_minutes_formatter()

    # this is how these functions format the ticks
    print(lat_degrees_formatter(lats[0], 0))
    print(lon_minutes_formatter(lons[0], 0), '\n')

    # now in a plot
    import matplotlib.pyplot as plt
    import numpy as np

    # grid points with decimal degrees latitudes and longitudes
    x = np.linspace(-33.21, -28.54, 15)
    y = np.linspace(-5.04, 15.66, 15)

    fig, ax = plt.subplots(1,1)

    ax.scatter(x, y)
    ax.xaxis.set_major_formatter(lat_degrees_formatter)
    ax.yaxis.set_major_formatter(lon_minutes_formatter)

    plt.show()



if __name__=='__main__':
    test()
    # import these classes using
    # from src.plot.ticker_formatters import LatDMS, LonDMS
    # or
    # import src.plot.ticker_formatters as tf
    # tf.LatDMS, tf.LonDMS