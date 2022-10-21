"""
Groups of coordinates in different systems
"""
from typing import NamedTuple
from math import floor


def degrees_to_dms(deg: float) -> str:
    """
    Converts float degree into degrees, minutes and seconds
    """
    degree = floor(deg)
    minute = floor((deg - degree) * 60)
    second = (deg - (degree + minute / 60)) * 3600
    return f"{degree}º {minute}' " + f'{second:.1f}"'


class Constants(NamedTuple):
    """
    Grouped constants
    """
    # Epsilon - angle of ecliptic inclanation
    epsilon: float
    cos_e: float
    sin_e: float
    tan_e: float

    # Phi - observer's geo latitude
    cos_p: float
    sin_p: float
    tan_p: float

    # Max allowable declination for Zodiac
    # degree to be ascendible at extreme latitudes
    dec_max: float


class Equatorial():
    """
    Vector in equatorial system
    """

    def __init__(self,
                 rasc: float,  # Right ascension
                 dec: float,  # Declination
                 ) -> None:
        self.rasc = rasc
        self.dec = dec

    def __str__(self) -> str:
        return (
            f'Equatorial(rasc={degrees_to_dms(self.rasc)}, '
            + f'dec={degrees_to_dms(self.dec)})'
        )


class Horizontal():
    """
    Vector in horizontal system
    """

    def __init__(self,
                 azm: float,  # Azimuth
                 alt: float,  # Altitude
                 ) -> None:
        self.azm = azm
        self.alt = alt

    def __str__(self) -> str:
        return (
            f'Horizontal(azm={degrees_to_dms(self.azm)}, '
            + f'alt={degrees_to_dms(self.alt)})'
        )


class Ecliptical():
    """
    Vector in ecliptical coordinate system
    """

    def __init__(self,
                 lon: float,  # Zodiac absolute degree
                 lat: float,  # Celestial Latitude
                 ) -> None:
        self.lon = lon
        self.lat = lat

    def __str__(self) -> str:
        return (
            f'Ecliptical(lon={degrees_to_dms(self.lon)}, '
            + f'lat={degrees_to_dms(self.lat)})'
        )


class Cartesian(NamedTuple):
    """
    Cartesian coordinates.
    """
    x: float
    y: float
    z: float
