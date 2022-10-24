"""
Groups of coordinates in different systems
"""
from typing import NamedTuple, Callable
from math import copysign


def degrees_to_dms(deg: float) -> str:
    """
    Converts float degree into degrees, minutes and seconds
    """
    degree = int(abs(deg * 100) // 100 * copysign(1, deg))
    minute = int(abs((deg - degree) * 60) * 10 // 10)
    second = round((abs(deg - degree) * 60 - minute) * 60)
    return f"{degree}º {minute:02d}'" + f'{second:02d}"'


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

    # Sid. time since Aries ascended, tsa
    sin_t: float
    cos_t: float


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
            f'Equatorial(rasc={degrees_to_dms(self.rasc)}|{self.rasc:.4f}º, '
            + f'dec={degrees_to_dms(self.dec)}|{self.dec:.4f}º)'
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
            f'Horizontal(azm={degrees_to_dms(self.azm)}|{self.azm:.4f}º, '
            + f'alt={degrees_to_dms(self.alt)}|{self.alt:.4f}º)'
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
            f'Ecliptical(lon={degrees_to_dms(self.lon)}|{self.lon:.4f}º, '
            + f'lat={degrees_to_dms(self.lat)}|{self.lat:.4f}º)'
        )


class Cartesian:
    """
    Cartesian coordinates.
    """

    # pylint: disable=C0103
    def __init__(self,
                 x: float,
                 y: float,
                 z: float,
                 ) -> None:
        self.x, self.y, self.z = x, y, z

    def __str__(self) -> str:
        return (
            f'Cartesian(x={self.x:.2f}, y={self.y:.2f}, z={self.z:.2f})'
        )

    def aslist(self):
        """
        Returns cartesian coordinates as a list
        """
        return [self.x, self.y, self.z]


class Vector(NamedTuple):
    """
    Once vector is set on the sphere,
    it becomes a member of this class
    """
    ecliptical: Callable[[], Ecliptical]
    equatorial: Callable[[], Equatorial]
    horizontal: Callable[[], Horizontal]
    ecliptical_xyz: Callable[[], Cartesian]
    equatorial_xyz: Callable[[], Cartesian]
    horizontal_xyz: Callable[[], Cartesian]
    obl_asc: Callable[[], float]
    asc_diff: Callable[[], float]
    umd: Callable[[], float]
    dsa: Callable[[], float]
