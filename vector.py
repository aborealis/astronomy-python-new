"""
Creates vector and transfers it between
spherical coordinate systems
"""
from cmath import tan
from datetime import datetime, timedelta
from math import asin, atan
from numpy import cos, sin, pi, tan
from typing import NamedTuple, Union
import components.time as tm


def atan_btw_xy(x: float, y: float) -> float:
    """
    Returns angle (in degrees) between cartesian x and y
    """
    alpha = atan(y / x) / pi * 180
    if x < 0 and y > 0:
        alpha += 180
    elif x < 0 and y < 0:
        alpha += 180
    elif x > 0 and y < 0:
        alpha += 360
    return alpha


def true_distance(abs_degree1: float, abs_degree2: float) -> float:
    """
    Returns minimum angle between two degrees.
    """
    distance = abs_degree2 - abs_degree1
    if abs(distance) > 180:
        if abs_degree2 < abs_degree1:
            distance += 360
        else:
            distance -= 360
    return abs(distance)


class __Constants__(NamedTuple):
    """
    Constants shared amomg functions
    """
    # Epsilon - angle of ecliptic inclanation
    cos_e: float
    sin_e: float
    tan_e: float

    # Phi - observer's geo latitude
    cos_p: float
    sin_p: float
    tan_p: float


class __Equatorial__(NamedTuple):
    """
    Vector in equatorial system
    """
    ra: float  # Real Ascention
    dec: float  # Declination


class __Horizontal__(NamedTuple):
    """
    Vector in horizontal system
    """
    azm: float  # Azimuth
    alt: float  # Altitude


class __Ecliptical__(NamedTuple):
    """
    Vector in ecliptical coordinate system
    """
    lon: float  # Zodiac absolute degree
    lat: float  # Celestial Latitude


class __Cartesian__(NamedTuple):
    """
    Cartesian coordinates.
    """
    x: float
    y: float
    z: float


class __Coords__(NamedTuple):
    """
    Types of coordinate system
    """
    ecliptical: __Ecliptical__
    horizontal: __Horizontal__


class Vector:
    """
    Requires local time and geoposition
    of the observer
    """

    def __init__(self,
                 naive_datetime: datetime,
                 time_zone: float,
                 geo_lon: float,
                 geo_lat: float) -> None:
        localtime = tm.__localtime__(naive_datetime, time_zone)
        localtime_utc = tm.__localtime_utc__(localtime)
        epsilon = tm.__inclination_ecliptic__(localtime_utc)

        self.__constants__ = __Constants__(
            cos_e=cos(epsilon * pi / 180),
            sin_e=sin(epsilon * pi / 180),
            tan_e=tan(epsilon * pi / 180),
            cos_p=cos(geo_lat * pi / 180),
            sin_p=sin(geo_lat * pi / 180),
            tan_p=tan(geo_lat * pi / 180),
        )

        self.coords = None
        self.lst = tm.__lst__(localtime_utc, geo_lon)

    def set_equatorial(self, ra: float, dec: float) -> __Equatorial__:
        """
        Sets equatorial coordinates
        """
        self.coords = __Equatorial__(ra, dec)

    def set_ecliptical(self, lon: float, lat: float) -> __Ecliptical__:
        """
        Sets equatorial coordinates
        """
        self.coords = __Ecliptical__(lon, lat)

    def set_horizontal(self, azm: float, alt: float) -> __Horizontal__:
        """
        Sets horizontal coordinates
        """
        self.coords = __Horizontal__(azm, alt)

    def cartesian(self,
                  vector: Union[
                      __Ecliptical__,
                      __Equatorial__,
                      __Horizontal__,
                      None,
                  ] = None) -> __Cartesian__:
        """
        Transforms spherical coordinates
        into cartesian (x, y, z).

        Directions:
        - Equatorial(ra=0, dec=0) -> (1, 0, 0)
        - Ecliptical(lon=0, lat=0) -> (1, 0, 0)
        - Horizontal(azm=0, alt=0) -> (0, 1, 0)
        """
        coords = self.coords if vector is None else vector
        if isinstance(coords, __Equatorial__):
            lon = coords.ra
            lat = coords.dec
        elif isinstance(coords, __Ecliptical__):
            lon = coords.lon
            lat = coords.lat
        else:
            lon = 90 - coords.azm
            lat = coords.alt

        cos_lat = cos(lat / 180 * pi)
        sin_lat = sin(lat / 180 * pi)
        cos_lon = cos(lon / 180 * pi)
        sin_lon = sin(lon / 180 * pi)

        x = cos_lat * cos_lon
        y = cos_lat * sin_lon
        z = sin_lat

        return __Cartesian__(x, y, z)

    def cartesian_horizontal(self):
        """
        This function is created for numpy
        visualization only
        """
        cartesian = self.cartesian()
        coords = self.coords
        if isinstance(coords, __Equatorial__):
            return self.__eqt_to_hrz_cartesian__(cartesian)
        if isinstance(coords, __Ecliptical__):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(cartesian)
            return self.__eqt_to_hrz_cartesian__(xyz_eqt)
        return cartesian

    @staticmethod
    def __cartesian_to_spherical__(coords: __Cartesian__) -> dict:
        """
        Transforms cartesian (x, y, z) coordinates
        into spherical angles
        """
        if coords.x == 0:
            lon = 270
        else:
            lon = atan_btw_xy(coords.x, coords.y)
        lat = asin(coords.z) / pi * 180

        return dict(horz_angle=lon, vert_angle=lat)

    def __ecl_to_eqt_cartesian__(self, ecl: __Cartesian__) -> __Cartesian__:
        """
        Transfers cartesian (x, y, z) from
        ecliptical to equatorial system
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        x_eqt = ecl.x
        y_eqt = cos_e * ecl.y - sin_e * ecl.z
        z_eqt = sin_e * ecl.y + cos_e * ecl.z

        return __Cartesian__(x_eqt, y_eqt, z_eqt)

    def __eqt_to_ecl_cartesian__(self, eqt: __Cartesian__) -> __Cartesian__:
        """
        Transfers cartesian (x, y, z) from
        equatorial to ecliptical system
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        x_ect = eqt.x
        y_ect = cos_e * eqt.y + sin_e * eqt.z
        z_ect = - sin_e * eqt.y + cos_e * eqt.z

        return __Cartesian__(x_ect, y_ect, z_ect)

    def __eqt_to_hrz_cartesian__(self, eqt: __Cartesian__) -> __Cartesian__:
        """
        Transfers cartesian (x, y, z) from
        equatorial to horizontal system
        """
        # time since Aries raising (tsa, in radians)
        tsa = (self.lst * 15 + 90) / 180 * pi
        sin_t = sin(tsa)
        cos_t = cos(tsa)
        cos_p = self.__constants__.cos_p
        sin_p = self.__constants__.sin_p

        # rotated vector with (1, 0, 0) pointing to (ra=0, dec=0)
        x_eqt_rot = cos_t * eqt.x + sin_t * eqt.y
        y_eqt_rot = - sin_t * eqt.x + cos_t * eqt.y
        z_eqt_rot = eqt.z

        x_hrz = x_eqt_rot
        y_hrz = sin_p * y_eqt_rot + cos_p * z_eqt_rot
        z_hrz = - cos_p * y_eqt_rot + sin_p * z_eqt_rot

        return __Cartesian__(x_hrz, y_hrz, z_hrz)

    def __htz_to_eqt_cartesian__(self, hrz: __Cartesian__) -> __Cartesian__:
        """
        Transfers cartesian (x, y, z) from
        horizontal to equatorial system
        """
        # time since Aries raising (tsa, in radians)
        tsa = (self.lst * 15 + 90) / 180 * pi
        sin_t = sin(tsa)
        cos_t = cos(tsa)
        cos_p = self.__constants__.cos_p
        sin_p = self.__constants__.sin_p

        # rotated vector with (1, 0, 0) pointing to (ra=0, dec=0)
        x_eqt_rot = hrz.x
        y_eqt_rot = sin_p * hrz.y - cos_p * hrz.z
        z_eqt_rot = cos_p * hrz.y + sin_p * hrz.z

        x_eqt = cos_t * x_eqt_rot - sin_t * y_eqt_rot
        y_eqt = sin_t * x_eqt_rot + cos_t * y_eqt_rot
        z_eqt = z_eqt_rot

        return __Cartesian__(x_eqt, y_eqt, z_eqt)

    def equatorial(self) -> __Equatorial__:
        """
        Transfers current vector to equatorial
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, __Ecliptical__):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz)
        elif isinstance(self.coords, __Horizontal__):
            xyz_eqt = self.__htz_to_eqt_cartesian__(xyz)
        else:
            xyz_eqt = xyz

        spherical = self.__cartesian_to_spherical__(xyz_eqt)
        return __Equatorial__(spherical['horz_angle'],
                              spherical['vert_angle'])

    def horizontal(self) -> __Horizontal__:
        """
        Transfers current vector to horizontal
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, __Ecliptical__):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz)
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
        elif isinstance(self.coords, __Equatorial__):
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz)
        else:
            xyz_hrz = xyz

        spherical = self.__cartesian_to_spherical__(xyz_hrz)
        return __Horizontal__(90 - spherical['horz_angle'],
                              spherical['vert_angle'])

    def ecliptical(self) -> __Ecliptical__:
        """
        Transfers current vector to ecliptical
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, __Equatorial__):
            xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz)
        elif isinstance(self.coords, __Horizontal__):
            xyz_eqt = self.__htz_to_eqt_cartesian__(xyz)
            xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz_eqt)
        else:
            xyz_ecl = xyz

        spherical = self.__cartesian_to_spherical__(xyz_ecl)
        return __Ecliptical__(spherical['horz_angle'],
                              spherical['vert_angle'])

    def oblique_asc(self):
        """
        Returns oblique ascention of <self> object
        """
        tan_p = self.__constants__.tan_p
        dec = self.equatorial().dec
        ra = self.equatorial().ra
        tan_d = tan(dec * pi / 180)
        ascention_diff = asin(tan_p * tan_d) * 180 / pi
        return (ra - ascention_diff) % 360

    def ramc(self):
        """
        Real ascention of medium coeli for
        given time in <self> object
        """
        return self.lst * 15

    def asc(self):
        """
        Returns ascending Zodiac degree
        """
        # time since Aries raising (tsa, in radians)
        tsa = (self.lst * 15 + 90) % 360
        cos_t = cos(tsa * pi / 180)
        sin_t = sin(tsa * pi / 180)
        tan_p = self.__constants__.tan_p
        tan_e = self.__constants__.tan_e
        y = sin_t * tan_p * tan_e
        x = 1 - cos_t * tan_p * tan_e
        diff = atan_btw_xy(x, y)

        ra_asc = (tsa + diff) % 360
        y = sin(diff * pi / 180)
        x = tan_p
        dec_asc = atan_btw_xy(x, y)

        cartesian = self.cartesian(__Equatorial__(ra_asc, dec_asc))
        xyz_ecl = self.__eqt_to_ecl_cartesian__(cartesian)

        return self.__cartesian_to_spherical__(xyz_ecl)['horz_angle']

    def umd(self):
        """
        Returns upper meridian distance, UMD
        of <self> object
        """
        return true_distance(
            self.equatorial().ra,
            self.ramc()
        )

    def dsa(self):
        """
        Returns diurnal semiarc
        """
        asc_diff = self.equatorial().ra - self.oblique_asc()
        return (90 + asc_diff) % 360


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':

    # Initiate a vector object
    vector = Vector(
        datetime(2022, 10, 10, 21, 58),
        time_zone=4,
        geo_lon=44 + 46/60,
        geo_lat=41 + 43/50,
    )

    # Get celestial information
    # Local Sidereal Time and real ascention of MC
    print('LST:', timedelta(hours=vector.lst))
    print('RAMC:', vector.ramc())

    # Ascendant Zodiac degree
    asc = vector.asc()
    print('ASC:', asc)

    # Set coordinates in a chosen system
    vector.set_horizontal(azm=90.0, alt=0.0)

    # Get the same vector in different systems
    print('True East:', vector.cartesian())
    print('True East:', vector.horizontal())
    print('True East:', vector.equatorial())
    print('True East:', vector.ecliptical())
    print('True East oblique ascention:', vector.oblique_asc())

    # Observe cartasian coordinates in any system
    print(vector.cartesian(vector.ecliptical()))

    # Let's check ascendant
    vector.set_ecliptical(lon=vector.asc(), lat=0.0)
    print('ASC:', vector.horizontal())
    print('ASC:', vector.equatorial())

    # The same as real ascention of True East above
    print('Oblique Ascention of ASC:', vector.oblique_asc())

    # Uppper Meridian Distance (UMD) and
    # Diurnal Semiarc (DSA) of a given vector
    print('UMD:', vector.umd())
    print('DSA:', vector.dsa())
