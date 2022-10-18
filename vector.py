"""
Creates vector and transfers it between
spherical coordinate systems
"""
from typing import NamedTuple, Optional, Union
from datetime import datetime, timedelta
from math import asin, atan, floor
from numpy import cos, sin, pi, tan
from components import time as tm


def atan_btw_xy(_x: float, _y: float) -> float:
    """
    Returns angle (in degrees) between cartesian x and y
    """
    alpha = atan(_y / _x) / pi * 180
    if _x < 0 and _y > 0:
        alpha += 180
    elif _x < 0 and _y < 0:
        alpha += 180
    elif _x > 0 and _y < 0:
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


class _Constants(NamedTuple):
    """
    Constants shared amomg functions
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


class _Equatorial(NamedTuple):
    """
    Vector in equatorial system
    """
    ra: float  # Real Ascention
    dec: float  # Declination


class _Horizontal(NamedTuple):
    """
    Vector in horizontal system
    """
    azm: float  # Azimuth
    alt: float  # Altitude


class _Ecliptical(NamedTuple):
    """
    Vector in ecliptical coordinate system
    """
    lon: float  # Zodiac absolute degree
    lat: float  # Celestial Latitude


class _Cartesian(NamedTuple):
    """
    Cartesian coordinates.
    """
    x: float
    y: float
    z: float


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

        self.__constants__ = _Constants(
            cos_e=cos(epsilon * pi / 180),
            sin_e=sin(epsilon * pi / 180),
            tan_e=tan(epsilon * pi / 180),
            cos_p=cos(geo_lat * pi / 180),
            sin_p=sin(geo_lat * pi / 180),
            tan_p=tan(geo_lat * pi / 180),
            epsilon=epsilon,
            dec_max=90 - geo_lat,
        )

        self.coords = None
        self.lst = tm.__lst__(localtime_utc, geo_lon)

    def set_equatorial(self, real_asc: float, dec: float) -> _Equatorial:
        """
        Sets equatorial coordinates
        """
        self.coords = _Equatorial(real_asc, dec)

    def set_ecliptical(self, lon: float, lat: float) -> _Ecliptical:
        """
        Sets equatorial coordinates
        """
        self.coords = _Ecliptical(lon, lat)

    def set_horizontal(self, azm: float, alt: float) -> _Horizontal:
        """
        Sets horizontal coordinates
        """
        self.coords = _Horizontal(azm, alt)

    def cartesian(self,
                  vector: Union[
                      _Ecliptical,
                      _Equatorial,
                      _Horizontal,
                      None,
                  ] = None) -> _Cartesian:
        """
        Transforms spherical coordinates
        into cartesian (x, y, z).

        Directions:
        - Equatorial(ra=0, dec=0) -> (1, 0, 0)
        - Ecliptical(lon=0, lat=0) -> (1, 0, 0)
        - Horizontal(azm=0, alt=0) -> (0, 1, 0)
        """
        coords = self.coords if vector is None else vector
        if isinstance(coords, _Equatorial):
            lon = coords.ra
            lat = coords.dec
        elif isinstance(coords, _Ecliptical):
            lon = coords.lon
            lat = coords.lat
        else:
            lon = 90 - coords.azm
            lat = coords.alt

        cos_lat = cos(lat / 180 * pi)
        sin_lat = sin(lat / 180 * pi)
        cos_lon = cos(lon / 180 * pi)
        sin_lon = sin(lon / 180 * pi)

        return _Cartesian(
            x=cos_lat * cos_lon,
            y=cos_lat * sin_lon,
            z=sin_lat,
        )

    def cartesian_horizontal(self):
        """
        This function is created for numpy
        visualization only
        """
        cartesian = self.cartesian()
        coords = self.coords
        if isinstance(coords, _Equatorial):
            return self.__eqt_to_hrz_cartesian__(cartesian)
        if isinstance(coords, _Ecliptical):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(cartesian)
            return self.__eqt_to_hrz_cartesian__(xyz_eqt)
        return cartesian

    @staticmethod
    def __cartesian_to_spherical__(coords: _Cartesian) -> dict:
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

    def __ecl_to_eqt_cartesian__(self, ecl: _Cartesian) -> _Cartesian:
        """
        Transfers cartesian (x, y, z) from
        ecliptical to equatorial system
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        x_eqt = ecl.x
        y_eqt = cos_e * ecl.y - sin_e * ecl.z
        z_eqt = sin_e * ecl.y + cos_e * ecl.z

        return _Cartesian(x_eqt, y_eqt, z_eqt)

    def __eqt_to_ecl_cartesian__(self, eqt: _Cartesian) -> _Cartesian:
        """
        Transfers cartesian (x, y, z) from
        equatorial to ecliptical system
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        x_ecl = eqt.x
        y_ecl = cos_e * eqt.y + sin_e * eqt.z
        z_ecl = - sin_e * eqt.y + cos_e * eqt.z

        return _Cartesian(x_ecl, y_ecl, z_ecl)

    def __eqt_to_hrz_cartesian__(self, eqt: _Cartesian) -> _Cartesian:
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

        return _Cartesian(x_hrz, y_hrz, z_hrz)

    def __htz_to_eqt_cartesian__(self, hrz: _Cartesian) -> _Cartesian:
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

        return _Cartesian(x_eqt, y_eqt, z_eqt)

    def equatorial(self) -> _Equatorial:
        """
        Transfers current vector to equatorial
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, _Ecliptical):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz)
        elif isinstance(self.coords, _Horizontal):
            xyz_eqt = self.__htz_to_eqt_cartesian__(xyz)
        else:
            xyz_eqt = xyz

        spherical = self.__cartesian_to_spherical__(xyz_eqt)
        return _Equatorial(spherical['horz_angle'],
                           spherical['vert_angle'])

    def horizontal(self) -> _Horizontal:
        """
        Transfers current vector to horizontal
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, _Ecliptical):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz)
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
        elif isinstance(self.coords, _Equatorial):
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz)
        else:
            xyz_hrz = xyz

        spherical = self.__cartesian_to_spherical__(xyz_hrz)
        return _Horizontal(90 - spherical['horz_angle'],
                           spherical['vert_angle'])

    def ecliptical(self) -> _Ecliptical:
        """
        Transfers current vector to ecliptical
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, _Equatorial):
            xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz)
        elif isinstance(self.coords, _Horizontal):
            xyz_eqt = self.__htz_to_eqt_cartesian__(xyz)
            xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz_eqt)
        else:
            xyz_ecl = xyz

        spherical = self.__cartesian_to_spherical__(xyz_ecl)
        return _Ecliptical(spherical['horz_angle'],
                           spherical['vert_angle'])

    def ascention_diff(self, dec: Optional[float] = None) -> float:
        """
        Returns ascention difference of <self>
        object or any given declination
        """
        decl = self.equatorial().dec if dec is None else dec
        tan_p = self.__constants__.tan_p
        tan_d = tan(decl * pi / 180)
        # If there is no ascention for a given declination on extreme latitude
        if abs(tan_p * tan_d) > 1:
            return None
        return asin(tan_p * tan_d) * 180 / pi

    def oblique_asc(self, dec: Optional[float] = None) -> float:
        """
        Returns oblique ascention of <self> object
        or any given declination
        """
        tan_p = self.__constants__.tan_p
        dec = self.equatorial().dec
        real_asc = self.equatorial().ra
        tan_d = tan(dec * pi / 180)

        # In case of no ascention
        if abs(tan_p * tan_d) > 1:
            return None
        ascention_diff = self.ascention_diff(dec)
        return (real_asc - ascention_diff) % 360

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
        _y = sin_t * tan_p * tan_e
        _x = 1 - cos_t * tan_p * tan_e
        diff = atan_btw_xy(_x, _y)

        ra_asc = (tsa + diff) % 360
        _y = sin(diff * pi / 180)
        _x = tan_p
        dec_asc = atan_btw_xy(_x, _y)

        xyz_eqt = self.cartesian(_Equatorial(ra_asc, dec_asc))
        xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz_eqt)
        possible_asc = self.__cartesian_to_spherical__(xyz_ecl)['horz_angle']

        # Now we need to ensure that asc is still on the eastern
        # horizon. If it is on the west, it is descendant (such
        # a situatin may happen when cross-point of ecliptic/horizon
        # planes travels through the South or North at extreme latitudes)
        xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
        if xyz_hrz.x < 0:
            possible_asc += 180

        return possible_asc % 360

    def medium_coeli(self):
        """
        Returns culminating Zodiac degree
        """
        ramc = self.ramc()
        _y = tan(ramc * pi / 180)
        _x = self.__constants__.cos_e
        if 90 < ramc <= 180:
            possible_mc = (atan_btw_xy(_x, _y) - 180) % 350
        elif 180 < ramc <= 270:
            possible_mc = (atan_btw_xy(_x, _y) + 180) % 350
        else:
            possible_mc = atan_btw_xy(_x, _y)

        # Now we need to ensure that asc is still on the eastern
        # horizon. If it is on the west, it is descendant (such
        # a situatin may happen when cross-point of ecliptic/horizon
        # planes travels through the South or North at extreme latitudes)
        xyz_ecl = self.cartesian(_Ecliptical(possible_mc, 0))
        xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz_ecl)
        xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
        if xyz_hrz.z < 0:
            possible_mc += 180

        return possible_mc % 360

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
        asc_diff = self.ascention_diff()
        if asc_diff is None:
            return None
        return (90 + asc_diff) % 360

    @classmethod
    def show_degrees_minutes(cls, deg: float) -> str:
        """
        Transforms float degree into degree and minute
        """
        degree = floor(deg)
        minute = (deg - degree) * 60
        return f"{degree}º {round(minute)}'"


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':

    # Initiate a vector object
    test_vector = Vector(
        datetime(2022, 10, 10, 7, 20),
        time_zone=4,
        geo_lon=44 + 46/60,
        geo_lat=76 + 43/50,
    )

    # Get celestial information
    # Local Sidereal Time and real ascention of MC
    print('LST:', timedelta(hours=test_vector.lst))
    print('RAMC:', test_vector.ramc())

    # Ascendant Zodiac degree
    asc = test_vector.asc()
    print('ASC:', asc)

    # Set coordinates in a chosen system
    test_vector.set_horizontal(azm=90.0, alt=0.0)

    # Get the same vector in different systems
    print('True East:', test_vector.cartesian())
    print('True East:', test_vector.horizontal())
    print('True East:', test_vector.equatorial())
    print('True East:', test_vector.ecliptical())
    print('True East oblique ascention:', test_vector.oblique_asc())

    # Observe cartasian coordinates in any system
    print(test_vector.cartesian(test_vector.ecliptical()))

    # Let's check ascendant
    test_vector.set_ecliptical(lon=test_vector.asc(), lat=0.0)
    print('ASC:', test_vector.horizontal())
    print('ASC:', test_vector.equatorial())

    # The same as real ascention of True East above
    print('Oblique Ascention of ASC:', test_vector.oblique_asc())

    # Uppper Meridian Distance (UMD) and
    # Diurnal Semiarc (DSA) of a given vector
    print('UMD:', test_vector.umd())
    print('DSA:', test_vector.dsa())
    print('MC:', test_vector.medium_coeli())
