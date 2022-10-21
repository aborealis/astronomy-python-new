"""
Creates vector and transfers it between
spherical coordinate systems
"""
from typing import Optional, Union
from datetime import datetime, timedelta
from math import asin, atan
from numpy import cos, sin, pi, tan
from components.vector import time as tm
from components.vector import coords as crd


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


class Vector:
    """
    Requires local time and geo-position
    of the observer at initialization
    """

    def __init__(self,
                 naive_datetime: datetime,
                 time_zone: float,
                 geo_lon: float,
                 geo_lat: float) -> None:
        localtime = tm.__localtime__(naive_datetime, time_zone)
        localtime_utc = tm.__localtime_utc__(localtime)
        epsilon = tm.__inclination_ecliptic__(localtime_utc)

        self.__constants__ = crd.Constants(
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

    def set_equatorial(self, rasc: float, dec: float) -> crd.Equatorial:
        """
        Sets equatorial coordinates
        """
        self.coords = crd.Equatorial(rasc, dec)

    def set_ecliptical(self, lon: float, lat: float) -> crd.Ecliptical:
        """
        Sets equatorial coordinates
        """
        self.coords = crd.Ecliptical(lon, lat)

    def set_horizontal(self, azm: float, alt: float) -> crd.Horizontal:
        """
        Sets horizontal coordinates
        """
        self.coords = crd.Horizontal(azm, alt)

    def cartesian(self,
                  vector: Union[
                      crd.Ecliptical,
                      crd.Equatorial,
                      crd.Horizontal,
                      None,
                  ] = None) -> crd.Cartesian:
        """
        Converts spherical coordinates
        into cartesian (x, y, z).

        Directions:
        - Equatorial(ra=0, dec=0) -> (1, 0, 0)
        - Ecliptical(lon=0, lat=0) -> (1, 0, 0)
        - Horizontal(azm=0, alt=0) -> (0, 1, 0)
        """
        coords = self.coords if vector is None else vector
        if isinstance(coords, crd.Equatorial):
            lon = coords.rasc
            lat = coords.dec
        elif isinstance(coords, crd.Ecliptical):
            lon = coords.lon
            lat = coords.lat
        else:
            lon = 90 - coords.azm
            lat = coords.alt

        cos_lat = cos(lat / 180 * pi)
        sin_lat = sin(lat / 180 * pi)
        cos_lon = cos(lon / 180 * pi)
        sin_lon = sin(lon / 180 * pi)

        return crd.Cartesian(
            x=cos_lat * cos_lon,
            y=cos_lat * sin_lon,
            z=sin_lat,
        )

    def cartesian_horizontal(self):
        """
        This function is created for numpy
        visualization only. It is completely
        equivalent to self.cartesian(self.horizontal())
        """
        cartesian = self.cartesian()
        coords = self.coords
        if isinstance(coords, crd.Equatorial):
            return self.__eqt_to_hrz_cartesian__(cartesian)
        if isinstance(coords, crd.Ecliptical):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(cartesian)
            return self.__eqt_to_hrz_cartesian__(xyz_eqt)
        return cartesian

    @staticmethod
    def __cartesian_to_spherical__(coords: crd.Cartesian) -> dict:
        """
        Converts cartesian (x, y, z) coordinates
        into spherical angles
        """
        if coords.x == 0:
            lon = 270
        else:
            lon = atan_btw_xy(coords.x, coords.y)
        lat = asin(coords.z) / pi * 180

        return dict(horz_angle=lon, vert_angle=lat)

    def __ecl_to_eqt_cartesian__(self, ecl: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        ecliptical to equatorial system
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        x_eqt = ecl.x
        y_eqt = cos_e * ecl.y - sin_e * ecl.z
        z_eqt = sin_e * ecl.y + cos_e * ecl.z

        return crd.Cartesian(x_eqt, y_eqt, z_eqt)

    def __eqt_to_ecl_cartesian__(self, eqt: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        equatorial to ecliptical system
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        x_ecl = eqt.x
        y_ecl = cos_e * eqt.y + sin_e * eqt.z
        z_ecl = - sin_e * eqt.y + cos_e * eqt.z

        return crd.Cartesian(x_ecl, y_ecl, z_ecl)

    def __eqt_to_hrz_cartesian__(self, eqt: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        equatorial to horizontal system
        """
        # Sid. time since Aries ascended (tsa, in radians)
        tsa = (self.lst * 15 + 90) / 180 * pi
        sin_t = sin(tsa)
        cos_t = cos(tsa)
        cos_p = self.__constants__.cos_p
        sin_p = self.__constants__.sin_p

        # Rotated vector with (1, 0, 0) pointing to (ra=0, dec=0)
        x_eqt_rot = cos_t * eqt.x + sin_t * eqt.y
        y_eqt_rot = - sin_t * eqt.x + cos_t * eqt.y
        z_eqt_rot = eqt.z

        x_hrz = x_eqt_rot
        y_hrz = sin_p * y_eqt_rot + cos_p * z_eqt_rot
        z_hrz = - cos_p * y_eqt_rot + sin_p * z_eqt_rot

        return crd.Cartesian(x_hrz, y_hrz, z_hrz)

    def __htz_to_eqt_cartesian__(self, hrz: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        horizontal to equatorial system
        """
        # Sid. time since Aries ascended (tsa, in radians)
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

        return crd.Cartesian(x_eqt, y_eqt, z_eqt)

    def equatorial(self) -> crd.Equatorial:
        """
        Transfers current vector (or any other
        given vector) to equatorial coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, crd.Ecliptical):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz)
        elif isinstance(self.coords, crd.Horizontal):
            xyz_eqt = self.__htz_to_eqt_cartesian__(xyz)
        else:
            xyz_eqt = xyz

        spherical = self.__cartesian_to_spherical__(xyz_eqt)
        return crd.Equatorial(spherical['horz_angle'],
                              spherical['vert_angle'])

    def horizontal(self) -> crd.Horizontal:
        """
        Transfers current vector to horizontal
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, crd.Ecliptical):
            xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz)
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
        elif isinstance(self.coords, crd.Equatorial):
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz)
        else:
            xyz_hrz = xyz

        spherical = self.__cartesian_to_spherical__(xyz_hrz)
        return crd.Horizontal(90 - spherical['horz_angle'],
                              spherical['vert_angle'])

    def ecliptical(self) -> crd.Ecliptical:
        """
        Transfers current vector to ecliptical
        coordinate system
        """
        xyz = self.cartesian()
        if isinstance(self.coords, crd.Equatorial):
            xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz)
        elif isinstance(self.coords, crd.Horizontal):
            xyz_eqt = self.__htz_to_eqt_cartesian__(xyz)
            xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz_eqt)
        else:
            xyz_ecl = xyz

        spherical = self.__cartesian_to_spherical__(xyz_ecl)
        return crd.Ecliptical(spherical['horz_angle'],
                              spherical['vert_angle'])

    def ascension_diff(self, dec: Optional[float] = None) -> float:
        """
        Returns ascension difference of <self>
        object or any given declination
        """
        decl = self.equatorial().dec if dec is None else dec
        tan_p = self.__constants__.tan_p
        tan_d = tan(decl * pi / 180)

        # In case of no ascension at extreme latitudes
        if abs(tan_p * tan_d) > 1:
            return None
        return asin(tan_p * tan_d) * 180 / pi

    def oblique_asc(self, dec: Optional[float] = None) -> float:
        """
        Returns oblique ascension of <self> object
        or any given declination
        """
        tan_p = self.__constants__.tan_p
        dec = self.equatorial().dec
        rasc = self.equatorial().rasc
        tan_d = tan(dec * pi / 180)

        # In case of no ascension at extreme latitudes
        if abs(tan_p * tan_d) > 1:
            return None
        ascension_diff = self.ascension_diff(dec)
        return (rasc - ascension_diff) % 360

    def ramc(self):
        """
        Right ascension of Medium Coeli for
        given time of <self> object
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

        xyz_eqt = self.cartesian(crd.Equatorial(ra_asc, dec_asc))
        xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz_eqt)
        possible_asc = self.__cartesian_to_spherical__(xyz_ecl)['horz_angle']

        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        # Now we need to ensure that ASC
        # is still on the eastern horizon.
        # If it is on the west, it is DSC
        # (such a situation may happen when
        # a cross-point of ecliptic/horizon
        # planes travels through the True
        # South/North at extreme latitudes)
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
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

        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        # Now we need to ensure that MC
        # is still above the horizon. If
        # it is below, it is IC (such a
        # situation may happen when a cross-
        # point of ecliptic/horizon planes
        # travels through the TrueSouth/
        # North at extreme latitudes)
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        xyz_ecl = self.cartesian(crd.Ecliptical(possible_mc, 0))
        xyz_eqt = self.__ecl_to_eqt_cartesian__(xyz_ecl)
        xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
        if xyz_hrz.z < 0:
            possible_mc += 180

        return possible_mc % 360

    def umd(self, rasc: Optional[float] = None) -> float:
        """
        Returns upper meridian distance, UMD
        for <self> object or a given RA
        """
        r_asc = self.equatorial().rasc if rasc is None else rasc
        return true_distance(r_asc, self.ramc())

    def dsa(self, dec: Optional[float] = None) -> float:
        """
        Returns diurnal semiarc for <self> object or
        a given declination
        """
        asc_diff = self.ascension_diff(dec)
        if asc_diff is None:
            return None
        return (90 + asc_diff) % 360

    @classmethod
    def dms(cls, deg: float) -> str:
        """
        Converts float degree into degrees, minutes and seconds
        """
        return crd.degrees_to_dms(deg)

    @classmethod
    def time_to_angle(cls, hours=0, minutes=0,
                      seconds=0,  milliseconds=0,
                      microseconds=0) -> float:
        """
        Converts timedelta into angle in degrees
        """
        delta = timedelta(
            hours=hours, minutes=minutes,
            seconds=seconds, milliseconds=milliseconds,
            microseconds=microseconds
        )
        return delta.total_seconds() / 3600 * 15

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Private methods:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    def __zodiac_ra(self, dec: float) -> float:
        """
        Returns zodiac RA fron the declination
        """
        tan_d = tan(dec * pi / 180)
        tan_e = self.__constants__.tan_e
        return asin(tan_d / tan_e) * 180 / pi

    def __zodiac_dec(self, rasc: float) -> float:
        """
        Returns zodiac declination from the RA
        """
        tan_e = self.__constants__.tan_e
        sin_ra = sin(rasc * pi / 180)
        return atan(sin_ra * tan_e) * 180 / pi

    def __regions_of_zodiac_ascension(self, eastern: bool = True) -> list[list]:
        """
        Returns ranges in RAs of the eastern/western horizon where
        Zodiac degrees are always ascending for a given latitude.
        There could be a few regions with always ascending degrees
        at extreme latitudes, so we return them all.
        """
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        # The eastern part of equator goes
        # from RAMC to RAIC and the eastern
        # one goes from RAIC to RAMC
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        ramc = self.ramc()
        start, finish = ramc, (ramc + 180) % 360
        if not eastern:
            start, finish = finish, start

        if self.__constants__.epsilon < self.__constants__.dec_max:
            # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
            # In case of moderate latitudes
            # there is only one region. It is
            # the whole eatern/western part
            # of the equator.
            # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
            return [[start, finish]]

        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        # Find 4 intersections between
        # ecliptical plane and hourly
        # circles at +-dec_max declination
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        dec_max = self.__constants__.dec_max
        point0 = self.__zodiac_ra(dec_max)
        intersections = [
            point0, 180 + point0, 180 - point0, 360 - point0,
        ]
        intersections = [
            crd.Equatorial(_ra, self.__zodiac_dec(_ra))
            for _ra in intersections
        ]

        # Choose only eastern/western intersections
        regions = []
        for point in intersections:
            xyz_eqt = self.cartesian(point)
            xyz_hrz = self.__eqt_to_hrz_cartesian__(xyz_eqt)
            if eastern and xyz_hrz.x > 0:
                regions.append(point.rasc)
            elif not eastern and xyz_hrz.x < 0:
                regions.append(point.rasc)

        # Order intersection by RA distance from RAMC
        regions = [start, *regions, finish]
        regions.sort(key=lambda item: true_distance(
            item, start), reverse=not eastern)

        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        # We will get rid of the situation
        # when adjacent degrees differ from
        # each other by more than 180 degrees.
        # It may happen arond the edge between
        # 360th and 0th degrees of the equator.
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        for _i in range(len(regions) - 1):
            if abs(regions[_i] - regions[_i+1]) > 180:
                if regions[_i] < 180:
                    regions[_i] += 360
                else:
                    regions[_i+1] += 360

        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        # Find the region with always-
        # ascending zodiac degrees. Zodiac
        # degree in the middle of such a
        # region will always have a "safe"
        # declination within [-d_max, d_max]
        # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
        always_ascending = []
        for _i in range(len(regions) - 1):
            region_center = (regions[_i] + regions[_i + 1]) / 2
            dec = self.__zodiac_dec(region_center)
            if abs(dec) < dec_max:
                if eastern:
                    always_ascending.append([regions[_i], regions[_i + 1]])
                else:
                    always_ascending.append([regions[_i + 1], regions[_i]])

        return always_ascending

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Astrology related methods:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    def placidus(self, cusp: int) -> float:
        """
        Returns zodiac degree of the cusp in Placidus system
        """
        if cusp in [4, 7, 5, 6, 2, 3]:
            opposite_cusp = 12 if cusp == 6 else (cusp + 6) % 12
            opposite_cusp = self.placidus(opposite_cusp)
            if opposite_cusp is None:
                return None
            return (opposite_cusp + 180) % 360
        if cusp == 1:
            return self.asc()
        if cusp == 10:
            return self.medium_coeli()

        alpha = 1/3 if cusp in [11, 9] else 2/3
        if cusp in [11, 12]:
            regions = self.__regions_of_zodiac_ascension()
        else:
            regions = self.__regions_of_zodiac_ascension(False)

        def bisected(rasc: float) -> float:
            """
            A function for bisection method
            """
            dec = self.__zodiac_dec(rasc)
            dsa = self.dsa(dec)
            umd = self.umd(rasc)
            return umd - dsa * alpha

        for region in regions:
            ra_min, ra_max = region
            if ra_max < ra_min:
                ra_max += 360

            # To avoid uncertainty on the edges
            ra_min += 1e-10
            ra_max -= 1e-10
            edges = [ra_min, ra_max]

            while True:
                ra_center = (ra_max + ra_min) / 2
                if (bisected(ra_min) * bisected(ra_center)) > 0:
                    ra_min = ra_center
                else:
                    ra_max = ra_center
                if (ra_max - ra_min) < 1e-10:
                    break

            # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
            # If the asymptotic bisectional method
            # ends up with the right/left side of
            # the initial region, it means there is
            # no solution in that region - we skip
            # it. Otherwise, we return the result.
            # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
            if not (
                abs(ra_max - edges[1]) < 1e-10 or
                abs(ra_min - edges[0]) < 1e-10
            ):
                # Converts RA and Dec into lon and lat
                xyz_eqt = self.cartesian(crd.Equatorial(
                    rasc=ra_center,
                    dec=self.__zodiac_dec(ra_center)))
                xyz_ecl = self.__eqt_to_ecl_cartesian__(xyz_eqt)
                return self.__cartesian_to_spherical__(xyz_ecl)['horz_angle']
        return None


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':

    # Initiate a vector object
    test_vector = Vector(
        datetime(2022, 10, 10, 4, 20),
        time_zone=4,
        geo_lon=44 + 46/60,
        geo_lat=66 + 43/60,
    )

    # Get celestial information
    # Local Sidereal Time and right ascension of MC
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
    print('True East oblique ascension:', test_vector.oblique_asc())

    # Observe cartasian coordinates in any system
    print(test_vector.cartesian(test_vector.horizontal()))

    # Let's check ascendant
    test_vector.set_ecliptical(lon=test_vector.asc(), lat=0.0)
    print('ASC:', test_vector.horizontal())
    print('ASC:', test_vector.equatorial())

    # The same as right ascension of True East above
    print('Oblique ascension of ASC:', test_vector.oblique_asc())

    # Uppper Meridian Distance (UMD) and
    # Diurnal Semiarc (DSA) of a given vector
    print('UMD:', test_vector.umd())
    print('DSA:', test_vector.dsa())
    print('MC:', test_vector.medium_coeli())

    # 11th house cusp in Placidus system
    cuspid = test_vector.placidus(11)
    print('H11:', cuspid)
    print('H5:', test_vector.placidus(5))
    test_vector.set_ecliptical(cuspid, 0)
    if cuspid is not None:
        print('UMD:', test_vector.umd())
        print('DSA/3:', test_vector.dsa()/3)

    # Set any star by known RA expressed in time units
    aldebaran_ra = test_vector.time_to_angle(
        hours=4, minutes=35, seconds=55
    )
    print('Aldebaran right ascension:', test_vector.dms(aldebaran_ra))
