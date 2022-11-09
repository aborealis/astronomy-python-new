"""
Initiates the celestial sphere, creates a point
on the sphere, and converts it to different
coordinate systems.
"""
from datetime import datetime, timedelta
from typing import Optional, Callable
from math import atan, tan, pi
from numpy import sin, cos, arcsin as asin, arctan2 as atan2, arccos as acos
from components.vector import time as tm
from components.vector import coords as crd


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


class Sphere:
    """
    Requires local time and geo-position
    of the observer at initialization.
    """

    def __init__(self,
                 naive_datetime: datetime,
                 time_zone: float,
                 geo_lon: float,
                 geo_lat: float) -> None:
        localtime = tm.__localtime__(naive_datetime, time_zone)
        localtime_utc = tm.__localtime_utc__(localtime)
        epsilon = tm.__inclination_ecliptic__(localtime_utc)

        self.lst = tm.__lst__(localtime_utc, geo_lon)
        self.__equatorial = None
        self.__equatorial_xyz = None
        self.ramc = self.lst * 15

        # Sid. time since Aries ascended (tsa, in radians)
        tsa = (self.lst * 15 + 90) / 180 * pi

        self.__constants__ = crd.Constants(
            cos_e=cos(epsilon * pi / 180),
            sin_e=sin(epsilon * pi / 180),
            tan_e=tan(epsilon * pi / 180),
            cos_p=cos(geo_lat * pi / 180),
            sin_p=sin(geo_lat * pi / 180),
            tan_p=tan(geo_lat * pi / 180),
            sin_t=sin(tsa),
            cos_t=cos(tsa),
            epsilon=epsilon,
            dec_max=90 - geo_lat,
        )

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Coordinates conversion methods
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def __ecl_to_eqt__(self, ecl: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        ecliptical to equatorial system.
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        return crd.Cartesian(
            x=ecl.x,
            y=cos_e * ecl.y - sin_e * ecl.z,
            z=sin_e * ecl.y + cos_e * ecl.z
        )

    def __eqt_to_ecl__(self, eqt: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        equatorial to ecliptical system.
        """
        sin_e = self.__constants__.sin_e
        cos_e = self.__constants__.cos_e

        return crd.Cartesian(
            x=eqt.x,
            y=cos_e * eqt.y + sin_e * eqt.z,
            z=- sin_e * eqt.y + cos_e * eqt.z
        )

    def __eqt_to_hrz__(self, eqt: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        equatorial to horizontal system.
        """
        # Sid. time since Aries ascended (tsa, in radians)
        sin_t = self.__constants__.sin_t
        cos_t = self.__constants__.cos_t
        cos_p = self.__constants__.cos_p
        sin_p = self.__constants__.sin_p

        # Rotated system with (1, 0, 0) pointing to (ra=0, dec=0)
        x_eqt_rot = cos_t * eqt.x + sin_t * eqt.y
        y_eqt_rot = - sin_t * eqt.x + cos_t * eqt.y
        z_eqt_rot = eqt.z

        return crd.Cartesian(
            x=x_eqt_rot,
            y=sin_p * y_eqt_rot + cos_p * z_eqt_rot,
            z=-cos_p * y_eqt_rot + sin_p * z_eqt_rot
        )

    def __hrz_to_eqt__(self, hrz: crd.Cartesian) -> crd.Cartesian:
        """
        Transfers cartesian (x, y, z) from
        horizontal to equatorial system.
        """
        # Sid. time since Aries ascended (tsa, in radians)
        sin_t = self.__constants__.sin_t
        cos_t = self.__constants__.cos_t
        cos_p = self.__constants__.cos_p
        sin_p = self.__constants__.sin_p

        # rotated vector with (1, 0, 0) pointing to (ra=0, dec=0)
        x_eqt_rot = hrz.x
        y_eqt_rot = sin_p * hrz.y - cos_p * hrz.z
        z_eqt_rot = cos_p * hrz.y + sin_p * hrz.z

        return crd.Cartesian(
            x=cos_t * x_eqt_rot - sin_t * y_eqt_rot,
            y=sin_t * x_eqt_rot + cos_t * y_eqt_rot,
            z=z_eqt_rot)

    @staticmethod
    def __cartesian_to_spherical__(xyz: crd.Cartesian) -> list:
        """
        Converts cartesian (x, y, z) coordinates
        into spherical angles.
        """
        _x, _y, _z = xyz.aslist()
        lon = (atan2(_y, _x) / pi * 180) % 360  # xy plane
        lat = asin(_z) / pi * 180  # z plane

        return [lon, lat]

    @staticmethod
    def __spherical_to_cartesian__(lon: float, lat: float) -> crd.Cartesian:
        """
        Converts spherical angles into
        cartesian (x, y, z).
        """
        cos_lat = cos(lat / 180 * pi)
        sin_lat = sin(lat / 180 * pi)
        cos_lon = cos(lon / 180 * pi)
        sin_lon = sin(lon / 180 * pi)

        return crd.Cartesian(
            x=cos_lat * cos_lon,
            y=cos_lat * sin_lon,
            z=sin_lat,
        )

    def __set_vector(self) -> crd.Vector:
        """
        Creates a vector objects from the <self> object
        which has an .__equatorial property.
        """
        equatorial = self.__equatorial
        equatorial_xyz = self.__equatorial_xyz
        rasc = self.__equatorial.rasc
        dec = self.__equatorial.dec

        # Instead of passing the links to the outer's class
        # methods, we will construct new methods which use
        # context environment. Otherwise each time the new
        # vector object is initiated, the functions below
        # will be called with a new .self.__equatorial property
        # for all previously created vectors.
        return crd.Vector(
            equatorial=lambda: self.__get_equatorial(equatorial),
            ecliptical=lambda: self.__get_ecliptical(equatorial_xyz),
            horizontal=lambda: self.__get_horizontal(equatorial_xyz),
            equatorial_xyz=lambda: self.__get_equatorial_xyz(equatorial_xyz),
            ecliptical_xyz=lambda: self.__get_ecliptical_xyz(equatorial_xyz),
            horizontal_xyz=lambda: self.__get_horizontal_xyz(equatorial_xyz),
            dsa=lambda: self.dsa(dec),
            umd=lambda: self.umd(rasc),
            obl_asc=lambda: self.obl_asc(rasc, dec),
            asc_diff=lambda: self.asc_diff(dec),
        )

    def set_equatorial(self, rasc: float, dec: float) -> crd.Vector:
        """
        Sets equatorial spherical coordinates.
        """
        self.__equatorial = crd.Equatorial(rasc, dec)
        self.__equatorial_xyz = self.__spherical_to_cartesian__(rasc, dec)
        return self.__set_vector()

    def set_ecliptical(self, lon: float, lat: float) -> crd.Vector:
        """
        Sets ecliptical spherical coordinates.
        """
        xyz_ecl = self.__spherical_to_cartesian__(lon, lat)
        xyz_eqt = self.__ecl_to_eqt__(xyz_ecl)
        coords = self.__cartesian_to_spherical__(xyz_eqt)
        self.__equatorial = crd.Equatorial(*coords)
        self.__equatorial_xyz = xyz_eqt
        return self.__set_vector()

    def set_horizontal(self, azm: float, alt: float) -> crd.Vector:
        """
        Sets horizontal spherical coordinates.
        """
        xyz_hrz = self.__spherical_to_cartesian__(90 - azm, alt)
        xyz_eqt = self.__hrz_to_eqt__(xyz_hrz)
        coords = self.__cartesian_to_spherical__(xyz_eqt)
        self.__equatorial = crd.Equatorial(*coords)
        self.__equatorial_xyz = xyz_eqt
        return self.__set_vector()

    def __get_equatorial_xyz(self,
                             equatorial_xyz: Optional[crd.Cartesian] = None
                             ) -> crd.Cartesian:
        """
        Converts current vector to equatorial cartesian coordinates.
        """
        if equatorial_xyz is None:
            return self.__equatorial_xyz
        return equatorial_xyz

    def __get_equatorial(self,
                         equatorial: Optional[crd.Equatorial] = None
                         ) -> crd.Equatorial:
        """
        Converts current vector to equatorial spherical coordinates.
        """
        if equatorial:
            return equatorial
        return self.__equatorial

    def __get_ecliptical_xyz(self,
                             equatorial_xyz: Optional[crd.Cartesian] = None
                             ) -> crd.Cartesian:
        """
        Converts current vector to ecliptical cartesian coordinates.
        """
        xyz_eqt = self.__equatorial_xyz if equatorial_xyz is None else equatorial_xyz
        return self.__eqt_to_ecl__(xyz_eqt)

    def __get_ecliptical(self,
                         equatorial_xyz: Optional[crd.Cartesian] = None
                         ) -> crd.Ecliptical:
        """
        Converts current vector to ecliptical spherical coordinates.
        """
        xyz_ecl = self.__get_ecliptical_xyz(equatorial_xyz)
        spherical = self.__cartesian_to_spherical__(xyz_ecl)
        return crd.Ecliptical(*spherical)

    def __get_horizontal_xyz(self,
                             equatorial_xyz: Optional[crd.Cartesian] = None
                             ) -> crd.Cartesian:
        """
        Converts current vector to horizontal spherical coordinates.
        """
        xyz_eqt = self.__equatorial_xyz if equatorial_xyz is None else equatorial_xyz
        return self.__eqt_to_hrz__(xyz_eqt)

    def __get_horizontal(self,
                         equatorial_xyz: Optional[crd.Cartesian] = None
                         ) -> crd.Horizontal:
        """
        Converts current vector to horizontal spherical coordinates.
        """
        xyz_hrz = self.__get_horizontal_xyz(equatorial_xyz)
        spherical = self.__cartesian_to_spherical__(xyz_hrz)
        return crd.Horizontal(90 - spherical[0], spherical[1])

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Vector-specific methods
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def asc_diff(self, dec: Optional[float] = None) -> float:
        """
        Returns ascension difference for a given declination
        """
        decl = dec if dec is not None else self.__equatorial.dec
        if decl is None:
            return None
        tan_p = self.__constants__.tan_p
        tan_d = tan(decl * pi / 180)

        # In case of no ascension at extreme latitudes
        if abs(tan_p * tan_d) > 1:
            return None
        return asin(tan_p * tan_d) * 180 / pi

    def obl_asc(self,
                rasc: Optional[float] = None,
                dec: Optional[float] = None) -> float:
        """
        Returns oblique ascension of a given point
        """
        r_asc = rasc if rasc is not None else self.__equatorial.rasc
        decl = dec if dec is not None else self.__equatorial.dec
        if None in [r_asc, decl]:
            return None
        ascension_diff = self.asc_diff(decl)

        # In case of no ascension at extreme latitudes
        if ascension_diff is None:
            return None
        return (r_asc - ascension_diff) % 360

    def umd(self, rasc: Optional[float] = None) -> float:
        """
        Returns upper meridian distance, UMD
        for a given RA
        """
        r_asc = rasc if rasc is not None else self.__equatorial.rasc
        if r_asc is None:
            return None
        return true_distance(r_asc, self.ramc)

    def dsa(self, dec: Optional[float] = None) -> float:
        """
        Returns diurnal semiarc for a given declination
        """
        decl = dec if dec is not None else self.__equatorial.dec
        if decl is None:
            return None
        asc_diff = self.asc_diff(decl)

        # In case of no ascension at extreme latitudes
        if asc_diff is None:
            return None
        return (90 + asc_diff) % 360

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Sphere properties methods
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    @property
    def asc(self):
        """
        Returns ascending Zodiac degree
        """
        # time since Aries raising (tsa, in degrees)
        tsa = (self.lst * 15 + 90) % 360
        cos_t = self.__constants__.cos_t
        sin_t = self.__constants__.sin_t
        tan_p = self.__constants__.tan_p
        tan_e = self.__constants__.tan_e
        _y = sin_t * tan_p * tan_e
        _x = 1 - cos_t * tan_p * tan_e
        diff = (atan2(_y, _x) * 180 / pi) % 360

        # Alternative formula
        # tan_ra = -cos(self.ramc * pi / 180) / (
        #     sin(self.ramc * pi / 180) + tan_e*tan_p
        # )
        # test_ra = atan(tan_ra) * 180 / pi % 360

        ra_asc = (tsa + diff) % 360
        _y = sin(diff * pi / 180)
        _x = tan_p
        dec_asc = (atan2(_y, _x) * 180 / pi) % 360

        xyz_eqt = self.__spherical_to_cartesian__(ra_asc, dec_asc)
        xyz_ecl = self.__eqt_to_ecl__(xyz_eqt)
        possible_asc = self.__cartesian_to_spherical__(xyz_ecl)[0]

        # Now we need to ensure that ASC is still on the eastern horizon.
        # If it is on the west, it is DSC (such a situation may happen when
        # a cross-point of ecliptic/horizon planes travels through the True
        # South/North at extreme latitudes)
        xyz_hrz = self.__eqt_to_hrz__(xyz_eqt)
        if xyz_hrz.x < 0:
            possible_asc += 180

        return possible_asc % 360

    @property
    def medium_coeli(self):
        """
        Returns culminating Zodiac degree
        """
        ramc = self.ramc
        _y = tan(ramc * pi / 180)
        _x = self.__constants__.cos_e
        atn = (atan2(_y, _x) * 180 / pi) % 360
        if 90 < ramc <= 180:
            possible_mc = (atn - 180) % 350
        elif 180 < ramc <= 270:
            possible_mc = (atn + 180) % 350
        else:
            possible_mc = atn

        # Now we need to ensure that MC is still above the horizon. If
        # it is below, it is IC (such a situation may happen when a cross-
        # point of ecliptic/horizon planes travels through the TrueSouth/
        # North at extreme latitudes)
        xyz_ecl = self.__spherical_to_cartesian__(possible_mc, 0)
        xyz_eqt = self.__ecl_to_eqt__(xyz_ecl)
        xyz_hrz = self.__eqt_to_hrz__(xyz_eqt)
        if xyz_hrz.z < 0:
            possible_mc += 180

        return possible_mc % 360

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Units conversion methods:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

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
    # Astrology houses related methods:
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

    @classmethod
    def bisect_root(cls, func: Callable, x_min: float, x_max: float) -> float:
        """
        Finds the root of the function in [x_min, x_max]
        """
        # If no single root in the interval
        if func(x_min) * func(x_max) > 0:
            return None

        while True:
            x_center = x_min + (x_max - x_min) * 0.618
            if (func(x_min) * func(x_center)) > 0:
                x_min = x_center
            else:
                x_max = x_center
            if (x_max - x_min) < 1e-10:
                break
        return x_center

    def placidus(self, cusp: int) -> list[float]:
        """
        Returns zodiac degree of the cusp in Placidus system.
        Since in subpolar extreme latitudes there may be more than
        one cusp for 9th and 11th houses, the function returns
        the list of cusp for a given house.
        """
        if cusp in [4, 7, 5, 6, 2, 3]:
            opposite_cusp = 12 if cusp == 6 else (cusp + 6) % 12
            opposite_cusp = self.placidus(opposite_cusp)
            if not opposite_cusp:
                return []
            return [(cusp + 180) % 360 for cusp in opposite_cusp]
        if cusp == 1:
            return [self.asc]
        if cusp == 10:
            return [self.medium_coeli]

        alpha = 1/3 if cusp in [11, 9] else 2/3
        if cusp in [8, 9]:
            alpha = -alpha
        tan_p = self.__constants__.tan_p
        tan_e = self.__constants__.tan_e

        def func(_x):
            """
            Equation of intersection of placiduc curve and
            ecliptic. The root of the equestion _x is the
            tan(dec) of intersecting point.
            """
            asc_diff = asin(_x * tan_p)
            dsa = pi / 2 + asc_diff
            r_asc = self.ramc * pi / 180 + dsa * alpha
            return sin(r_asc) * tan_e - _x

        # We cannot exceed this declination while search for
        # root of the function
        decl_max = min(self.__constants__.dec_max,
                       self.__constants__.epsilon)
        x_max = tan(decl_max * pi / 180)

        # Since we do not knoe the regions with roots,
        # we will split interval [-x_max, x_max] into
        # this number of subintervals
        number_of_search_intervals = 8
        search_intervals = [
            [-x_max + _i * 2 * x_max/number_of_search_intervals,
             -x_max + (_i + 1) * 2 * x_max/number_of_search_intervals
             ]
            for _i in range(number_of_search_intervals)
        ]

        # Make adjustments on the edges of
        # [-x_max, x_max] interval to avoid uncertainty
        # on the edges at extreme latitudes
        search_intervals[0][0] += 1e-10
        search_intervals[len(search_intervals) - 1][1] -= 1e-10

        roots = []
        for interval in search_intervals:
            root = self.bisect_root(
                func, interval[0], interval[1]
            )
            if root is not None:
                roots.append(root)

        if not roots:
            return []

        cusps = []
        for root in roots:
            dec_cusp = atan(root) * 180 / pi
            ra_cusp = self.__zodiac_ra(dec_cusp) % 360

            # Check which ra to choose
            part_of_dsa = abs(self.dsa(dec_cusp) * alpha)
            merid_dist = true_distance(ra_cusp, self.ramc)
            if abs(merid_dist - part_of_dsa) > 1e-3:
                ra_cusp = 180 - ra_cusp

            # Converts RA and Dec into lon and lat
            xyz_eqt = self.__spherical_to_cartesian__(
                ra_cusp,
                dec_cusp,
            )
            xyz_ecl = self.__eqt_to_ecl__(xyz_eqt)
            cusps.append(
                self.__cartesian_to_spherical__(xyz_ecl)[0]
            )
        return cusps

    def regiomontanus(self, cusp: int) -> float:
        """
        Returns zodiac degree of the cusp in Regiomontanus system
        """
        if cusp in [4, 7, 5, 6, 2, 3]:
            opposite_cusp = 12 if cusp == 6 else (cusp + 6) % 12
            return (self.regiomontanus(opposite_cusp) + 180) % 360
        if cusp == 1:
            return self.asc
        if cusp == 10:
            return self.medium_coeli

        tan_p = self.__constants__.tan_p
        tan_e = self.__constants__.tan_e
        cos_e = self.__constants__.cos_e

        oa_asc = self.ramc + 90
        cos_oa = cos(oa_asc * pi / 180)
        sin_oa = sin(oa_asc * pi / 180)

        ra_w = self.ramc + (cusp - 10) * 30
        tan_ra_w = tan(ra_w * pi / 180)
        ra_w *= pi / 180

        _y = tan_ra_w
        _x = 1 - tan_e * tan_p * (
            cos_oa + sin_oa * tan_ra_w
        )

        possible_cusp = atan(_y/_x / cos_e) * 180 / pi % 360
        xyz_ecl = self.__spherical_to_cartesian__(possible_cusp, 0)
        xyz_eqt = self.__ecl_to_eqt__(xyz_ecl)
        xyz_hor = self.__eqt_to_hrz__(xyz_eqt)
        if xyz_hor.z < 0:
            possible_cusp += 180

        return possible_cusp % 360


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':

    # Initiate a celestial sphere object
    sphere = Sphere(
        datetime(2022, 10, 10, 21, 20),
        time_zone=4,
        geo_lon=44 + 46/60,
        geo_lat=76 + 43/60,
    )

    # You can access properties of the sphere
    print('Sphere properties:')
    print('- LST (hours):', sphere.lst)
    print('- LST (h:m:s):', timedelta(hours=sphere.lst))
    print('- RAMC:', sphere.ramc)
    print('- ASC:', sphere.asc)
    print('- MC:', sphere.medium_coeli)

    # True East
    true_east = sphere.set_horizontal(azm=90, alt=0)
    print('\nTrue East:')
    print('-', true_east.horizontal_xyz())
    print('-', true_east.horizontal())  # Alt = 0 for True East
    print('-', true_east.equatorial())
    print('-', true_east.ecliptical())
    print('- right ascension:', true_east.equatorial().rasc)
    print('- oblique ascension', true_east.obl_asc())  # OA = RA for True Esat

    # Ascendant
    asc = sphere.set_ecliptical(lon=sphere.asc, lat=0)
    print('\nASC:')
    print('-', asc.equatorial())
    print('-', asc.horizontal())  # Altitude = 0 for ASC
    print('-', asc.ecliptical())
    print('- oblique ascension:', asc.obl_asc())
    print('- ascension diff:', asc.asc_diff())
    print('- UMD:', asc.umd())
    print('- DSA:', asc.dsa())  # UMD = DSA for ASC

    # MC
    mc = sphere.set_ecliptical(lon=sphere.medium_coeli, lat=0)
    print('MC')
    print('-', mc.horizontal())
    print('-', mc.horizontal_xyz())

    # Eleventh house cusp
    cuspid = sphere.placidus(11)
    print('\nH11 (Placidus):')
    print('-', cuspid)
    if cuspid:
        h11 = sphere.set_ecliptical(lon=cuspid[0], lat=0)
        print('- UMD:', h11.umd())
        print('- DSA/3:', h11.dsa() / 3)  # UMD = DSA/3 for H11

    cuspid = sphere.regiomontanus(11)
    print('\nH11 (Regiomontanus):')
    print('-', cuspid)

    # Aldebaran star
    aldebaran = sphere.set_equatorial(
        rasc=sphere.time_to_angle(hours=4, minutes=35, seconds=54.32),
        dec=16 + 30/60 + 26.1/3600
    )
    print('aldebaran', aldebaran.equatorial())
