"""
Test for astronomical calculations
"""

from datetime import datetime
from typing import NamedTuple
import components.vector.time as tm
from sphere import Sphere
import tests.sweph as swe


class Data(NamedTuple):
    """
    King Charles' birt data
    """
    naive_datetime: datetime = datetime(1948, 11, 14, 21, 14, 39)
    time_zone: float = 0
    geo_lon: float = -1/6
    geo_lat: float = 51.5


data = Data()
localtime_utc = tm.__localtime_utc__(
    tm.__localtime__(data.naive_datetime, data.time_zone)
)
all_planets = swe.scan_planets(localtime_utc)
all_houses_placidus = swe.scan_houses(
    localtime_utc,
    data.geo_lon,
    data.geo_lat,
    system='P'
)

all_houses_regio = swe.scan_houses(
    localtime_utc,
    data.geo_lon,
    data.geo_lat,
    system='R'
)

sphere = Sphere(data.naive_datetime,
                data.time_zone,
                data.geo_lon,
                data.geo_lat)


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Assertions:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

def test_julian_time_is_correct():
    """
    Compare self-programmed jt with one from the swisseph module
    """
    assert tm.__julian_day__(localtime_utc) == swe.julian_day(localtime_utc)


def test_asc_is_on_horizon():
    """
    ASC should always be on the horizon
    """
    asc = sphere.set_ecliptical(sphere.asc, 0)
    assert abs(asc.horizontal().alt) < 1e-10


def test_mc_is_on_prime_meridian():
    """
    MC should always be on the prime meridian
    """
    _mc = sphere.set_ecliptical(sphere.medium_coeli, 0)
    assert abs(_mc.horizontal_xyz().x) < 1e-10
    assert _mc.horizontal_xyz().z >= 0


def test_placidus_5th_house_cusp():
    """
    Compare my calculations with swisseph module
    for placidus house system
    """
    diff = abs(sphere.placidus(5)[0] - all_houses_placidus[4])
    assert diff < 1e-2


def test_regiomontanus_5th_house_cusp():
    """
    Compare my calculations with swisseph module
    for placidus house system
    """
    cusp5 = sphere.regiomontanus(5)
    cusp5swe = all_houses_regio[4]
    diff = abs(cusp5 - cusp5swe)
    assert diff < 1e-2


def test_placidus_5th_lmd_is_third_of_nsa():
    """
    lower meridian distance for 5th house should be
    one third of nocturnal semiarc
    """
    cusp5 = sphere.set_ecliptical(sphere.placidus(5)[0], 0)
    nsa = 180 - cusp5.dsa()
    lmd = 180 - cusp5.umd()
    assert abs(lmd - nsa / 3) < 1e-9


def test_pluto_equatorial_coordinates():
    """
    Compare Pluto's coordinates with results of Bob Makransky
    """
    pluto = sphere.set_ecliptical(
        all_planets[9].lon,
        all_planets[9].lat
    ).equatorial()

    assert abs(pluto.rasc - 141.49) < 1e-1
    assert abs(pluto.dec - 23 - 4/60) < 1e-2


def test_sun_asc_diff():
    """
    Compare Pluto's asc difference with result of Bob Makransky
    """
    sun = sphere.set_ecliptical(
        all_planets[0].lon,
        all_planets[0].lat
    )

    assert abs(sun.asc_diff() + 24.7) < 1e-1


def test_moon_umd():
    """
    Compare Moon's upper meridian distance with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        all_planets[1].lon,
        all_planets[1].lat
    )

    assert abs(moon.umd() - 16.1) < 1e-1
