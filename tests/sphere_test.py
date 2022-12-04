"""
Test for astronomical calculations for
King Charles' birt data
"""

from datetime import datetime
import swisseph as swe
import components.vector.time as tm
from sphere import Sphere
from components.vector.planets_and_stars import StellarObject


test_sphere = Sphere(
    datetime(1948, 11, 14, 21, 14, 39),
    time_zone=0,
    geo_lon=-1/6,
    geo_lat=51.5
)

stellar_objects = StellarObject(test_sphere)


def scan_houses(sphere: Sphere, system: str) -> list[float]:
    """
    Returns a list of house cusp in a given system
    """
    return swe.houses(
        sphere.jday, sphere.geo_lat, sphere.geo_lon, bytes(system, 'ascii')
    )[0]


all_planets = stellar_objects.planets
all_houses_placidus = scan_houses(test_sphere, system='P')
all_houses_regio = scan_houses(test_sphere, system='R')


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Assertions:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

def test_julian_time_is_correct():
    """
    Compare self-programmed jt with one from the swisseph module
    """
    ltime_utc = test_sphere.localtime_utc
    assert tm.__julian_day__(ltime_utc) == swe.julday(
        ltime_utc.year,
        ltime_utc.month,
        ltime_utc.day,
        (ltime_utc.hour
         + ltime_utc.minute/60
         + ltime_utc.second/3600)
    )


def test_asc_is_on_horizon():
    """
    ASC should always be on the horizon
    """
    asc = test_sphere.set_ecliptical(test_sphere.asc, 0)
    assert abs(asc.horizontal().alt) < 1e-10


def test_mc_is_on_prime_meridian():
    """
    MC should always be on the prime meridian
    """
    _mc = test_sphere.set_ecliptical(test_sphere.medium_coeli, 0)
    assert abs(_mc.horizontal_xyz().x) < 1e-10
    assert _mc.horizontal_xyz().z >= 0


def test_placidus_5th_house_cusp():
    """
    Compare my calculations with swisseph module
    for placidus house system
    """
    diff = abs(test_sphere.placidus(5)[0] - all_houses_placidus[4])
    assert diff < 1e-2


def test_regiomontanus_5th_house_cusp():
    """
    Compare my calculations with swisseph module
    for placidus house system
    """
    cusp5 = test_sphere.regiomontanus(5)
    cusp5swe = all_houses_regio[4]
    diff = abs(cusp5 - cusp5swe)
    assert diff < 1e-2


def test_placidus_5th_lmd_is_third_of_nsa():
    """
    lower meridian distance for 5th house should be
    one third of nocturnal semiarc
    """
    cusp5 = test_sphere.set_ecliptical(test_sphere.placidus(5)[0], 0)
    nsa = 180 - cusp5.dsa()
    lmd = 180 - cusp5.umd()
    assert abs(lmd - nsa / 3) < 1e-9


def test_pluto_equatorial_coordinates():
    """
    Compare Pluto's coordinates with results of Bob Makransky
    """
    pluto = test_sphere.set_ecliptical(
        all_planets[9].lon,
        all_planets[9].lat
    ).equatorial()

    assert abs(pluto.rasc - 141.49) < 1e-1
    assert abs(pluto.dec - 23 - 4/60) < 1e-2


def test_sun_asc_diff():
    """
    Compare Pluto's asc difference with result of Bob Makransky
    """
    sun = test_sphere.set_ecliptical(
        all_planets[0].lon,
        all_planets[0].lat
    )

    assert abs(sun.asc_diff() + 24.7) < 1e-1


def test_moon_umd():
    """
    Compare Moon's upper meridian distance with result of Bob Makransky
    """
    moon = test_sphere.set_ecliptical(
        all_planets[1].lon,
        all_planets[1].lat
    )

    assert abs(moon.umd() - 16.1) < 1e-1
