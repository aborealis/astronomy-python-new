"""
Compare primary directions with the
results of Bob Makransky
"""
from datetime import datetime
from typing import NamedTuple
import components.vector.time as tm
from sphere import Sphere
from primary_directions import Directions
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
all_houses = swe.scan_houses(
    localtime_utc,
    data.geo_lon,
    data.geo_lat,
    system='P'
)
sphere = Sphere(data.naive_datetime,
                data.time_zone,
                data.geo_lon,
                data.geo_lat)
directions = Directions(sphere)


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Assertions:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

def test_placidus_mundi_moon_conj_mc():
    """
    Compare Moon-MC mundane conjunction with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        all_planets[1].lon,
        all_planets[1].lat
    )
    _mc = sphere.set_ecliptical(
        sphere.medium_coeli,
        0
    )

    all_directions = directions.placidus_mundane(moon, _mc, aspect=0)
    conj = all_directions[0]['dist']

    assert abs(conj - 16.1) < 1e-1


def test_placidus_zodiac_moon_conj_mc():
    """
    Compare Moon-MC zodiacal conjunction with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        all_planets[1].lon,
        0
    )
    _mc = sphere.set_ecliptical(
        sphere.medium_coeli,
        0
    )

    all_directions = directions.placidus_mundane(moon, _mc, aspect=0)
    conj = all_directions[0]['dist']

    assert abs(conj - 15.95) < 1e-1


def test_placidus_mundi_jup_conj_dsc():
    """
    Compare Jupiter-DSC mundane conjunction with result of Bob Makransky
    """
    jupiter = sphere.set_ecliptical(
        all_planets[5].lon,
        all_planets[5].lat
    )
    dsc = sphere.set_ecliptical(
        (sphere.asc + 180) % 360,
        0
    )

    all_directions = directions.placidus_mundane(jupiter, dsc, aspect=0)
    conj = all_directions[0]['dist']
    assert abs(conj + 45.43) < 1e-1


def test_placidus_zodiac_jup_conj_dsc():
    """
    Compare Jupiter-DSC zodiacal conjunction with result of Bob Makransky
    """
    jupiter = sphere.set_ecliptical(
        all_planets[5].lon,
        0
    )
    dsc = sphere.set_ecliptical(
        (sphere.asc + 180) % 360,
        0
    )

    all_directions = directions.placidus_mundane(jupiter, dsc, aspect=0)
    conj = all_directions[0]['dist']
    assert abs(conj + 45.54) < 1e-1


def test_placidus_mundi_sun_conj_mer():
    """
    Compare Sun-Merc mundane conjunction with result of Bob Makransky
    """
    sun = sphere.set_ecliptical(
        all_planets[0].lon,
        all_planets[0].lat
    )
    mer = sphere.set_ecliptical(
        all_planets[2].lon,
        all_planets[2].lat
    )

    all_directions = directions.placidus_mundane(sun, mer, aspect=0)
    conj = all_directions[0]['dist']
    assert abs(conj - 12.85) < 1e-1


def test_placidus_mundi_mon_tri_sat():
    """
    Compare Moon-Saturn mundane trigone with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        all_planets[1].lon,
        all_planets[1].lat
    )
    sat = sphere.set_ecliptical(
        all_planets[6].lon,
        all_planets[6].lat
    )

    all_directions = directions.placidus_mundane(moon, sat, aspect=120)
    all_directions.sort(key=lambda item: abs(item['dist']))
    shortest = all_directions[0]['dist']
    assert abs(shortest + 5.77) < 1e-1


def test_placidus_zodiac_mon_tri_sat():
    """
    Compare Moon-Saturn zodiacal trigone with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        (all_planets[1].lon + 120) % 360,
        0
    )
    sat = sphere.set_ecliptical(
        all_planets[6].lon,
        all_planets[6].lat
    )

    all_directions = directions.placidus_mundane(moon, sat)
    aspect = [
        item['dist']
        for item in all_directions
        if item['aspect'] == 0
    ][0]

    assert abs(aspect + 5.38) < 1e-1


def test_placidus_zodiac_mon_tri_sat2():
    """
    Compare Moon-Saturn zodiacal trigone with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        all_planets[1].lon,
        0
    )
    sat = sphere.set_ecliptical(
        all_planets[6].lon,
        all_planets[6].lat
    )

    all_directions = directions.placidus_zodiac(moon, sat, 120)
    all_directions.sort(key=lambda item: abs(item['dist']))
    shortest = all_directions[0]['dist']
    assert abs(shortest + 5.38) < 1e-1


def test_placidus_field_plane_mon_sqr_asc():
    """
    Compare Moon-ASC zodiacal quadrature on field plane
    latitude 5 N 11 with result of Bob Makransky
    """
    moon = sphere.set_ecliptical(
        all_planets[1].lon,
        0
    )
    asc = sphere.set_equatorial(
        sphere.ramc + 90,
        0
    )

    all_directions = directions.placidus_zodiac(
        moon, asc, 90, field_plane_lat=5 + 11/60)
    all_directions.sort(key=lambda item: abs(item['dist']))
    shortest = all_directions[0]['dist']
    assert abs(shortest + 14.63) < 1e-1


def test_regio_mundi_sun_conj_mer():
    """
    Compare Sun-Mercury mundane conj in Regiomontanus system
    with result of Bob Makransky
    """
    sun = sphere.set_ecliptical(
        all_planets[0].lon,
        0
    )
    mer = sphere.set_ecliptical(
        all_planets[2].lon,
        all_planets[2].lat
    )
    arc = directions.regio_mundane(sun, mer, 0)[0]['dist']
    assert abs(arc - 12.17) < 1e-1


def test_regio_mundi_sat_conj_ven():
    """
    Compare Sat-Venus mundane conj in Regiomontanus system
    with result of Bob Makransky
    """
    sat = sphere.set_ecliptical(
        all_planets[6].lon,
        all_planets[6].lat
    )
    ven = sphere.set_ecliptical(
        all_planets[3].lon,
        all_planets[3].lat
    )
    arc = directions.regio_mundane(sat, ven, 0)[0]['dist']
    assert abs(arc + 37.09) < 1e-1


def test_regio_mundi_mon_tri_sat():
    """
    Compare Sat-Venus mundane conj in Regiomontanus system
    with result of Bob Makransky
    """
    sat = sphere.set_ecliptical(
        all_planets[6].lon,
        all_planets[6].lat
    )
    mon = sphere.set_ecliptical(
        all_planets[1].lon,
        all_planets[1].lat
    )
    all_directions = directions.regio_mundane(mon, sat, 120)
    all_directions.sort(key=lambda item: abs(item['dist']))
    shortest = all_directions[0]['dist']
    assert abs(shortest + 3.19) < 1e-1
