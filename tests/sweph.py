"""
Wrappers for swisseph module's functions
"""
from typing import NamedTuple
from datetime import datetime
import swisseph as swe


class Planet(NamedTuple):
    """
    Planetary coordinates
    """
    lon: float
    lat: float


def julian_day(ltime_utc: datetime) -> float:
    """
    Возвращает юлианский день по дате.
    """
    return swe.julday(
        ltime_utc.year,
        ltime_utc.month,
        ltime_utc.day,
        (ltime_utc.hour
         + ltime_utc.minute/60
         + ltime_utc.second/3600)
    )


def scan_planets(ltime_utc: datetime) -> list[Planet]:
    """
    Returns a list of 9 visible planets [Sun...Pluto]
    """
    result = []
    for planet in range(10):
        raw_data = swe.calc_ut(julian_day(ltime_utc), planet)
        result.append(Planet(
            lon=float(raw_data[0][0]),
            lat=float(raw_data[0][1])
        ))
    return result


def scan_houses(ltime_utc: datetime,
                geo_lon: float,
                geo_lat: float,
                system: str) -> list[float]:
    """
    Returns a list of house cusp in a given system
    """
    return swe.houses(
        julian_day(ltime_utc), geo_lat, geo_lon, bytes(system, 'ascii')
    )[0]
