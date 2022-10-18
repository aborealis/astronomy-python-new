"""
Time-related functions
"""
from datetime import datetime, timedelta, timezone
from math import floor
from numpy import sin, cos, pi


def __localtime__(naive_datetime: datetime,
                  time_zone: float) -> datetime:
    """
    Creates time-aware localtime object
    """
    offset = timedelta(hours=time_zone)
    return datetime(naive_datetime.year, naive_datetime.month,
                    naive_datetime.day, naive_datetime.hour,
                    naive_datetime.minute, naive_datetime.second,
                    naive_datetime.microsecond,
                    tzinfo=timezone(offset))


def __localtime_utc__(localtime: datetime):
    """
    Express time aware localtime object in
    universal time
    """
    return localtime.astimezone(timezone.utc)


def __julian_day__(timeaware_utc: datetime) -> float:
    """
    Returns Julian day for a given local time
    """
    grigorian_era = datetime(1582, 10, 15, tzinfo=timezone.utc)

    # Between Julian and Grigorian calendars
    leapyears = (1582 + 4713 - 1) // 4
    nonleapyears = 1582 + 4713 - 1 - leapyears
    days_btw_cals = (
        leapyears * 366 + nonleapyears * 365
        + (grigorian_era - datetime(1582, 1, 1, tzinfo=timezone.utc)).days
        # 10 days were extracted by Grigory
        # 0.5 days were added to switch from noon-to-noon to
        # midnight-to-midnight day calculation
        - 10 + 0.5
    )

    days_from_grig_era = (
        timeaware_utc - grigorian_era
    ).total_seconds() / 3600 / 24

    return days_btw_cals + days_from_grig_era


def __gst__(timeaware_utc: datetime):
    """
    Returns greenwich sidereal apparent time
    """
    jday = __julian_day__(timeaware_utc)

    # 2000 January 1, 12h
    j2000 = 2451545.0

    if jday - floor(jday) < 0.5:
        # from noon to midnight
        past_midnight = floor(jday) - 0.5
    else:
        # from midnight to noon
        past_midnight = floor(jday) + 0.5

    days_since_j2000 = jday - j2000
    centuries_since_j2000 = days_since_j2000 / 36525
    delta_jdays_ut = past_midnight - j2000
    hours_since_midnight = (jday - past_midnight) * 24
    gmst = (6.697375
            + 0.065707485828 * delta_jdays_ut
            + 1.0027379 * hours_since_midnight
            + 0.0854103 * centuries_since_j2000
            + 0.0000258 * centuries_since_j2000**2) % 24

    # The equation of the equinoxes
    epsilon = __inclination_ecliptic__(timeaware_utc)
    sun_mean_long = (280.47 + 0.98565 * days_since_j2000) % 360
    nnode_long = (125.04 - 0.052954 * days_since_j2000) % 360
    eqeq = (-0.000319 * sin(nnode_long * pi / 180)
            - 0.000024 * sin(2 * sun_mean_long / 180)) * cos(epsilon * pi / 180)

    return gmst + eqeq


def __lst__(timeaware_utc: datetime, geo_lon: float):
    """
    Returns local sidereal time, LST (in hours)
    """
    gst = __gst__(timeaware_utc)
    offset = geo_lon / 360 * 24
    return (gst + offset) % 24


def __inclination_ecliptic__(timeaware_utc: datetime) -> float:
    """
    Returns the angle of ecliptic inclanation
    in degrees
    """
    jday = __julian_day__(timeaware_utc)
    centuries_since_epoch = (jday - 2415020.0) / 36525
    return (
        23 + 27 / 60 + 8.26 / 3600 -
        46.845 / 3600 * centuries_since_epoch +
        0.0059 / 3600 * centuries_since_epoch ** 2 +
        0.001811 / 3600 * centuries_since_epoch ** 3
    )
