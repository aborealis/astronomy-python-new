"""
Time-related functions
"""
from math import floor
from datetime import datetime, timedelta, timezone


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
    GRIGORIAN_ERA = datetime(1582, 10, 15, tzinfo=timezone.utc)

    # Between Julian and Grigorian calendars
    leapyears = (1582 + 4713 - 1) // 4
    nonleapyears = 1582 + 4713 - 1 - leapyears
    days_btw_cals = (
        leapyears * 366 + nonleapyears * 365
        + (GRIGORIAN_ERA - datetime(1582, 1, 1, tzinfo=timezone.utc)).days
        # 10 days were extracted by Grigory
        # 0.5 days were added to switch from noon-to-noon to
        # midnight-to-midnight day calculation
        - 10 + 0.5
    )

    days_from_grig_era = (
        timeaware_utc - GRIGORIAN_ERA
    ).total_seconds() / 3600 / 24

    return days_btw_cals + days_from_grig_era


def __gst__(timeaware_utc: datetime):
    """
    Returns greenwich sidereal mean time
    """
    jday = __julian_day__(timeaware_utc)
    J2000 = 2451545.0
    midnight = floor(jday) + 0.5
    days_since_midnight = jday - midnight
    hours_since_midnight = days_since_midnight * 24
    days_since_epoch = jday - J2000
    centuries_since_epoch = days_since_epoch / 36525
    whole_days_since_epoch = midnight - J2000

    GMST = (6.697374558
            + 0.06570982441908 * whole_days_since_epoch
            + 1.00273790935 * hours_since_midnight
            + 0.000026 * centuries_since_epoch**2) % 24

    return GMST


def __lst__(timeaware_utc: datetime, geo_lon: float):
    """
    Returns local sidereal time, LST (in hours)
    """
    GST = __gst__(timeaware_utc)
    offset = geo_lon / 360 * 24
    return (GST + offset) % 24


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
