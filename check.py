"""
Compare results with swiss ephemeris module
"""

from datetime import datetime, timedelta
from email.utils import localtime
import components.time as tm
import swisseph as swe


def julian_day(naive_datetime: datetime,
               time_zone: float,) -> float:
    """
    Возвращает юлианский день по дате.
    """
    localtime = tm.__localtime__(naive_datetime, time_zone)
    localtime_utc = tm.__localtime_utc__(localtime)
    return swe.julday(
        localtime_utc.year,
        localtime_utc.month,
        localtime_utc.day,
        localtime_utc.hour + localtime_utc.minute/60
    )


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':
    jday = julian_day(
        datetime(2022, 10, 10, 21, 58),
        time_zone=4,
    )

    ltime = tm.__localtime__(datetime(2022, 10, 10, 21, 58),
                             time_zone=4,)
    ltime_utc = tm.__localtime_utc__(ltime)

    print('Julian day swe:', jday)
    print('Julian day own:', tm.__julian_day__(ltime))
    print('Julian day own:', tm.__julian_day__(ltime_utc))

    print('GST swe:', timedelta(hours=swe.sidtime(jday)))
    print('GST own:', timedelta(hours=tm.__gst__(ltime_utc)))
