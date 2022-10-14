"""
Compare results with swiss ephemeris module
"""

from datetime import datetime, timedelta
from email.utils import localtime
import components.time as tm
import swisseph as swe
from vector import Vector


def julian_day(aware_datetime: datetime) -> float:
    """
    Возвращает юлианский день по дате.
    """
    return swe.julday(
        aware_datetime.year,
        aware_datetime.month,
        aware_datetime.day,
        aware_datetime.hour + aware_datetime.minute/60
    )


def get_asc(jday: float, geo_lon: float, geo_lat: float) -> list[float]:
    """
    Извлекает информацию по куспидам домов.
    """
    cusp_degrees = swe.houses(
        jday, geo_lat, geo_lon, bytes('R', 'ascii')
    )[0]
    return cusp_degrees[0]


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':
    dt = datetime(2022, 10, 10, 17, 40)
    geo = dict(
        geo_lon=44 + 46/60,
        geo_lat=25 + 41.713664,
    )

    ltime = tm.__localtime__(dt,
                             time_zone=4,)
    ltime_utc = tm.__localtime_utc__(ltime)

    jday = julian_day(ltime_utc)
    print('Julian day swe:', jday)
    print('Julian day own:', tm.__julian_day__(ltime))
    print('Julian day own:', tm.__julian_day__(ltime_utc))

    print('GST swe:', timedelta(hours=swe.sidtime(jday)))
    print('GST own:', timedelta(hours=tm.__gst__(ltime_utc)))

    swe_asc = get_asc(
        jday, **geo)
    vector = Vector(
        naive_datetime=dt,
        time_zone=4,
        **geo,
    )
    own_asc = vector.asc()

    print('ASC swe:', swe_asc)
    print('ASC_own:', own_asc)

    vector.set_ecliptical(swe_asc, 0)
    print('altitude (swe asc):', vector.horizontal().alt)
    vector.set_ecliptical(own_asc, 0)
    print('altitude (own asc):', vector.horizontal().alt)
