"""
Functions for matplotlib visualisation
"""
import numpy as np
from matplotlib.axes import Axes
from vector import Vector
from primary_directions import Directions


def xyz_hrz(vector: Vector) -> list:
    """
    Return x, y, z of a given vector in
    horizontal system
    """
    return [
        vector.cartesian_horizontal().x,
        vector.cartesian_horizontal().y,
        vector.cartesian_horizontal().z,
    ]


def surface(vector: Vector, axs: Axes):
    """
    Draw celestial surface with horizon
    """
    rasc = np.linspace(0, 360, 30)
    dec = np.linspace(0, 360, 30)
    rasc, dec = np.meshgrid(rasc, dec)
    vector.set_equatorial(rasc, dec)
    axs.plot_wireframe(*xyz_hrz(vector), alpha=0.05)

    # Horizon
    azimuth = np.linspace(0, 2 * np.pi, 100)
    _z = 0
    _x = np.sin(azimuth)
    _y = np.cos(azimuth)
    axs.plot(_x, _y, _z, label='Horizon')

    # Equator
    rasc = np.linspace(0, 360, 100)
    vector.set_equatorial(rasc, dec=0)
    axs.plot(*xyz_hrz(vector), label='Equator', linewidth=0.7, color="#c4d6e7")

    # Edges of region of never-ascending stars
    rasc = np.linspace(0, 360, 100)
    vector.set_equatorial(rasc, dec=vector.__constants__.dec_max)
    axs.plot(*xyz_hrz(vector), linewidth=0.7, color="#c4d6e7")

    rasc = np.linspace(0, 360, 100)
    vector.set_equatorial(rasc, dec=-vector.__constants__.dec_max)
    axs.plot(*xyz_hrz(vector), linewidth=0.7, color="#c4d6e7")

    # Equator degrees
    for rasc in range(0, 360, 30):
        vector.set_equatorial(rasc, dec=0)
        axs.text(*xyz_hrz(vector), s=str(rasc),
                 fontsize=6, color="cornflowerblue")

    # North - South Axis
    axs.plot([0, 0], [1, -1], color='grey', linewidth=1, linestyle='dashed')
    axs.text(0, -1, 0, s="S")
    axs.text(0, 1, 0, s="N")

    # East - West Axis
    axs.plot([-1, 1], [0, 0], color='grey', linewidth=1, linestyle='dashed')
    axs.text(-1, 0, 0, s="W")
    axs.text(1, 0, 0, s=f"E {int((vector.ramc() + 90) % 360)}º")

    # Polar Axis
    cos_p = vector.__constants__.cos_p
    sin_p = vector.__constants__.sin_p
    axs.plot([0, 0], [0, cos_p], [0, sin_p],
             color='grey', linewidth=1, linestyle='dashed')

    # Main Meridian
    alpha = np.linspace(0, np.pi, 100)
    _z = np.sin(alpha)
    _x = np.sin(alpha) * 0
    _y = np.cos(alpha)
    axs.plot(_x, _y, _z, color='grey', linewidth=1,  linestyle='dashed')

    vector.set_equatorial(rasc=vector.ramc(), dec=0)
    axs.text(*xyz_hrz(vector), s="RAMC")


def ecliptic(vector: Vector, axs: Axes):
    """
    Draw ecliptic
    """
    # Ecliptic
    lon = np.linspace(0, 360, 100)
    vector.set_ecliptical(lon, lat=0)
    axs.plot(*xyz_hrz(vector), label='Ecliptic', color=(.7, .2, .7))

    # Ecliptic degrees
    for lon in range(0, 360, 30):
        vector.set_ecliptical(lon, lat=0)
        axs.text(*xyz_hrz(vector), s=str(lon), fontsize=5, color=(.5, 0, .5))

    # ASC degree
    vector.set_ecliptical(lon=vector.asc(), lat=0.0)
    axs.text(*xyz_hrz(vector), s=f"ASC {int(vector.asc())}º",
             fontsize=8, color=(.5, 0, .5))


def point(vector: Vector, point_data: dict, axs: Axes, line: bool = True):
    """
    Draw chosen point
    """
    # Draw pont
    vector.set_ecliptical(point_data['lon'], point_data['lat'])
    _x, _y, _z = xyz_hrz(vector)
    if line:
        axs.plot([_x, 0], [_y, 0], [_z, 0], color='red',
                 linewidth=1, linestyle='dotted')
    axs.scatter(
        _x, _y, _z,
        color=point_data['color'],
        label=point_data['label']
    )


def semiarc(vector: Vector, point_data: dict, axs: Axes, s_type: str):
    """
    Draw diurnal semiarc
    s_type is DSA or NSA for diurnal/nocturnal semiarc
    """
    vector.set_ecliptical(point_data['lon'], point_data['lat'])

    # Point parameters
    point_dec = vector.equatorial().dec
    point_dsa = vector.dsa()
    if point_dsa is None:
        return None
    point_umd = vector.umd()
    eastern = vector.cartesian_horizontal().x > 0

    # Draw DSA
    delta = point_dsa if s_type == 'DSA' else point_umd
    if eastern:
        rasc = np.linspace(vector.ramc(), vector.ramc() + delta, 50)
    else:
        rasc = np.linspace(vector.ramc() - delta, vector.ramc(), 50)
    vector.set_equatorial(rasc, dec=point_dec)
    axs.plot(*xyz_hrz(vector),
             label=s_type,
             linewidth=0.7 if s_type == 'DSA' else 1.2,
             color="red" if s_type == 'DSA' else "blue")
    return None


def eqt_projection(vector: Vector, point_data: dict, axs: Axes):
    """
    Draw a projection line from the point to the equator
    """
    vector.set_ecliptical(point_data['lon'], point_data['lat'])

    # Point parameters
    point_ra = vector.equatorial().rasc
    point_dec = vector.equatorial().dec
    dec = np.linspace(0, point_dec, 50)
    vector.set_equatorial(point_ra, dec)
    axs.plot(*xyz_hrz(vector), color="#c4d6e7")


def ecl_projection(vector: Vector, point_data: dict, axs: Axes):
    """
    Draw a projection line from the point to the ecliptic
    """
    vector.set_ecliptical(point_data['lon'], point_data['lat'])

    # Point parameters
    point_lon = vector.ecliptical().lon
    point_lat = vector.ecliptical().lat
    lat = np.linspace(0, point_lat, 50)
    vector.set_ecliptical(point_lon, lat)
    axs.plot(*xyz_hrz(vector), color=(.5, 0, .5))


def placidus_schema(vector: Vector, axs: Axes) -> None:
    """
    Illustrates the principle behind the
    Placidus house system
    """
    ramc = vector.ramc()
    for dec in range(-90, 90):
        dsa = vector.dsa(dec)
        if dsa is not None:
            rasc = np.linspace(ramc, ramc + dsa/3)
            vector.set_equatorial(rasc, dec)
            axs.plot(*xyz_hrz(vector), linewidth=0.7, color='lightgray')


def placidus(vector: Vector, axs: Axes, under_horizon: bool = False):
    """
    Draw placidus lines on sphere
    """
    ramc = vector.ramc()
    _x = [[], [], [], [], [], ]
    _y = [[], [], [], [], [], ]
    _z = [[], [], [], [], [], ]
    for dec in range(-90, 90):
        vector.set_equatorial(0, dec)
        dsa = vector.dsa()
        if dsa is not None:
            for i in range(0, 5):
                if under_horizon:
                    vector.set_equatorial(ramc + 180 + (i - 2) * dsa/3, -dec)
                else:
                    vector.set_equatorial(ramc + (i - 2) * dsa/3, dec)
                _xx, _yy, _zz = xyz_hrz(vector)
                _x[i].append(_xx)
                _y[i].append(_yy)
                _z[i].append(_zz)

    for i in range(0, 5):
        axs.plot(_x[i], _y[i], _z[i], linewidth=0.7, color='green')


def zodiac_2d(vector: Vector, points: list[dict], axs: Axes) -> None:
    """
    Draws 2d Zodiac circle
    """
    alpha = np.linspace(0, 2 * np.pi)
    _x = np.cos(alpha)
    _y = np.sin(alpha)
    axs.plot(_x, _y, color=(.7, .2, .7))
    axs.text(-0.95, -0.1, s='0', fontsize=6)
    axs.text(-0.1, -0.9, s='90', fontsize=6)
    for _p in points:
        vector.set_ecliptical(lon=_p['lon'], lat=_p['lat'])
        xyz = vector.cartesian()
        axs.plot([-xyz.x, 0], [-xyz.y, 0], color='red',
                 linewidth=1, linestyle='dotted')
        axs.scatter(-xyz.x, -xyz.y, color=_p['color'])


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Primary direction related
# functions:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

def mundane_positions_placidus(vector: Vector,
                               acceptor_data: dict,
                               axs: Axes) -> None:
    """
    Draws mundane positions of the acceptor
    """
    # Set acceptor
    directions = Directions(vector)
    vector.set_ecliptical(acceptor_data['lon'], acceptor_data['lat'])
    directions.acceptor = vector.equatorial()

    # Draw acceptor point
    point(vector, acceptor_data, axs)

    # Find mundane positions of the acceptor
    mundane_positions = directions.mundane_positions_placidus()
    if not mundane_positions:
        return None
    conjunction = [
        item['rasc'] for item in mundane_positions
        if item['aspect'] == 0
    ][0]
    mundane_positions = [
        item['rasc'] for item in mundane_positions
    ]

    # Draw acceptor's mundane position
    vector.set_equatorial(conjunction, 0)
    _x0, _y0, _z0 = xyz_hrz(vector)

    for _m in mundane_positions:
        vector.set_equatorial(_m, 0)
        _x, _y, _z = xyz_hrz(vector)
        axs.plot([_x, _x0], [_y, _y0], [_z, _z0], color="#c4d6e7",
                 linewidth=1.5, linestyle='solid')
        axs.scatter(
            _x, _y, _z,
            color="#c4d6e7",
            label=None
        )
    return None


def meridian_distance_portions(vector: Vector,
                               acceptor_data: dict,
                               axs: Axes) -> None:
    """
    Draws meridian distance portions
    of the acceptor
    """
    # Set acceptor
    directions = Directions(vector)
    vector.set_ecliptical(acceptor_data['lon'], acceptor_data['lat'])
    directions.acceptor = vector.equatorial()
    acc_rasc = directions.acceptor.rasc
    acc_dec = directions.acceptor.dec

    # Draw MDPs
    acceptor_quadrant = directions.quadrant(acc_rasc, acc_dec)
    acceptor_mdp = directions.md_portion(acc_rasc, acc_dec)
    if acceptor_mdp is None:
        return None

    if acceptor_quadrant == 0:
        start_point = vector.ramc()
        ratio = acceptor_mdp
    elif acceptor_quadrant == 1:
        start_point = vector.ramc() + 180
        ratio = -1 * acceptor_mdp
    elif acceptor_quadrant == 2:
        start_point = vector.ramc() + 180
        ratio = acceptor_mdp
    else:
        start_point = vector.ramc()
        ratio = -1 * acceptor_mdp

    for dec in range(-90, 90):
        dsa = vector.dsa(dec)
        if dsa is None:
            continue
        path = dsa if acceptor_quadrant in [0, 3] else 180 - dsa
        rasc = np.linspace(start_point, start_point + path * ratio)
        vector.set_equatorial(rasc, dec)
        axs.plot(*xyz_hrz(vector), linewidth=0.7, color='lightgray')
    return None


def direction_arc(vector: Vector,
                  promissor_data: dict,
                  acceptor_data: dict,
                  aspect: int,
                  axs: Axes,
                  axs2: Axes) -> None:
    """
    Draws directional arc between promissor
    and acceptor's aspect
    """
    # Set promissor and significator
    directions = Directions(vector)
    vector.set_ecliptical(promissor_data['lon'], promissor_data['lat'])
    directions.promissor = vector.equatorial()
    vector.set_ecliptical(acceptor_data['lon'], acceptor_data['lat'])
    directions.acceptor = vector.equatorial()

    prom_rasc = directions.promissor.rasc
    prom_dec = directions.promissor.dec

    # Draw acceptor and promissor points
    point(vector, promissor_data, axs, line=False)

    # Get directional arc
    all_directions = directions.direction_placidus()
    if all_directions is None:
        return None

    # We'll take a positive arc to exclude reverse direction
    # Alternatively, we may choose the arc with min abs(len)
    arc = sorted([
        item['dist'] for item in all_directions
        if item['aspect'] == aspect
    ], reverse=True)[0]
    print(arc)

    # Draw the arc
    rasc = np.linspace(prom_rasc - arc, prom_rasc)
    vector.set_equatorial(rasc, prom_dec)
    axs.plot(*xyz_hrz(vector),
             label='direction',
             linewidth=1.2,
             color="red")

    # Draw end point of directional arc on sphere
    vector.set_equatorial(prom_rasc - arc, prom_dec)
    _x, _y, _z = xyz_hrz(vector)
    axs.plot([_x, 0], [_y, 0], [_z, 0], color='red',
             linewidth=1, linestyle='dotted')
    axs.scatter(
        _x, _y, _z,
        color="grey",
        label=None
    )

    # Draw end point of directional arc on Zodiac
    end_point = vector.ecliptical()
    end_point_data = dict(
        lon=end_point.lon,
        lat=end_point.lat,
        label='Aspect',
        color="red",
    )
    zodiac_2d(vector, [end_point_data, acceptor_data], axs2)

    meridian_distance_portions(vector, end_point_data, axs)
    return None
