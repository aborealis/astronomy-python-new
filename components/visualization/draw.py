"""
Functions for matplotlib visualisation.
"""
from typing import Callable
import numpy as np
from matplotlib.axes import Axes
from sphere import Sphere
from components.vector.coords import Vector
from primary_directions import Directions


def make_data_for_plotting(function: Callable,
                           hrz_angles: np.ndarray,
                           vrt_angles: np.ndarray) -> list[list]:
    """
    Returns three 1d arrays for plotting a parametric
    curve on the celestial sphere.
    """
    _x, _y, _z = [], [], []
    for lat in vrt_angles:
        for lon in hrz_angles:
            vector: Vector = function(lon, lat)
            _x.append(vector.horizontal_xyz().x)
            _y.append(vector.horizontal_xyz().y)
            _z.append(vector.horizontal_xyz().z)
    return [_x, _y, _z]


def surface(sphere: Sphere, axs: Axes):
    """
    Draw celestial surface.
    """
    rasc = np.linspace(0, 360, 37)
    dec = np.linspace(-90, 90, 19)
    rasc, dec = np.meshgrid(rasc, dec)
    xyz_hrz = sphere.set_equatorial(rasc, dec).horizontal_xyz()
    axs.plot_wireframe(*xyz_hrz.aslist(), alpha=0.05)


def horizon(sphere: Sphere, axs: Axes):
    """
    Draws a circle of the horizon.
    """
    azimuth = np.linspace(0, 2 * np.pi, 100)
    _z = 0
    _x = np.sin(azimuth)
    _y = np.cos(azimuth)
    axs.plot(_x, _y, _z, linewidth=0.7, label='Horizon')

    # North - South Axis
    axs.plot([0, 0], [1, -1], color='grey', linewidth=1, linestyle='dashed')
    axs.text(0, -1, 0, s="S")
    axs.text(0, 1, 0, s="N")

    # East - West Axis
    axs.plot([-1, 1], [0, 0], color='grey', linewidth=1, linestyle='dashed')
    axs.text(-1, 0, 0, s="W")
    axs.text(1, 0, 0, s=f"E {int((sphere.ramc + 90) % 360)}º")

    # Polar Axis
    cos_p = sphere.__constants__.cos_p
    sin_p = sphere.__constants__.sin_p
    axs.plot([0, 0], [0, cos_p], [0, sin_p],
             color='grey', linewidth=1, linestyle='dashed')

    # Main Meridian
    alpha = np.linspace(0, np.pi, 100)
    _z = np.sin(alpha)
    _x = np.sin(alpha) * 0
    _y = np.cos(alpha)
    axs.plot(_x, _y, _z, color='grey', linewidth=1,  linestyle='dashed')


def equator(sphere: Sphere, axs: Axes):
    """
    Draws equatorial circle.
    """
    rasc = np.linspace(0, 360, 100)
    dec = 0
    plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz().aslist()
    axs.plot(*plot_data, label='Equator', linewidth=0.7, color="#c4d6e7")

    # Edges of region of never-ascending stars
    rasc = np.linspace(0, 360, 100)
    dec = sphere.__constants__.dec_max
    plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz().aslist()
    axs.plot(*plot_data, linewidth=0.7, color="#c4d6e7")

    rasc = np.linspace(0, 360, 100)
    dec = -sphere.__constants__.dec_max
    plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz().aslist()
    axs.plot(*plot_data, linewidth=0.7, color="#c4d6e7")

    # Equator degrees
    for rasc in range(0, 360, 30):
        xyz = sphere.set_equatorial(rasc, dec=0).horizontal_xyz().aslist()
        axs.text(*xyz, s=str(rasc), fontsize=6, color="cornflowerblue")

    xyz = sphere.set_equatorial(rasc=sphere.ramc, dec=0).horizontal_xyz()
    axs.text(*xyz.aslist(), s="RAMC")


def ecliptic(sphere: Sphere, axs: Axes):
    """
    Draw ecliptic
    """
    # Ecliptic
    lon = np.linspace(0, 360, 100)
    lat = 0
    plot_data = sphere.set_ecliptical(lon, lat).horizontal_xyz().aslist()
    axs.plot(*plot_data, label='Ecliptic', color=(.7, .2, .7))

    # Ecliptic degrees
    for lon in range(0, 360, 30):
        xyz = sphere.set_ecliptical(lon, lat=0).horizontal_xyz().aslist()
        axs.text(*xyz, s=str(lon), fontsize=5, color=(.5, 0, .5))

    # ASC degree
    xyz = sphere.set_ecliptical(lon=sphere.asc, lat=0.0).horizontal_xyz()
    axs.text(*xyz.aslist(), s=f"ASC {int(sphere.asc)}º",
             fontsize=8, color=(.5, 0, .5))


def point(sphere: Sphere, point_data: dict, axs: Axes, line: bool = True):
    """
    Draw a chosen point on celestial sphere.
    """
    vector = sphere.set_ecliptical(
        point_data['lon'],
        point_data['lat']).horizontal_xyz()
    if line:
        axs.plot([vector.x, 0], [vector.y, 0], [vector.z, 0],
                 color='red', linewidth=1, linestyle='dotted')
    axs.scatter(
        vector.x, vector.y, vector.z,
        color=point_data['color'],
        label=point_data['label']
    )


def semiarc(sphere: Sphere, point_data: dict, axs: Axes, s_type: str):
    """
    Draw diurnal semiarc.
    Parameters:
        s_type is DSA or NSA for diurnal/nocturnal semiarc
    """
    vector = sphere.set_ecliptical(point_data['lon'], point_data['lat'])

    if vector.dsa is None:
        return None

    eastern = vector.horizontal_xyz().x > 0
    delta = vector.dsa if s_type == 'DSA' else vector.umd
    if eastern:
        rasc = np.linspace(sphere.ramc, sphere.ramc + delta, 50)
    else:
        rasc = np.linspace(sphere.ramc - delta, sphere.ramc, 50)
    dec = vector.equatorial().dec

    plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz().aslist()
    axs.plot(*plot_data,
             label=s_type,
             linewidth=0.7 if s_type == 'DSA' else 1.2,
             color="red" if s_type == 'DSA' else "blue")
    return None


def eqt_projection(sphere: Sphere, point_data: dict, axs: Axes):
    """
    Draw a projection line from the point to the equator
    """
    vector = sphere.set_ecliptical(
        point_data['lon'],
        point_data['lat']).equatorial()

    rasc = vector.rasc
    dec = np.linspace(0, vector.dec, 50)
    plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz().aslist()
    axs.plot(*plot_data, color="#c4d6e7")


def ecl_projection(sphere: Sphere, point_data: dict, axs: Axes):
    """
    Draw a projection line from the point to the ecliptic
    """
    lon = point_data['lon']
    lat = np.linspace(0, point_data['lat'], 50)
    plot_data = sphere.set_ecliptical(lon, lat).horizontal_xyz().aslist()
    axs.plot(*plot_data, color=(.5, 0, .5))


def placidus_schema(sphere: Sphere, axs: Axes) -> None:
    """
    Illustrates the principle behind the
    Placidus house system.
    """
    ramc = sphere.ramc
    for dec in range(-90, 90):
        dsa = sphere.dsa(dec)
        if dsa is not None:
            rasc = np.linspace(ramc, ramc + dsa/3)
            plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz()
            axs.plot(*plot_data.aslist(), linewidth=0.7, color='lightgray')


def placidus(sphere: Sphere, axs: Axes, under_horizon: bool = False):
    """
    Draw placidus lines on sphere
    """
    ramc = sphere.ramc
    _x = [[], [], [], [], [], ]
    _y = [[], [], [], [], [], ]
    _z = [[], [], [], [], [], ]
    for dec in range(-90, 90):
        dsa = sphere.dsa(dec)
        if dsa is not None:
            for i in range(0, 5):
                if under_horizon:
                    vector = sphere.set_equatorial(
                        ramc + 180 + (i - 2) * dsa/3,
                        -dec).horizontal_xyz
                else:
                    vector = sphere.set_equatorial(
                        ramc + (i - 2) * dsa/3,
                        dec).horizontal_xyz()
                _x[i].append(vector.x)
                _y[i].append(vector.y)
                _z[i].append(vector.z)

    for i in range(0, 5):
        axs.plot(_x[i], _y[i], _z[i], linewidth=0.7, color='green')


def zodiac_2d(sphere: Sphere, points: list[dict], axs: Axes) -> None:
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
        vector = sphere.set_ecliptical(
            lon=_p['lon'],
            lat=_p['lat']).ecliptical_xyz()
        axs.plot([-vector.x, 0], [-vector.y, 0], color='red',
                 linewidth=1, linestyle='dotted')
        axs.scatter(-vector.x, -vector.y, color=_p['color'])


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Primary direction related
# functions:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

def mundane_positions_placidus(sphere: Sphere,
                               acceptor_data: dict,
                               axs: Axes) -> None:
    """
    Draws mundane positions of the acceptor
    """
    # Set acceptor
    directions = Directions(sphere)
    vector = sphere.set_ecliptical(acceptor_data['lon'], acceptor_data['lat'])
    directions.acceptor = vector.equatorial()

    # Draw acceptor point
    point(sphere, acceptor_data, axs)

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
    xyz_hrz = sphere.set_equatorial(conjunction, 0).horizontal_xyz()
    _x0, _y0, _z0 = xyz_hrz.aslist()

    for _m in mundane_positions:
        xyz_hrz = sphere.set_equatorial(_m, 0).horizontal_xyz()
        _x, _y, _z = xyz_hrz.aslist()
        axs.plot([_x, _x0], [_y, _y0], [_z, _z0], color="#c4d6e7",
                 linewidth=1.5, linestyle='solid')
        axs.scatter(
            _x, _y, _z,
            color="#c4d6e7",
            label=None
        )
    return None


def meridian_distance_portions(sphere: Sphere,
                               acceptor_data: dict,
                               axs: Axes) -> None:
    """
    Draws meridian distance portions
    of the acceptor
    """
    # Set acceptor
    directions = Directions(sphere)
    vector = sphere.set_ecliptical(acceptor_data['lon'], acceptor_data['lat'])
    directions.acceptor = vector.equatorial()
    acc_rasc = directions.acceptor.rasc
    acc_dec = directions.acceptor.dec

    # Draw MDPs
    acceptor_quadrant = directions.quadrant(acc_rasc, acc_dec)
    acceptor_mdp = directions.md_portion(acc_rasc, acc_dec)
    if acceptor_mdp is None:
        return None

    if acceptor_quadrant == 0:
        start_point = sphere.ramc
        ratio = acceptor_mdp
    elif acceptor_quadrant == 1:
        start_point = sphere.ramc + 180
        ratio = -1 * acceptor_mdp
    elif acceptor_quadrant == 2:
        start_point = sphere.ramc + 180
        ratio = acceptor_mdp
    else:
        start_point = sphere.ramc
        ratio = -1 * acceptor_mdp

    for dec in range(-90, 90):
        dsa = sphere.dsa(dec)
        if dsa is None:
            continue
        path = dsa if acceptor_quadrant in [0, 3] else 180 - dsa
        rasc = np.linspace(start_point, start_point + path * ratio)
        plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz()
        axs.plot(*plot_data.aslist(), linewidth=0.7, color='lightgray')
    return None


def direction_arc(sphere: Sphere,
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
    directions = Directions(sphere)
    vector = sphere.set_ecliptical(
        promissor_data['lon'], promissor_data['lat'])
    directions.promissor = vector.equatorial()
    vector = sphere.set_ecliptical(
        acceptor_data['lon'], acceptor_data['lat'])
    directions.acceptor = vector.equatorial()

    prom_rasc = directions.promissor.rasc
    prom_dec = directions.promissor.dec

    # Draw acceptor and promissor points
    point(sphere, promissor_data, axs, line=False)

    # Get directional arc
    all_directions = directions.placidus_mundane()
    if all_directions is None:
        return None

    # We'll take a positive arc to exclude reverse direction
    # Alternatively, we may choose the arc with min abs(len)
    arc = sorted([
        item['dist'] for item in all_directions
        if item['aspect'] == aspect
    ], reverse=True)[0]

    # Draw the arc
    rasc = np.linspace(prom_rasc - arc, prom_rasc)
    plot_data = sphere.set_equatorial(rasc, prom_dec).horizontal_xyz()
    axs.plot(*plot_data.aslist(),
             label='direction',
             linewidth=1.2,
             color="red")

    # Draw end point of directional arc on sphere
    xyz_hrz = sphere.set_equatorial(prom_rasc - arc, prom_dec).horizontal_xyz()
    _x, _y, _z = xyz_hrz.aslist()
    axs.plot([_x, 0], [_y, 0], [_z, 0], color='red',
             linewidth=1, linestyle='dotted')
    axs.scatter(
        _x, _y, _z,
        color="grey",
        label=None
    )

    # Draw end point of directional arc on Zodiac
    end_point = sphere.set_equatorial(
        prom_rasc - arc,
        prom_dec
    ).ecliptical()
    end_point_data = dict(
        lon=end_point.lon,
        lat=end_point.lat,
        label='Aspect',
        color="grey",
    )
    print(end_point_data)
    zodiac_2d(sphere, [end_point_data, acceptor_data], axs2)

    meridian_distance_portions(sphere, end_point_data, axs)
    return None
