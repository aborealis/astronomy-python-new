"""
Functions for matplotlib visualisation of stellar
objects on celestial sphere.
"""
from typing import Callable
import numpy as np
from matplotlib.axes import Axes
from sphere import Sphere
from components.vector.coords import Vector


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

    if vector.dsa() is None:
        return None

    eastern = vector.horizontal_xyz().x > 0
    delta = vector.dsa() if s_type == 'DSA' else vector.umd()
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
