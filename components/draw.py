import numpy as np
from matplotlib.axes import Axes
from vector import Vector


def surface(vector: Vector, ax: Axes):
    """
    Draw celestial surface with horizon
    """
    ra = np.linspace(0, 360, 30)
    dec = np.linspace(0, 360, 30)
    ra, dec = np.meshgrid(ra, dec)

    vector.set_equatorial(ra, dec)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.plot_wireframe(x, y, z, alpha=0.05)

    # Horizon
    azimuth = np.linspace(0, 2 * np.pi, 100)
    z = 0
    x = np.sin(azimuth)
    y = np.cos(azimuth)
    ax.plot(x, y, z, label='Horizon')

    # Equator
    ra = np.linspace(0, 360, 100)
    vector.set_equatorial(ra, dec=0)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.plot(x, y, z, label='Equator', linewidth=0.7, color="#c4d6e7")

    # Equator degrees
    for ra in range(0, 360, 30):
        vector.set_equatorial(ra, dec=0)
        x = vector.cartesian_horizontal().x
        y = vector.cartesian_horizontal().y
        z = vector.cartesian_horizontal().z
        ax.text(x, y, z, s=str(ra), fontsize=6, color="cornflowerblue")

    # North - South Axis
    ax.plot([0, 0], [1, -1], color='grey', linewidth=1, linestyle='dashed')
    ax.text(0, -1, 0, s="S")
    ax.text(0, 1, 0, s="N")

    # East - West Axis
    ax.plot([-1, 1], [0, 0], color='grey', linewidth=1, linestyle='dashed')
    ax.text(-1, 0, 0, s="W")
    ax.text(1, 0, 0, s=f"E {int((vector.ramc() + 90) % 360)}º")

    # Polar Axis
    cos_p = vector.__constants__.cos_p
    sin_p = vector.__constants__.sin_p
    ax.plot([0, 0], [0, cos_p], [0, sin_p],
            color='grey', linewidth=1, linestyle='dashed')

    # Main Meridian
    alpha = np.linspace(0, np.pi, 100)
    z = np.sin(alpha)
    x = np.sin(alpha) * 0
    y = np.cos(alpha)
    ax.plot(x, y, z, color='grey', linewidth=1,  linestyle='dashed')

    vector.set_equatorial(ra=vector.ramc(), dec=0)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.text(x, y, z, s="RAMC")


def ecliptic(vector: Vector, ax: Axes):
    """
    Draw ecliptic
    """
    # Ecliptic
    lon = np.linspace(0, 360, 100)
    vector.set_ecliptical(lon, lat=0)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.plot(x, y, z, label='Ecliptic', color=(.7, .2, .7))

    # Ecliptic degrees
    for lon in range(0, 360, 30):
        vector.set_ecliptical(lon, lat=0)
        x = vector.cartesian_horizontal().x
        y = vector.cartesian_horizontal().y
        z = vector.cartesian_horizontal().z
        ax.text(x, y, z, s=str(lon), fontsize=5, color=(.5, 0, .5))

    # ASC degree
    vector.set_ecliptical(lon=vector.asc(), lat=0.0)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.text(x, y, z, s=f"ASC {int(vector.asc())}º",
            fontsize=8, color=(.5, 0, .5))


def point(vector: Vector, vector_coordinates: dict, ax: Axes):
    """
    Draw chosen point
    """
    # Draw pont
    vector.set_ecliptical(**vector_coordinates)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.text(x, y, z, s='P', color="red")
    ax.plot([x, 0], [y, 0], [z, 0], color='red',
            linewidth=1, linestyle='dotted')


def semiarc(vector: Vector, vector_coordinates: dict, ax: Axes, s_type: str):
    """
    Draw diurnal semiarc
    """
    vector.set_ecliptical(**vector_coordinates)
    # Point parameters
    point_ra = vector.equatorial().ra
    point_dec = vector.equatorial().dec
    point_dsa = vector.dsa()
    point_umd = vector.umd()
    eastern = vector.cartesian_horizontal().x > 0

    # Draw DSA
    delta = point_dsa if s_type == 'DSA' else point_umd
    if eastern:
        ra = np.linspace(vector.ramc(), vector.ramc() + delta, 100)
    else:
        ra = np.linspace(vector.ramc() - delta, vector.ramc(), 100)
    vector.set_equatorial(ra, dec=point_dec)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.plot(x, y, z,
            label=s_type,
            linewidth=0.7 if s_type == 'DSA' else 1.2,
            color="red" if s_type == 'DSA' else "blue")
