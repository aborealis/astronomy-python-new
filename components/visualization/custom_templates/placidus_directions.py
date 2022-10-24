"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw


def points(sphere: Sphere) -> dict[dict]:
    """
    Returns a dictionary of points to display on
    celestial_sphere and on Zodiac circle
    """

    acceptor = dict(
        lon=140,
        lat=0,
        label='Acceptor',
        color="blue",
    )

    promissor = dict(
        lon=150,
        lat=0,
        label='Promissor',
        color="red",
    )

    point_mc = dict(
        lon=sphere.medium_coeli,
        lat=0,
        label="MC",
        color="orange",
    )

    return dict(
        acceptor=acceptor,
        promissor=promissor,
        mc=point_mc,
        label="MC"
    )


def figures(sphere: Sphere, figure3d: Axes, figure2d: Axes) -> None:
    """
    Call draw functions to include points
    and objects to visualization on celestial
    sphere and Zodiac circle
    """

    pnt = points(sphere)
    draw.surface(sphere, figure3d)
    draw.ecliptic(sphere, figure3d)
    draw.horizon(sphere, figure3d)
    draw.equator(sphere, figure3d)
    draw.mundane_positions_placidus(sphere, pnt['acceptor'], figure3d)
    draw.meridian_distance_portions(sphere, pnt['acceptor'], figure3d)
    draw.direction_arc(
        sphere,
        pnt['promissor'],
        pnt['acceptor'],
        60, figure3d, figure2d)
    draw.point(sphere, pnt['mc'], figure3d, line=False)
