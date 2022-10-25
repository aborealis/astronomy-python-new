"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw

ASPECT = 60


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

    promissor_aspect = dict(
        lon=150 + ASPECT,
        lat=0,
        label=None,
        color=(.7, .2, .7)
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
        promissor_aspect=promissor_aspect,
        mc=point_mc,
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
    draw.zodiac_positions_placidus(sphere, pnt['promissor'], figure3d)
    draw.meridian_distance_portions(sphere, pnt['acceptor'], figure3d)
    draw.point(sphere, pnt['promissor_aspect'], figure3d, line=True)
    draw.point(sphere, pnt['promissor'], figure3d, line=False)
    draw.point(sphere, pnt['acceptor'], figure3d, line=False)
    end_point_data = draw.direction_arc(
        sphere,
        pnt['promissor_aspect'],
        pnt['acceptor'],
        0, figure3d)
    # draw.ecl_projection(sphere, end_point_data, figure3d)
    draw.zodiac_2d(sphere, [end_point_data, pnt['promissor_aspect']], figure2d)
