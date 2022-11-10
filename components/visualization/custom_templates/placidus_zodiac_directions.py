"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw_sphere
from components.visualization import draw_placidus
from components.visualization.draw_astrology_common import zodiac_2d

ASPECT = 60
FIELD_PLANE_LAT = -10


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
        lat=FIELD_PLANE_LAT,
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
    draw_sphere.surface(sphere, figure3d)
    draw_sphere.ecliptic(sphere, figure3d)
    draw_sphere.horizon(sphere, figure3d)
    draw_sphere.equator(sphere, figure3d)
    draw_placidus.zodiac_positions_placidus(
        sphere, pnt['promissor'], figure3d, FIELD_PLANE_LAT)
    draw_placidus.meridian_distance_portions(sphere, pnt['acceptor'], figure3d)
    draw_sphere.point(sphere, pnt['promissor_aspect'], figure3d, line=True)
    draw_sphere.point(sphere, pnt['promissor'], figure3d, line=False)
    draw_sphere.point(sphere, pnt['acceptor'], figure3d, line=False)
    end_point_data = draw_placidus.direction_arc(
        sphere,
        pnt['promissor_aspect'],
        pnt['acceptor'],
        0, figure3d)
    zodiac_2d(
        sphere, [end_point_data, pnt['promissor_aspect']], figure2d)
