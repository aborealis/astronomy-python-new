"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw_sphere
from components.visualization import draw_astrology_common
from components.visualization import draw_placidus

ASPECT = 90


def points() -> dict[dict]:
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

    return dict(
        acceptor=acceptor,
        promissor=promissor,
        label="MC"
    )


def figures(sphere: Sphere, figure3d: Axes, figure2d: Axes) -> None:
    """
    Call draw functions to include points
    and objects to visualization on celestial
    sphere and Zodiac circle
    """

    pnt = points()
    draw_sphere.surface(sphere, figure3d)
    draw_sphere.ecliptic(sphere, figure3d)
    draw_sphere.horizon(sphere, figure3d)
    draw_sphere.equator(sphere, figure3d)
    draw_placidus.mundane_positions(sphere, pnt['acceptor'], figure3d)
    draw_placidus.meridian_distance_portions(sphere, pnt['acceptor'], figure3d)
    draw_sphere.point(sphere, pnt['promissor'], figure3d, line=False)
    end_point_data = draw_placidus.direction_arc(
        sphere,
        pnt['promissor'],
        pnt['acceptor'],
        ASPECT, figure3d)
    if end_point_data is not None:
        draw_placidus.meridian_distance_portions(
            sphere, end_point_data, figure3d)
        draw_astrology_common.zodiac_2d(
            sphere, [end_point_data, pnt['acceptor']], figure2d)
