"""
Define points and figures to display
a Regiomontanus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw_sphere
from components.visualization import draw_regiomontanus
from components.visualization.draw_astrology_common import zodiac_2d

ASPECT = 60


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
    draw_sphere.point(sphere, pnt['promissor'], figure3d, line=False)
    draw_sphere.point(sphere, pnt['acceptor'], figure3d, line=False)
    draw_regiomontanus.mundane_positions(
        sphere, pnt['acceptor'], figure3d)
    end_point_data = draw_regiomontanus.direction_arc(
        sphere,
        pnt['promissor'],
        pnt['acceptor'],
        ASPECT, figure3d)
    draw_regiomontanus.regiomontanus_curve(sphere, end_point_data, figure3d)
    zodiac_2d(sphere, [end_point_data, pnt['acceptor']], figure2d)
