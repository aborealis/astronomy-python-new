"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw_sphere
from components.visualization import draw_placidus


def points(sphere: Sphere) -> dict[dict]:
    """
    Returns a dictionary of points to display on
    celestial_sphere and on Zodiac circle
    """
    point_asc = dict(
        lon=sphere.asc,
        lat=0,
        label="ASC",
        color="red",
    )

    point_mc = dict(
        lon=sphere.medium_coeli,
        lat=0,
        label="MC",
        color="orange",
    )

    return dict(
        asc=point_asc,
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
    draw_sphere.point(sphere, pnt['mc'], figure3d)
    draw_sphere.point(sphere, pnt['asc'], figure3d)

    draw_placidus.placidus(sphere, figure3d)
    draw_placidus.placidus_schema(sphere, figure3d)
