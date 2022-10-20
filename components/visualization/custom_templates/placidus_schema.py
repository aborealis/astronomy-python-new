"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from vector import Vector
from components.visualization import draw


def points(vector: Vector) -> dict[dict]:
    """
    Returns a dictionary of points to display on
    celestial_sphere and on Zodiac circle
    """
    point_asc = dict(
        lon=vector.asc(),
        lat=0,
        label="ASC",
        color="red",
    )

    point_mc = dict(
        lon=vector.medium_coeli(),
        lat=0,
        label="MC",
        color="orange",
    )

    return dict(
        asc=point_asc,
        mc=point_mc,
    )


def figures(vector: Vector, sphere: Axes, zodiac: Axes) -> None:
    """
    Call draw functions to include points
    and objects to visualization on celestial
    sphere and Zodiac circle
    """

    pnt = points(vector)
    draw.surface(vector, sphere)
    draw.ecliptic(vector, sphere)
    draw.point(vector, pnt['mc'], sphere)
    draw.point(vector, pnt['asc'], sphere)

    draw.placidus(vector, sphere)
    draw.placidus_schema(vector, sphere)
