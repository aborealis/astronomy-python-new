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
        lon=vector.medium_coeli(),
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


def figures(vector: Vector, sphere: Axes, zodiac: Axes) -> None:
    """
    Call draw functions to include points
    and objects to visualization on celestial
    sphere and Zodiac circle
    """

    pnt = points(vector)
    draw.surface(vector, sphere)
    draw.ecliptic(vector, sphere)
    draw.mundane_positions_placidus(vector, pnt['acceptor'], sphere)
    draw.meridian_distance_portions(vector, pnt['acceptor'], sphere)
    draw.direction_arc(
        vector,
        pnt['promissor'],
        pnt['acceptor'],
        60, sphere, zodiac)
    draw.point(vector, pnt['mc'], sphere, line=False)
