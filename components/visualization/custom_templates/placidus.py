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

    cusps = [
        dict(
            lon=vector.placidus(num),
            lat=0,
            label=None,
            color="grey",
        ) for num in [2, 3, 4, 5, 6, 7, 8, 9, 11, 12]
    ]

    cusps = [
        item for item in cusps
        if item['lon'] is not None
    ]

    return dict(
        asc=point_asc,
        mc=point_mc,
        cusps=cusps,
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

    for cusp in pnt['cusps']:
        draw.point(vector, cusp, sphere)

    # draw.semiarc(vector, point_asc, ax, 'UMD')
    # draw.semiarc(vector, point_asc, ax, 'DSA')
    draw.eqt_projection(vector, pnt['asc'], sphere)
    draw.ecl_projection(vector, pnt['asc'], sphere)
    draw.placidus(vector, sphere)
    # draw.placidus(vector, sphere, True)
    draw.zodiac_2d(vector, [
        pnt['asc'], pnt['mc'], *pnt['cusps'],
    ], zodiac)
