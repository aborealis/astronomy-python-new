"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw as draw


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

    cusps = [
        dict(
            lon=sphere.placidus(num),
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


def figures(sphere: Sphere, figure3d: Axes, figure2d: Axes) -> None:
    """
    Call draw functions to include points
    and objects to visualization on celestial
    sphere and Zodiac circle
    """

    pnt = points(sphere)
    draw.surface(sphere, figure3d)
    draw.horizon(sphere, figure3d)
    draw.equator(sphere, figure3d)
    draw.ecliptic(sphere, figure3d)
    draw.point(sphere, pnt['mc'], figure3d)
    draw.point(sphere, pnt['asc'], figure3d)

    for cusp in pnt['cusps']:
        draw.point(sphere, cusp, figure3d)

    # # draw.semiarc(sphere, pnt['asc'], figure3d, 'UMD')
    # # draw.semiarc(sphere, pnt['asc'], figure3d, 'DSA')
    draw.eqt_projection(sphere, pnt['asc'], figure3d)
    draw.ecl_projection(sphere, pnt['asc'], figure3d)
    draw.placidus(sphere, figure3d, under_horizon=False)
    draw.zodiac_2d(sphere, [
        pnt['asc'], pnt['mc'], *pnt['cusps'],
    ], figure2d)
