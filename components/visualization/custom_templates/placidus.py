"""
Define points and figures to display
a Placidus house system.
"""
from matplotlib.axes import Axes
from sphere import Sphere
from components.visualization import draw_sphere
from components.visualization import draw_placidus
from components.visualization.draw_astrology_common import zodiac_2d


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

    cusps = []
    for num in [2, 3, 4, 5, 6, 7, 8, 9, 11, 12]:
        cusp = sphere.placidus(num)
        if cusp:
            cusps += cusp

    cusps = [
        dict(
            lon=lon,
            lat=0,
            label=None,
            color="grey",
        ) for lon in cusps
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
    draw_sphere.surface(sphere, figure3d)
    draw_sphere.horizon(sphere, figure3d)
    draw_sphere.equator(sphere, figure3d)
    draw_sphere.ecliptic(sphere, figure3d)
    draw_sphere.point(sphere, pnt['mc'], figure3d)
    draw_sphere.point(sphere, pnt['asc'], figure3d)

    for cusp in pnt['cusps']:
        draw_sphere.point(sphere, cusp, figure3d)

    # # draw.semiarc(sphere, pnt['asc'], figure3d, 'UMD')
    # # draw.semiarc(sphere, pnt['asc'], figure3d, 'DSA')
    draw_sphere.eqt_projection(sphere, pnt['asc'], figure3d)
    draw_sphere.ecl_projection(sphere, pnt['asc'], figure3d)
    draw_placidus.placidus(sphere, figure3d, under_horizon=False)
    zodiac_2d(sphere, [
        pnt['asc'], pnt['mc'], *pnt['cusps'],
    ], figure2d)
