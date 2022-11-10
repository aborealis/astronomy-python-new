"""
Functions for matplotlib visualisation of
Regiomontanus house system and primary directions.
"""
import numpy as np
from matplotlib.axes import Axes
from sphere import Sphere
from primary_directions import Directions
from components.visualization.draw_astrology_common import direction_arc as arc


def __regiomontanus_curve_from_oa_utp(sphere: Sphere,
                                      axs: Axes,
                                      oa_utp: float,
                                      under_horizon: bool = False):
    """
    Draw regiomontanus curve for a given
    oblique ascension under the pole (oa_utp)
    """
    cos_p = sphere.__constants__.cos_p
    # angle delta between oa_asc and oa_utp
    delta = sphere.ramc + 90 - oa_utp
    # angle alpha betwen dividing-plane and eastern horizon
    alpha = np.arctan2(
        cos_p * np.sin(delta * np.pi / 180),
        np.cos(delta * np.pi / 180)
    )
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    # Parameter angle for the curve
    if under_horizon and alpha > 0:
        angle = np.linspace(np.pi / 2, np.pi * 3/2, 50)
    else:
        angle = np.linspace(-np.pi / 2, np.pi / 2, 50)
    _x = np.cos(angle) * cos_a
    _y = np.sin(angle)
    _z = np.cos(angle) * sin_a
    axs.plot(_x, _y, _z, linewidth=0.7, color='green')


def __regiomontanus_curve_from_eqt_coords(sphere: Sphere,
                                          axs: Axes,
                                          rasc: float,
                                          dec: float):
    """
    Draw regiomontanus curve for a given
    RA and dec
    """
    directions = Directions(sphere)
    oa_utp = directions.__oblique_asc_utp__(rasc, dec)

    if directions.quadrant(rasc, dec) in [1, 2]:
        __regiomontanus_curve_from_oa_utp(
            sphere, axs, oa_utp, under_horizon=True)
    else:
        __regiomontanus_curve_from_oa_utp(sphere, axs, oa_utp)


def regiomontanus_curve(sphere: Sphere,
                        point_data: dict,
                        axs: Axes):
    """
    Draw regiomontanus curve for a given
    RA and dec
    """
    point_eqt = sphere.set_ecliptical(
        point_data['lon'], point_data['lat']
    ).equatorial()
    __regiomontanus_curve_from_eqt_coords(
        sphere, axs,
        point_eqt.rasc, point_eqt.dec
    )


def regiomontanus(sphere: Sphere, axs: Axes, under_horizon: bool = False):
    """
    Draw regiomontanus house lines on sphere
    """
    for i in range(1, 6):
        oa_utp = sphere.ramc + 90 - 30 * i
        __regiomontanus_curve_from_oa_utp(sphere, axs, oa_utp, under_horizon)


def mundane_positions(sphere: Sphere,
                      acceptor_data: dict,
                      axs: Axes) -> None:
    """
    Draws mundane positions of the acceptor in
    Regiomontanus house system
    """
    # Set acceptor's oblique ascension under the pole
    directions = Directions(sphere)
    acceptor = sphere.set_ecliptical(
        acceptor_data['lon'], acceptor_data['lat'])
    acceptor_equatorial = acceptor.equatorial()

    __regiomontanus_curve_from_eqt_coords(
        sphere, axs, acceptor_equatorial.rasc, acceptor_equatorial.dec)

    # Find mundane positions of the acceptor
    mund_positions = directions.aspect_positions_regio_mundane(acceptor)
    conjunction = [
        item['rasc'] for item in mund_positions
        if item['aspect'] == 0
    ][0]
    mund_positions = [
        item['rasc'] for item in mund_positions
    ]

    # Draw acceptor's mundane position
    xyz_hrz = sphere.set_equatorial(conjunction, 0).horizontal_xyz()
    _x0, _y0, _z0 = xyz_hrz.aslist()

    for _m in mund_positions:
        xyz_hrz = sphere.set_equatorial(_m, 0).horizontal_xyz()
        _x, _y, _z = xyz_hrz.aslist()
        axs.plot([_x, _x0], [_y, _y0], [_z, _z0], color="#c4d6e7",
                 linewidth=1.5, linestyle='solid')
        axs.scatter(
            _x, _y, _z,
            color="#c4d6e7",
            label=None
        )


def direction_arc(sphere: Sphere,
                  promissor_data: dict,
                  acceptor_data: dict,
                  aspect: int,
                  axs: Axes) -> None:
    """
    Draws directional arc between promissor
    and acceptor's aspect in Regiomontanus system
    """
    return arc(
        sphere, promissor_data, acceptor_data, aspect, axs, system='R'
    )
