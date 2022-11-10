"""
Functions for matplotlib visualisation of
Placidus house system and primary directions.
"""
import numpy as np
from matplotlib.axes import Axes
from sphere import Sphere
from primary_directions import Directions
from components.visualization import draw_sphere as draw
from components.visualization.draw_astrology_common import direction_arc as arc


def placidus_schema(sphere: Sphere, axs: Axes) -> None:
    """
    Illustrates the principle behind the
    Placidus house system.
    """
    ramc = sphere.ramc
    for dec in range(-90, 90):
        dsa = sphere.dsa(dec)
        if dsa is not None:
            rasc = np.linspace(ramc, ramc + dsa/3)
            plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz()
            axs.plot(*plot_data.aslist(), linewidth=0.7, color='lightgray')


def placidus(sphere: Sphere, axs: Axes, under_horizon: bool = False):
    """
    Draw placidus house lines on sphere
    """
    ramc = sphere.ramc
    _x = [[], [], [], [], [], ]
    _y = [[], [], [], [], [], ]
    _z = [[], [], [], [], [], ]
    for dec in range(-90, 90):
        dsa = sphere.dsa(dec)
        if dsa is not None:
            for i in range(0, 5):
                if under_horizon:
                    vector = sphere.set_equatorial(
                        ramc + 180 + (i - 2) * dsa/3,
                        -dec).horizontal_xyz
                else:
                    vector = sphere.set_equatorial(
                        ramc + (i - 2) * dsa/3,
                        dec).horizontal_xyz()
                _x[i].append(vector.x)
                _y[i].append(vector.y)
                _z[i].append(vector.z)

    for i in range(0, 5):
        axs.plot(_x[i], _y[i], _z[i], linewidth=0.7, color='green')


def mundane_positions(sphere: Sphere,
                      acceptor_data: dict,
                      axs: Axes) -> None:
    """
    Draws mundane positions of the acceptor in
    Placidus house system
    """
    # Set acceptor
    directions = Directions(sphere)
    acceptor = sphere.set_ecliptical(
        acceptor_data['lon'], acceptor_data['lat'])

    # Draw acceptor point
    draw.point(sphere, acceptor_data, axs)

    # Find mundane positions of the acceptor
    mund_positions = directions.aspect_positions_placidus_mundane(acceptor)
    if not mund_positions:
        return None
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
    return None


def meridian_distance_portions(sphere: Sphere,
                               acceptor_data: dict,
                               axs: Axes) -> None:
    """
    Draws meridian distance portions
    of the acceptor
    """
    # Set acceptor
    directions = Directions(sphere)
    acceptor_eqt = sphere.set_ecliptical(
        acceptor_data['lon'],
        acceptor_data['lat']
    ).equatorial()
    acc_rasc = acceptor_eqt.rasc
    acc_dec = acceptor_eqt.dec

    # Draw MDPs
    acceptor_quadrant = directions.quadrant(acc_rasc, acc_dec)
    acceptor_mdp = directions.md_portion(acc_rasc, acc_dec)
    if acceptor_mdp is None:
        return None

    if acceptor_quadrant == 0:
        start_point = sphere.ramc
        ratio = acceptor_mdp
    elif acceptor_quadrant == 1:
        start_point = sphere.ramc + 180
        ratio = -1 * acceptor_mdp
    elif acceptor_quadrant == 2:
        start_point = sphere.ramc + 180
        ratio = acceptor_mdp
    else:
        start_point = sphere.ramc
        ratio = -1 * acceptor_mdp

    for dec in range(-90, 90):
        dsa = sphere.dsa(dec)
        if dsa is None:
            continue
        path = dsa if acceptor_quadrant in [0, 3] else 180 - dsa
        rasc = np.linspace(start_point, start_point + path * ratio)
        plot_data = sphere.set_equatorial(rasc, dec).horizontal_xyz()
        axs.plot(*plot_data.aslist(), linewidth=0.7, color='lightgray')
    return None


def direction_arc(sphere: Sphere,
                  promissor_data: dict,
                  acceptor_data: dict,
                  aspect: int,
                  axs: Axes) -> None:
    """
    Draws directional arc between promissor
    and acceptor's aspect in Placidus system
    """
    return arc(
        sphere, promissor_data, acceptor_data, aspect, axs, system='P'
    )
