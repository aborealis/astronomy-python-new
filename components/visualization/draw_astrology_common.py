"""
Shared functions for matplotlib visualisation
of primary directions and zodiacal positions
not connected to specific house system.
"""
from typing import Optional
import numpy as np
from matplotlib.axes import Axes
from sphere import Sphere
from primary_directions import Directions


def zodiac_2d(sphere: Sphere, points: list[dict], axs: Axes) -> None:
    """
    Draws 2d Zodiac circle
    """
    alpha = np.linspace(0, 2 * np.pi)
    _x = np.cos(alpha)
    _y = np.sin(alpha)
    axs.plot(_x, _y, color=(.7, .2, .7))
    axs.text(-0.95, -0.1, s='0', fontsize=6)
    axs.text(-0.1, -0.9, s='90', fontsize=6)
    for _p in points:
        vector = sphere.set_ecliptical(
            lon=_p['lon'],
            lat=_p['lat']).ecliptical_xyz()
        axs.plot([-vector.x, 0], [-vector.y, 0], color='red',
                 linewidth=1, linestyle='dotted')
        axs.scatter(-vector.x, -vector.y, color=_p['color'])


def direction_arc(sphere: Sphere,
                  promissor_data: dict,
                  acceptor_data: dict,
                  aspect: int,
                  axs: Axes,
                  system: str = 'P') -> None:
    """
    Draws directional arc between promissor
    and acceptor's aspect in chosen house system
    """
    # Set promissor and significator
    directions = Directions(sphere)
    promissor = sphere.set_ecliptical(
        promissor_data['lon'],
        promissor_data['lat'])
    acceptor = sphere.set_ecliptical(
        acceptor_data['lon'],
        acceptor_data['lat'])

    promissor_eqt = promissor.equatorial()
    prom_rasc = promissor_eqt.rasc
    prom_dec = promissor_eqt.dec

    # Get directional arc
    if system == 'P':
        all_directions = directions.placidus_mundane(
            promissor,
            acceptor,
            aspect
        )
    else:
        all_directions = directions.regio_mundane(
            promissor,
            acceptor,
            aspect
        )

    if all_directions is None:
        return None

    # We'll take a positive arc to exclude reverse direction
    # Alternatively, we may choose the arc with min abs(len)
    arc = sorted([
        item['dist'] for item in all_directions
    ], reverse=True)[0]

    # Draw the arc
    rasc = np.linspace(prom_rasc - arc, prom_rasc)
    plot_data = sphere.set_equatorial(rasc, prom_dec).horizontal_xyz()
    axs.plot(*plot_data.aslist(),
             label='direction',
             linewidth=1.2,
             color="red")

    # Draw end point of directional arc on sphere
    xyz_hrz = sphere.set_equatorial(prom_rasc - arc, prom_dec).horizontal_xyz()
    _x, _y, _z = xyz_hrz.aslist()
    axs.plot([_x, 0], [_y, 0], [_z, 0], color='red',
             linewidth=1, linestyle='dotted')
    axs.scatter(
        _x, _y, _z,
        color="grey",
        label=None
    )

    # Draw end point of directional arc on Zodiac
    end_point = sphere.set_equatorial(
        prom_rasc - arc,
        prom_dec
    ).ecliptical()
    end_point_data = dict(
        lon=end_point.lon,
        lat=end_point.lat,
        label='Aspect',
        color="grey",
    )

    return end_point_data


def zodiac_positions(sphere: Sphere,
                     promissor_data: dict,
                     axs: Axes,
                     field_plane_lat: Optional[float] = 0,
                     ) -> None:
    """
    Draws zodiac positions of the promissor
    """
    # Set acceptor
    directions = Directions(sphere)
    promissor = sphere.set_ecliptical(
        promissor_data['lon'], promissor_data['lat'])

    # Find aspect positions of the acceptor
    aspect_positions = directions.aspect_positions_zodiac(promissor)
    if not aspect_positions:
        return None
    aspect_positions = [
        item['lon'] for item in aspect_positions
        if item['aspect']
    ]

    # Draw promissor's aspect positions
    xyz_hrz = sphere.set_ecliptical(
        promissor_data['lon'], field_plane_lat).horizontal_xyz()
    _x0, _y0, _z0 = xyz_hrz.aslist()

    for _m in aspect_positions:
        xyz_hrz = sphere.set_ecliptical(_m, field_plane_lat).horizontal_xyz()
        _x, _y, _z = xyz_hrz.aslist()
        axs.plot([_x, _x0], [_y, _y0], [_z, _z0], color=(.7, .2, .7),
                 linewidth=0.8, linestyle='solid')
        axs.scatter(
            _x, _y, _z,
            color=(.7, .2, .7),
            label=None
        )
    return None
