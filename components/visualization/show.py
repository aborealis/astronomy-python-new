"""
Draw static or animated picture of
celestial sphere and Zodiac circle
"""
from typing import Callable
from datetime import timedelta
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from sphere import Sphere
import components.vector.time as tm


def visualize(data: dict, include_figures: Callable):
    """
    Visualize the celestial sphere and Zodiac circle
    """
    fig = plt.figure()
    figure3d = fig.add_subplot(111, projection='3d')
    zodiac = fig.add_subplot(444)
    for _m in range(0, 24 * 60 if data['animate'] else 2, 2):
        draw_frame(data, include_figures, _m, figure3d, zodiac)


def draw_frame(data: dict,
               include_figures: Callable,
               minute: int,
               ax1: Axes, ax2: Axes) -> None:
    """
    Draws celestial sphere for a specific minute
    """
    naive_datetime = data['init_datetime'] + timedelta(minutes=minute)
    sphere = Sphere(
        naive_datetime, data['time_zone'], data['geo_lon'], data['geo_lat'])

    # Include figures from your custom template
    include_figures(sphere, ax1, ax2)

    ax1.set_axis_off()
    ax2.set_axis_off()
    zoom = 1.5
    ax1.set_xlim3d(-1/zoom, 1/zoom)
    ax1.set_zlim3d(-1/zoom, 1/zoom)
    ax1.set_ylim3d(-1/zoom, 1/zoom)
    ax1.set_title(
        (f'{tm.__localtime__(naive_datetime, data["time_zone"])} '
         + f'lon: {sphere.dms(data["geo_lon"])} '
         + f'lat: {sphere.dms(data["geo_lat"])}')
    )
    ax1.set_aspect('equal')
    ax2.set_aspect('equal')
    ax1.legend(loc=3)

    if data['animate']:
        plt.pause(1e-4)
        ax1.cla()
        ax2.cla()
    else:
        plt.show()
