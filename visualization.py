from datetime import datetime, timedelta
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from vector import Vector
import components.time as tm
import components.draw as draw

# Set spacial data
init_datetime = datetime(2022, 10, 10, 0, 0)
time_zone = 4
geo = dict(
    geo_lon=44 + 46/60,
    geo_lat=66 + 43/50,
)

# Set mode for the output figure
ANIMATE = True


def draw_frame(minute: int, ax: Axes, ax2: Axes) -> None:
    """
    Draws celestial sphere for a specific minute
    """
    naive_datetime = init_datetime + timedelta(minutes=minute)
    vector = Vector(naive_datetime, time_zone, **geo)

    # Set as many points as you wish
    # in ecliptical coordinates
    point_asc = dict(
        lon=vector.asc(),
        lat=0,
        label="ASC",
        color="red",
    )
    point_mc = dict(
        lon=vector.mc(),
        lat=0,
        label="MC",
        color="orange",
    )

    # Uncomment the shapes you want to display
    draw.surface(vector, ax)
    draw.ecliptic(vector, ax)
    draw.point(vector, point_mc, ax)
    draw.point(vector, point_asc, ax)
    # draw.semiarc(vector, point_asc, ax, 'UMD')
    # draw.semiarc(vector, point_asc, ax, 'DSA')
    draw.eqt_projection(vector, point_asc, ax)
    draw.ecl_projection(vector, point_asc, ax)
    # draw.placidus(vector, ax)
    draw.circle_2d(vector, [point_asc, point_mc], ax2)

    ax.set_axis_off()
    ax2.set_axis_off()
    zoom = 1.5
    ax.set_xlim3d(-1/zoom, 1/zoom)
    ax.set_zlim3d(-1/zoom, 1/zoom)
    ax.set_ylim3d(-1/zoom, 1/zoom)
    ax.set_title(
        (f'{tm.__localtime__(naive_datetime, time_zone)} '
         + f'lon: {vector.show_degrees_minutes(geo["geo_lon"])} '
         + f'lat: {vector.show_degrees_minutes(geo["geo_lat"])}')
    )
    ax.set_aspect('equal')
    ax2.set_aspect('equal')
    ax.legend()

    if ANIMATE:
        plt.pause(1e-4)
        ax.cla()
        ax2.cla()
    else:
        plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax2 = fig.add_subplot(444)
for minute in range(0, 24 * 60 if ANIMATE else 2, 2):
    draw_frame(minute, ax, ax2)
