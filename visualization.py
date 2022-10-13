from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from vector import Vector
import components.draw as draw

local_coordinates = dict(
    naive_datetime=datetime(2022, 10, 10, 21, 58),
    time_zone=4,
    geo_lon=44.761088,
    geo_lat=41.713664
)

vector_coordinates = dict(
    lon=286,
    lat=0,
)

vector = Vector(**local_coordinates)


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_axis_off()

draw.surface(vector, ax)
draw.ecliptic(vector, ax)
draw.point(vector, vector_coordinates, ax)
draw.semiarc(vector, vector_coordinates, ax, 'UMD')
draw.semiarc(vector, vector_coordinates, ax, 'DSA')

# Set an equal aspect ratio
ax.set_aspect('equal')
ax.legend()

plt.show()
