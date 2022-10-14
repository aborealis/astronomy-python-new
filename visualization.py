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

# In ecliptical coordinates
point_coordinates = dict(
    lon=86,
    lat=10,
)

vector = Vector(**local_coordinates)

# Create canvas
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_axis_off()

# Draw specific elements
draw.surface(vector, ax)
draw.ecliptic(vector, ax)
draw.point(vector, point_coordinates, ax)
draw.semiarc(vector, point_coordinates, ax, 'UMD')
draw.semiarc(vector, point_coordinates, ax, 'DSA')
draw.eqt_projection(vector, point_coordinates, ax)
draw.ecl_projection(vector, point_coordinates, ax)

# Set an equal aspect ratio
ax.set_aspect('equal')
ax.legend()

plt.show()
