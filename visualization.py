from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from vector import Vector


vector = Vector(
    datetime(2022, 10, 10, 21, 58),
    time_zone=4,
    geo_lon=44.761088,
    geo_lat=41.713664
)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_axis_off()

# Celestial sphere surface
ra = np.linspace(0, 360, 30)
dec = np.linspace(0, 360, 30)
ra, dec = np.meshgrid(ra, dec)

vector.set_equatorial(ra, dec)
x = vector.cartesian_horizontal().x
y = vector.cartesian_horizontal().y
z = vector.cartesian_horizontal().z
ax.plot_wireframe(x, y, z, alpha=0.05)

# Horizon
azimuth = np.linspace(0, 2 * np.pi, 100)
z = 0
x = np.sin(azimuth)
y = np.cos(azimuth)
ax.plot(x, y, z, label='Horizon')

# Equator
ra = np.linspace(0, 360, 100)
vector.set_equatorial(ra, dec=0)
x = vector.cartesian_horizontal().x
y = vector.cartesian_horizontal().y
z = vector.cartesian_horizontal().z
ax.plot(x, y, z, label='Equator', linewidth=0.7, color="#c4d6e7")

# Equator degrees
for ra in range(0, 360, 30):
    vector.set_equatorial(ra, dec=0)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.text(x, y, z, s=str(ra), fontsize=6, color="cornflowerblue")

# Ecliptic
lon = np.linspace(0, 360, 100)
vector.set_ecliptical(lon, lat=0)
x = vector.cartesian_horizontal().x
y = vector.cartesian_horizontal().y
z = vector.cartesian_horizontal().z
ax.plot(x, y, z, label='Ecliptic', color=(.7, .2, .7))

# Ecliptic degrees
for lon in range(0, 360, 30):
    vector.set_ecliptical(lon, lat=0)
    x = vector.cartesian_horizontal().x
    y = vector.cartesian_horizontal().y
    z = vector.cartesian_horizontal().z
    ax.text(x, y, z, s=str(lon), fontsize=5, color=(.5, 0, .5))


cos_p = vector.__constants__.cos_p
sin_p = vector.__constants__.sin_p

# North - South Axis
ax.plot([0, 0], [1, -1], color='grey', linewidth=1, linestyle='dashed')
ax.text(0, -1, 0, s="S")
ax.text(0, 1, 0, s="N")

# East - West Axis
ax.plot([-1, 1], [0, 0], color='grey', linewidth=1, linestyle='dashed')
ax.text(-1, 0, 0, s="W")
ax.text(1, 0, 0, s="E")

# Polar Axis
ax.plot([0, 0], [0, cos_p], [0, sin_p],
        color='grey', linewidth=1, linestyle='dashed')

# Main Meridian
alpha = np.linspace(0, np.pi, 100)
z = np.sin(alpha)
x = np.sin(alpha) * 0
y = np.cos(alpha)
ax.plot(x, y, z, color='grey', linewidth=1,  linestyle='dashed')

vector.set_equatorial(ra=vector.ramc(), dec=0)
x = vector.cartesian_horizontal().x
y = vector.cartesian_horizontal().y
z = vector.cartesian_horizontal().z
ax.text(x, y, z, s="RAMC")

# Set an equal aspect ratio
ax.set_aspect('equal')
ax.legend()

plt.show()
