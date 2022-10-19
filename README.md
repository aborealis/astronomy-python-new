# Astronomy of the celestial sphere

This package calculates and visualizes astronomical phenomena on the celestial sphere, like the real ascension of MC, the oblique ascension of a given stellar point, even the cusps of the popular astrological Placidus house system, and much more.

It is a handy tool for astronomy teachers, scientists, and astrologers.

## Quick start

1. Clone this repository in a separate directory.
2. Within that directory, run `pipenv install`. It will install the required packages into a virtual environment. Note that you should have a python pipenv module installed on your computer.
3. Run `pipenv shell` to enter the newly created virtual environment.

## The usage

Import the Vector class and init the vector object with the current observer's local time and geo coordinates. Remember to put negative values for western longitudes and southern latitudes.

```python
from vector import Vector

# Initiate a vector object
test_vector = Vector(
    datetime(2022, 10, 10, 4, 20),
    time_zone=4,
    geo_lon=44 + 46/60,
    geo_lat=66 + 43/50,
)
```

Now you can access to general information about celestial sphere

```python
# Local Sidereal Time, LST in hours
print(test_vector.lst)

# Human-readible format of LST
from datetime import timedelta
print(timedelta(hours=test_vector.lst))

# Real Ascention of Medium Coeli, RAMC
print(test_vector.ramc())

# Ascending Zodiac degree, ASC
print(test_vector.asc())

# Culminating Zodiac degree, MC
print(test_vector.medium_coeli())

# The Zodiac degree of astrological
# House Cusp in Placidus System. It
# works even at extreme latitudes for
# the regions of the ascendible part
# of the ecliptic plane.
print(test_vector.placidus(cusp=11))
```

Now you can specify the stellar point in any coordibate system:

```python
# Set the ASC in ecliptic plane
test_vector.set_ecliptical(lon=test_vector.asc(), lat=0.0)

# Set the True East point in horizontal coordinates
test_vector.set_horizontal(azm=90.0, alt=0.0)

# Set the RAMC on the celestial equator
test_vector.set_ecliptical(ra=test_vector.ramc(), dec=0.0)
```

After the stellar point is specified, you can get extra information on it:

```python
# Express the same point in different coordinate systems
test_vector.set_horizontal(azm=90.0, alt=0.0)
print('True East:', test_vector.horizontal())
print('True East:', test_vector.equatorial())
print('True East:', test_vector.ecliptical())

# Transfer spherical coordinates into cartesian system
# in the same plane
test_vector.set_equatorial(ra=124.3, dec=-3.2)
print('Equatorial XYZ:', test_vector.cartesian())

# Or you can get cartesian coordinates in any other plane
zodiac_plane = test_vector.ecliptical()
print('Zodiac XYZ:', 
       test_vector.cartesian(zodiac_plane))

# Get an oblique ascention (OA), ascention difference (AD),
# diurnal semiarc (DSA) and upper meridian distance (UMD)
# of a given point
oblique_asc = test_vector.oblique_asc()
ascention_diff = test_vector.ascension_diff()
diurnal_semiarc = test_vector.dsa()
upper_meridian_distance = test_vector.umd()


# Alternatively, you can get the same values
# for any declination
test_vector.ascension_diff(dec=10.2)
```

## Visualization

See the sample [https://youtu.be/mKVUZVK_6NM](https://youtu.be/mKVUZVK_6NM)

You can use `matplotlib` or any other library to visualize these values. I have created a short module `visualization.py`. To use it:

### Set initial data

Set spacial data and choose if you want an animation or a static 3d-figure

```python
# Set spacial data
init_datetime = datetime(2022, 10, 10, 0, 0)
TIME_ZONE = 4
geo = dict(
    geo_lon=44 + 46/60,
    geo_lat=76 + 43/50,
)

# Set mode for the output figure
ANIMATE = True
```

### Add celestial points

Within `draw_frame()` function add as many points as you wish **in ecliptical coordinates**, like
```python
point_asc = dict(
    lon=vector.asc(),
    lat=0,
    label="ASC",
    color="red",
)
```

### Draw/hide element

There are several valuable functions to draw elements in 3d-axes (referred to as `ax1`) and in 2d-axes (a projection of a zodiac circle, referred to as `ax2`)

```python

# Draw the surface of the sphere in 3d-axes
draw.surface(vector, ax1)

# Draw ecliptic plane in 3d-axes
draw.ecliptic(vector, ax1)

# Draw a chosen point in 3d-axes
draw.point(vector, point_asc, ax1)

# Draw projection of a chosen point to
# ecliptic and equator
draw.eqt_projection(vector, point_asc, ax1)
draw.ecl_projection(vector, point_asc, ax1)

# Draw Placidus dome curves
draw.placidus(vector, ax1)

# Represent the points in 2d-axes
draw.circle_2d(vector, [
    point_asc, point2, point3,
], ax2)
```

Change the list of called functions within the `draw_frame()` to set up the visualization of your choice.