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

### Set initial data

Set spacial data and choose if you want an animation or a static 3d-figure. In `visualization.py` set the dictionary

```python
# Set spacial data
DATA = dict(
    init_datetime=datetime(2022, 10, 10, 0, 0),
    time_zone=4,
    geo_lon=44 + 46/60,
    geo_lat=76 + 43/50,
    # Set the mode for the output figure
    animate=True
)
```

### Choose what to include into the picture

Template is a short file which lists the elements to be included into the picture.

You can import one of the pre-defined templates from `visualization/custom_templates` directory or create your own template, and load it:

```python
# Include predefined template
from components.visualization.templates.placidus import figures
# Or include your own
from components.visualization.templates.my template import figures
```

### Create your own template

To create your own template, copy one of predefined templates and modify it.

First, modify the `points()` function to include your own points (planets, stars, galaxies) into the picture. All points are set in ecliptical coordinates:

```python
my_custom_point1 = dict(
    lon=22, # Ecliptical longitude
    lat=10, # Ecliptical latitude
    label="ASC", # Can be None
    color="red",
)

return dict(p1=my_custom_point1)
```

Second, modify the `figures()` function to include elements you want to show on sidereal sphere or Zodiac circle. You can use one of the following fumctions:

```python

# Draw the surface of the sphere on sphere
draw.surface(vector, sphere)

# Draw ecliptic plane on sphere
draw.ecliptic(vector, sphere)

# Draw a chosen point on sphere
draw.point(vector, pnt['p1'], sphere)

# Draw projection of a chosen point to
# ecliptic and equator
draw.eqt_projection(vector, pnt['p1'], sphere)
draw.ecl_projection(vector, pnt['p1'], sphere)

# Draw Placidus dome curves
draw.placidus(vector, sphere, under_horizon=False)

# Draw points on Zodiac circle on the top-right corner
draw.zodiac_2d(vector, [
    pnt['p1'], pnt['p2'], pnt['p3'],
], ax2)
```