# Astronomy of the celestial sphere

This package calculates and visualizes astronomical phenomena on the celestial sphere, like the right ascension of MC, the oblique ascension of a given stellar point, even the cusps of the popular astrological Placidus house system, and much more.

It is a handy tool for astronomy teachers, scientists, and astrologers.

## Quick start

1. Clone this repository in a separate directory.
2. Within that directory, run `pipenv install`. It will install the required packages into a virtual environment. Note that you should have a python pipenv module installed on your computer.
3. Run `pipenv shell` to enter the newly created virtual environment.

## The usage

Import the Sphere class and init the sphere object with the current observer's local time and geo coordinates. Remember to put negative values for western longitudes and southern latitudes.

```python
from sphere import Sphere
sphere = Sphere(
    datetime(2022, 10, 10, 4, 20),
    time_zone=4,
    geo_lon=44 + 46/60,
    geo_lat=66 + 43/60,
)
```

Now you can access to general information about celestial sphere

```python
# Local Sidereal Time, LST in hours
sphere.lst
>>> 4.5573066546517

# The same in hours and minutes
from datetime import timedelta
timedelta(hours=sphere.lst)
>>> 4:33:26.303957

# Right Ascention of Medium Coeli, RAMC, un degrees
sphere.ramc
>>> 68.3595998197755

# Ascending Zodiac degree, ASC
sphere.asc
>>> 168.27644339901482

# Culminating Zodiac degree, MC
sphere.medium_coeli
>>> 69.99780721838862

# The Zodiac degree of astrological
# House Cusp in Placidus System. It
# works even at extreme latitudes for
# the regions of the ascendible part
# of the ecliptic plane.
sphere.placidus(11)
>>> 116.92134279023504
```

Now you can specify the stellar point in any coordibate system:

```python
true_east = sphere.set_horizontal(azm=90, alt=0)
asc = sphere.set_ecliptical(lon=sphere.asc, lat=0)
aldebaran = sphere.set_equatorial(rasc=68.976, dec=16.507)
```

After the stellar point is specified, you can get extra information on it:

```python
# Express the same point in different coordinate systems
true_east.equatorial()
>>> Equatorial(rasc=158º 21'35"|158.3596º, dec=0º 00'00"|0.0000º)

true_east.ecliptical()
>>> Ecliptical(lon=159º 59'52"|159.9978º, lat=-8º 26'03"|-8.4343º)

# Right ascension of the True East
true_east.equatorial().rasc
>>> 158.3595998197755

# Oblique ascension of ascendant (the same as RA of the East)
asc.obl_asc()
>>> 158.35959981977547

# Ascension difference
asc.asc_diff()
>>> 10.86036453571822

# Upper Meridian Distance, UMD
cusp11 = sphere.set_ecliptical(lon=sphere.placidus(11), lat=0)
cusp11.umd()
>>> 50.603001275330115

# Diurnal semiarc, DSA. The .dsa() property.
# Check 1/3 of cusp11'th DSA:
cusp11.dsa() / 3
>>> 50.603001275330115
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

# Draw the surface of the sphere on 3d-figure
draw.surface(sphere, figure3d)

# Draw horizon plane on figure3d
draw.horizon(sphere, figure3d)

# Draw equator plane on figure3d
draw.equator(sphere, figure3d)

# Draw ecliptic plane on figure3d
draw.ecliptic(sphere, figure3d)

# Draw a chosen point on figure3d
draw.point(sphere, pnt['p1'], figure3d)

# Draw projection of a chosen point to
# ecliptic and equator
draw.eqt_projection(sphere, pnt['p1'], figure3d)
draw.ecl_projection(sphere, pnt['p1'], figure3d)

# Draw Placidus dome curves
draw.placidus(sphere, figure3d, under_horizon=False)

# Draw points on Zodiac circle on the top-right corner
draw.zodiac_2d(sphere, [
    pnt['p1'], pnt['p2'], pnt['p3'],
], figure2d)
```

### Unit converters

The Sphere class contains few methods to convert units. The `.dms()` converts float degree into the degrees minutes and seconds. The `.time_to_angle()` converts right ascention from time units into absolute degrees.

Example:
```python
# Set any star by known RA expressed in time units
aldebaran = sphere.set_equatorial(
    rasc=sphere.time_to_angle(hours=4, minutes=35, seconds=54.32),
    dec=16 + 30/60 + 26.1/3600
)
```