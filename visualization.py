"""
Visualize celestial sphere
"""
from datetime import datetime
from components.visualization.show import visualize
from components.visualization.custom_templates.regio_zodiac_directions import figures

# Set spacial data
DATA = dict(
    init_datetime=datetime(2022, 10, 10, 23, 20),
    time_zone=4,
    geo_lon=44 + 46/60,
    geo_lat=56 + 43/60,
    # Set mode for the output figure
    animate=False
)

visualize(DATA, figures)
