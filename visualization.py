"""
Visualize celestial sphere
"""
from datetime import datetime
from components.visualization.show import visualize
from components.visualization.custom_templates.placidus import figures

# Set spacial data
DATA = dict(
    init_datetime=datetime(2022, 10, 10, 0, 0),
    time_zone=4,
    geo_lon=44 + 46/60,
    geo_lat=76 + 43/50,
    # Set mode for the output figure
    animate=True
)

visualize(DATA, figures)
