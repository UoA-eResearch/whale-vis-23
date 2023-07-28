import numpy as np
from bokeh.io import show
from bokeh.models import GeoJSONDataSource, CategoricalColorMapper
from bokeh.plotting import figure
from bokeh.palettes import Viridis256


def plot_traces(whale_df, vessel_df, basemap, bounds):
    whale_source = GeoJSONDataSource(geojson=whale_df.to_crs(2193).to_json())
    vessel_source = GeoJSONDataSource(geojson=vessel_df.to_crs(2193).to_json())
    basemap_source = GeoJSONDataSource(geojson=basemap.to_crs(2193).to_json())

    whale_names = whale_df['name'].unique()
    inds = np.floor(np.linspace(0, 255, len(whale_names))).astype(int)
    colors = [Viridis256[i] for i in inds]
    cmapper = CategoricalColorMapper(factors=whale_names, palette=colors)

    p = figure(width=1200, height=1200)

    p.patches('xs', 'ys', source=basemap_source, fill_color='lightgreen', line_alpha=1)

    p.multi_line('xs', 'ys', source=vessel_source, color='gray', line_width=1, alpha=0.2)

    w_lines = p.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, line_width=2)

    # Zoom to bounds
    p.x_range.start = bounds[0]
    p.x_range.end = bounds[2]
    p.y_range.start = bounds[1]
    p.y_range.end = bounds[3]

    show(p)