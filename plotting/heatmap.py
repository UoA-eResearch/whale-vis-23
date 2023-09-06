import numpy as np
from bokeh.models import GeoJSONDataSource, CategoricalColorMapper, LinearColorMapper, CDSView, BooleanFilter
from bokeh.plotting import figure
from bokeh.palettes import Viridis256, Inferno256


def plot_traces(whale_df, vessel_df, protected_areas, basemap, bounds, timestamp=None):
    whale_source = GeoJSONDataSource(geojson=whale_df.to_json(default=str))
    vessel_source = GeoJSONDataSource(geojson=vessel_df.to_json())
    protected_source = GeoJSONDataSource(geojson=protected_areas.to_json())
    basemap_source = GeoJSONDataSource(geojson=basemap.to_json())

    whale_names = whale_df['name'].unique()
    inds = np.floor(np.linspace(0, 255, len(whale_names))).astype(int)
    colors = [Viridis256[i] for i in inds]
    cmapper = CategoricalColorMapper(factors=whale_names, palette=colors)

    fig = figure(width=1200, height=1200)

    fig.patches('xs', 'ys', source=protected_source, fill_color='lightblue', line_alpha=1, fill_alpha=0.5)

    fig.patches('xs', 'ys', source=basemap_source, fill_color='lightgreen', line_alpha=1, fill_alpha=0.5)

    fig.multi_line('xs', 'ys', source=vessel_source, color='gray', line_width=1, line_alpha=0.1)

    if timestamp:
        whale_view = CDSView(filter=BooleanFilter((whale_df['date'] < timestamp).to_list()))
        fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, fill_alpha=0.5,
                  view=whale_view)
    else:
        fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, fill_alpha=0.5)

    # Zoom to bounds
    fig.x_range.start = bounds[0]
    fig.x_range.end = bounds[2]
    fig.y_range.start = bounds[1]
    fig.y_range.end = bounds[3]

    return fig


def plot_encounters(vessel_points, fig, max_dist=20000):
    vessel_data = vessel_points[~vessel_points['whale_dist'].isna()]
    vessel_data = vessel_data[vessel_data['whale_dist'] < max_dist]
    vessel_data.drop(columns=['timestamp'], inplace=True)
    vessel_source = GeoJSONDataSource(geojson=vessel_data.to_json())

    cmap = LinearColorMapper(Inferno256, low=vessel_data['whale_dist'].max(), high=0)

    fig.scatter('x', 'y', source=vessel_source, color={'field': 'whale_dist', 'transform': cmap},
                fill_alpha=0.2, size=10, line_color=None)

    return fig


def plot_location(fig, vessel_points, whale_points, timestamp):
    """Plot the current location of vessels and whales"""
    vessel_data = vessel_points[vessel_points['timestamp'] == timestamp]
    whale_data = whale_points[whale_points['timestamp'] == timestamp]

    vessel_source = GeoJSONDataSource(geojson=vessel_data.to_json(default=str))
    whale_source = GeoJSONDataSource(geojson=whale_data.to_json(default=str))

    whale_names = whale_points['name'].unique()
    inds = np.floor(np.linspace(0, 255, len(whale_names))).astype(int)
    colors = [Viridis256[i] for i in inds]
    cmapper = CategoricalColorMapper(factors=whale_names, palette=colors)

    fig.scatter('x', 'y', source=vessel_source, color='gray', fill_alpha=1, size=15, line_color=None)
    fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, fill_alpha=1, size=15)

    return fig