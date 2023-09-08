from datetime import datetime
import numpy as np
from geopandas import GeoDataFrame
from bokeh.models import GeoJSONDataSource, CategoricalColorMapper, LinearColorMapper, CDSView, BooleanFilter
from bokeh.plotting import figure
from bokeh.palettes import Viridis256, Inferno256


def whale_colormap(whale_df):
    """Color whales by name"""
    whale_names = whale_df['name'].unique()
    inds = np.floor(np.linspace(0, 255, len(whale_names))).astype(int)
    colors = [Viridis256[i] for i in inds]
    return CategoricalColorMapper(factors=whale_names, palette=colors)


def plot_whale_pts(fig: figure, whale_df: GeoDataFrame, timestamp: datetime = None):
    """Add whale points to plot, optionally up to a given timestamp"""
    whale_source = GeoJSONDataSource(geojson=whale_df.to_json(default=str))
    cmapper = whale_colormap(whale_df)

    if timestamp:
        whale_view = CDSView(filter=BooleanFilter((whale_df['date'] < timestamp).to_list()))
        fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, fill_alpha=0.5,
                    view=whale_view)
    else:
        fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, fill_alpha=0.5)


def plot_whale_lines(fig: figure, whale_df: GeoDataFrame):
    """Add coloured whale traces to plot"""
    whale_source = GeoJSONDataSource(geojson=whale_df.to_json(default=str))
    cmapper = whale_colormap(whale_df)

    fig.multi_line('xs', 'ys', source=whale_source, color={'field': 'name', 'transform': cmapper}, line_alpha=0.5,
                   line_width=3)


def plot_vessel_traces(fig: figure, vessel_df: GeoDataFrame):
    """Add vessel traces to plot"""
    vessel_source = GeoJSONDataSource(geojson=vessel_df.to_json(default=str))

    fig.multi_line('xs', 'ys', source=vessel_source, color='gray', line_width=1, line_alpha=0.05)


def plot_basemap(fig: figure, basemap: GeoDataFrame):
    """Add basemap (coast outline) to plot"""
    basemap_source = GeoJSONDataSource(geojson=basemap.to_json())
    fig.patches('xs', 'ys', source=basemap_source, fill_color='lightgreen', line_alpha=1, fill_alpha=1)


def plot_protected_areas(fig: figure, protected_areas: GeoDataFrame):
    """Add marine protected areas to plot"""
    protected_source = GeoJSONDataSource(geojson=protected_areas.to_json())
    fig.patches('xs', 'ys', source=protected_source, fill_color='lightblue', line_alpha=1, fill_alpha=0.5)


def plot_traces(whale_df: GeoDataFrame, vessel_df: GeoDataFrame, protected_areas: GeoDataFrame,
                basemap: GeoDataFrame, bounds, timestamp=None):
    """Produce plot of basemap, MPAs, vessel traces and whale points"""
    fig = figure(width=1200, height=1200, output_backend='webgl')

    # Add layers
    plot_protected_areas(fig, protected_areas)
    plot_vessel_traces(fig, vessel_df)
    plot_whale_pts(fig, whale_df, timestamp)
    plot_basemap(fig, basemap)

    # Zoom to bounds
    fig.x_range.start = bounds[0]
    fig.x_range.end = bounds[2]
    fig.y_range.start = bounds[1]
    fig.y_range.end = bounds[3]

    return fig


def plot_encounters(vessel_points: GeoDataFrame, fig: figure, max_dist=20000):
    mask = (~vessel_points['whale_dist'].isna()) & (vessel_points['whale_dist'] < max_dist)
    vessel_data = (
        vessel_points[mask]           # Only points with encounters
        .drop(columns=['timestamp'])  # Drop timestamp (not json serializable)
        .sort_values('whale_dist')    # Sort by distance so that closest points are plotted last
        .iloc[::-1]                   # Reverse order
    )
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