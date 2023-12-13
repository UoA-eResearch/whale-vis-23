from datetime import datetime, timedelta
import numpy as np
from geopandas import GeoDataFrame
from bokeh.models import GeoJSONDataSource, CategoricalColorMapper, LinearColorMapper, CDSView, BooleanFilter, ColorBar
from bokeh.plotting import figure
from bokeh.palettes import Viridis256, Inferno256, Bright5

from plotting.annotations import north_arrow, scale_bar, date_annotation
from utils import timer


def _fig_size(bounds, plot_width=1200):
    """Calculates figure size required to give square axes"""
    data_width, data_height = bounds[2] - bounds[0], bounds[3] - bounds[1]
    x_margin, y_margin = 164, 28

    frame_width = plot_width - x_margin
    frame_height = int(frame_width * data_height / data_width)

    return plot_width, frame_height + y_margin


def _fade(ts, plot_ts, cutoff):
    """Quad ease out fade function, based on how long since point was plotted"""
    diff = (ts - plot_ts).total_seconds()
    if diff > cutoff:
        return 0
    else:
        return 1 - (diff / cutoff) ** 4


def whale_colormap(whale_df):
    """Color whales by name"""
    whale_names = whale_df['name'].unique()
    inds = np.floor(np.linspace(0, 255, len(whale_names))).astype(int)
    colors = [Viridis256[i] for i in inds]
    return CategoricalColorMapper(factors=whale_names, palette=colors)


def plot_whale_pts(fig: figure, whale_df: GeoDataFrame, timestamp: datetime = None):
    """Add whale points to plot, optionally up to a given timestamp"""
    # whale_df['fade'] = whale_df['timestamp'].apply(_fade, plot_ts=timestamp, cutoff=14*24*3600)
    whale_source = GeoJSONDataSource(geojson=whale_df.to_json(default=str))
    cmapper = whale_colormap(whale_df)

    if timestamp:
        mask = whale_df['timestamp'] <= timestamp
        mask &= (timestamp - timedelta(days=14)) < whale_df['timestamp']
        if mask.sum() == 0:
            return
        whale_view = CDSView(filter=BooleanFilter(mask.to_list()))
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
    imma_view = CDSView(filter=BooleanFilter((protected_areas['ptype'] == 'imma').to_list()))
    mpa_view = CDSView(filter=BooleanFilter((protected_areas['ptype'] == 'mpa').to_list()))
    fig.patches('xs', 'ys', source=protected_source, fill_color='#ddd', line_alpha=1, fill_alpha=0.8,
                view=imma_view)
    fig.patches('xs', 'ys', source=protected_source, fill_color='lightskyblue', line_alpha=1, fill_alpha=0.8,
                view=mpa_view)


def zoom_to_bounds(fig, bounds):
    """Crop a figure to given bounds"""
    fig.x_range.start = bounds[0]
    fig.x_range.end = bounds[2]
    fig.y_range.start = bounds[1]
    fig.y_range.end = bounds[3]


def traces_map(whale_df: GeoDataFrame, vessel_df: GeoDataFrame, protected_areas: GeoDataFrame,
               basemap: GeoDataFrame, bounds, timestamp=None):
    """Produce plot of basemap, MPAs, vessel traces and whale points"""
    fig = figure(width=1200, height=1200, output_backend='webgl', toolbar_location=None)

    # Add layers
    plot_protected_areas(fig, protected_areas)
    plot_vessel_traces(fig, vessel_df)
    plot_whale_pts(fig, whale_df, timestamp)
    plot_basemap(fig, basemap)

    # Add annotations
    north_arrow(fig)
    scale_bar(fig, convert_from_deg=whale_df.crs.equals(4326))

    zoom_to_bounds(fig, bounds)

    return fig


def plot_encounters(fig: figure, vessel_pts: GeoDataFrame, max_dist=20000):
    """Add encounter scatter/heatmap to figure"""
    mask = (~vessel_pts['whale_dist'].isna()) & (vessel_pts['whale_dist'] < max_dist)
    vessel_data = (
        vessel_pts[mask]  # Only points with encounters
        .sort_values('whale_dist')  # Sort by distance so that closest points are plotted last
        .iloc[::-1]  # Reverse order
    )
    vessel_source = GeoJSONDataSource(geojson=vessel_data.to_json())

    cmap = LinearColorMapper(Inferno256, low=vessel_data['whale_dist'].max(), high=0)

    fig.scatter('x', 'y', source=vessel_source, color={'field': 'whale_dist', 'transform': cmap},
                fill_alpha=0.4, size=10, line_color=None)

    # Add colorbar
    color_bar = ColorBar(color_mapper=cmap, title='Encounter distance')
    fig.add_layout(color_bar, 'right')


def encounters_map(whale_df: GeoDataFrame, vessel_df: GeoDataFrame, vessel_pts: GeoDataFrame,
                    protected_areas: GeoDataFrame, basemap: GeoDataFrame, bounds, max_dist=20000):
    """Produce a plot showing encounters between vessels & whales"""
    # Base figure containing basemap, MPAs, vessel & whale traces
    plot_width, plot_height = _fig_size(bounds)
    fig = figure(width=plot_width, height=plot_height, output_backend='webgl', toolbar_location=None)

    # Add layers
    plot_protected_areas(fig, protected_areas)
    plot_vessel_traces(fig, vessel_df)
    plot_whale_lines(fig, whale_df)
    plot_encounters(fig, vessel_pts, max_dist)
    plot_basemap(fig, basemap)

    # Add annotations
    north_arrow(fig)
    scale_bar(fig, convert_from_deg=whale_df.crs.equals(4326))

    zoom_to_bounds(fig, bounds)

    return fig


def plot_location(fig, vessel_pts_df, whale_pts_df, timestamp):
    if (whale_pts_df['timestamp'] == timestamp).sum() > 0:
        cmapper = whale_colormap(whale_pts_df)
        whale_data = whale_pts_df[whale_pts_df['timestamp'] == timestamp]
        whale_source = GeoJSONDataSource(geojson=whale_data.to_json(default=str))
        fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': cmapper}, fill_alpha=1, size=10)

    if (vessel_pts_df['timestamp'] == timestamp).sum() > 0:
        vessel_data = vessel_pts_df[vessel_pts_df['timestamp'] == timestamp]
        vessel_source = GeoJSONDataSource(geojson=vessel_data.to_json(default=str))
        fig.scatter('x', 'y', source=vessel_source, color='gray', fill_alpha=1, size=5, line_color=None)


def plot_partial_vessel_traces(fig, vessels_pts_df, timestamp):
    """Plot vessel traces up to a given timestamp"""
    mask = vessels_pts_df['timestamp'] <= timestamp
    mask &= vessels_pts_df['timestamp'] > (timestamp - timedelta(days=14))
    if mask.sum() == 0:
        return
    vessel_data = vessels_pts_df[mask]
    # vessel_data['fade'] = vessel_data['timestamp'].apply(_fade, plot_ts=timestamp, cutoff=14 * 24 * 3600)

    # TODO: consider building CDS for each callsign and streaming data in
    vessel_colors = {t: c for t, c in zip(vessels_pts_df['type'].unique(), Bright5)}
    for callsign, group in vessel_data.groupby('callsign'):
        vtype = group.iloc[0]['type']
        fig.line(group.geometry.x, group.geometry.y, color=vessel_colors[vtype], line_width=1, line_alpha=0.5,
                 legend_label=vtype)

    fig.legend.location = 'bottom_left'


def animation_frame(whales_df, vessels_pts_df, protected_areas, basemap, bounds, timestamp):
    """Produce a plot showing the current location of vessels and whales"""
    plot_width, plot_height = _fig_size(bounds)
    fig = figure(width=plot_width, height=plot_height, output_backend='webgl', toolbar_location=None)

    # Add layers
    with timer('plot_protected_areas'):
        plot_protected_areas(fig, protected_areas)
    with timer('plot_partial_vessel_traces'):
        plot_partial_vessel_traces(fig, vessels_pts_df, timestamp)
    with timer('plot_whale_pts'):
        plot_whale_pts(fig, whales_df, timestamp)
    with timer('plot_basemap'):
        plot_basemap(fig, basemap)
    with timer('plot_location'):
        plot_location(fig, vessels_pts_df, whales_df, timestamp)

    # Add annotations
    north_arrow(fig)
    scale_bar(fig, convert_from_deg=whales_df.crs.equals(4326))
    date_annotation(fig, timestamp)

    zoom_to_bounds(fig, bounds)

    return fig
