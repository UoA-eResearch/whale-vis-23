from datetime import datetime, timedelta

import pandas as pd
from bokeh.models import GeoJSONDataSource, CategoricalColorMapper, LinearColorMapper, CDSView, BooleanFilter, ColorBar, \
    Legend, ColorMapper
from bokeh.palettes import Inferno256, Bright5
from bokeh.plotting import figure
from colorcet import glasbey_cool
from geopandas import GeoDataFrame

from plotting.annotations import north_arrow, scale_bar, date_annotation, add_logo
from utils import timer

pd.options.mode.chained_assignment = None


def _fig_size(bounds, plot_width=1200):
    """Calculates figure size required to give square axes"""
    data_width, data_height = bounds[2] - bounds[0], bounds[3] - bounds[1]

    return plot_width, int(plot_width * data_height / data_width)


def _fade(ts, plot_ts, cutoff, max_alpha=0.5):
    """Quad ease out fade function, based on how long since point was plotted"""
    diff = (plot_ts - ts).total_seconds()
    if 0 <= diff <= cutoff:
        return max_alpha * (1 - (diff / cutoff) ** 4)
    else:
        return 0


def whale_colormap(whale_df):
    """Color whales by name"""
    whale_names = whale_df['name'].unique()
    # inds = np.floor(np.linspace(0, 255, len(whale_names))).astype(int)
    # colors = [Viridis256[i] for i in inds]
    colors = [glasbey_cool[i] for i in range(len(whale_names))]
    return CategoricalColorMapper(factors=sorted(whale_names), palette=colors)


def vessel_colormap():
    """Color vessels by type"""
    vessel_types = ['Fishing', 'Passenger', 'Cargo', 'Tanker', 'Other']
    vessel_colors = {t: c for t, c in zip(vessel_types, Bright5)}
    return CategoricalColorMapper(factors=vessel_types, palette=Bright5), vessel_colors


def add_whale_legend(fig: figure, cmapper: ColorMapper, names: list[str]) -> int:
    def _get_from_cmap(name):
        idx = cmapper.factors.index(name)
        return cmapper.palette[idx]

    # Split legend into multiple lines
    per_line = 6
    num_lines = int(len(names) / per_line)
    start_y = 30 * num_lines + 2
    for idx in range(num_lines+1):
        start = idx * per_line
        stop = min((idx + 1) * per_line, len(names))
        sub_names = names[start:stop]
        # Create legend line
        legend_dummies = {
            name: fig.line([0, 0], [0, 0], color=_get_from_cmap(name), line_width=8)
            for name in sub_names
        }
        whale_legend = Legend(items=[(name, [legend_dummies[name]]) for name in sub_names],
                              location=(2, start_y - idx * 30), orientation='horizontal')
        fig.add_layout(whale_legend)

    # Space taken up by whale legends
    return num_lines * 30

def plot_whale_pts(fig: figure, whale_df: GeoDataFrame, timestamp: datetime = None) -> int:
    """Add whale points to plot, optionally up to a given timestamp"""
    cmapper = whale_colormap(whale_df)
    cmap = {t: c for t, c in zip(whale_df['name'].unique(), cmapper.palette)}

    if timestamp:
        mask = whale_df['timestamp'] <= timestamp
        mask &= (timestamp - timedelta(days=14)) < whale_df['timestamp']
        if mask.sum() == 0:
            return
        for name, group in whale_df[mask].groupby('name'):
            fig.line(group.geometry.x, group.geometry.y, color=cmap[name], line_width=2, line_alpha=0.8)

        return add_whale_legend(fig, cmapper, sorted(whale_df[mask]['name'].unique()))
    else:
        for name, group in whale_df.groupby('name'):
            fig.line(group.geometry.x, group.geometry.y, color=cmap[name], line_width=2, line_alpha=0.8)

        return add_whale_legend(fig, cmapper, sorted(whale_df['name'].unique()))


def plot_whales_fade(fig: figure, whale_seg_df: GeoDataFrame, timestamp: datetime) -> int:
    """Plot whale traces, fading out after a cutoff"""
    mask = whale_seg_df['timestamp'] <= timestamp
    mask &= whale_seg_df['timestamp'] > (timestamp - timedelta(days=14))
    if mask.sum() == 0:
        return

    whale_data = whale_seg_df[mask]
    whale_data['fade'] = whale_data['timestamp'].apply(_fade, plot_ts=timestamp, cutoff=14 * 24 * 3600, max_alpha=0.8)

    cmapper = whale_colormap(whale_seg_df)
    src = GeoJSONDataSource(geojson=whale_data.drop(columns=['timestamp']).to_json())
    fig.multi_line('xs', 'ys', source=src, color={'field': 'name', 'transform': cmapper}, line_width=3,
                   line_alpha='fade')

    return add_whale_legend(fig, cmapper, sorted(whale_data['name'].unique()))


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


def plot_basemap(fig: figure, basemap: GeoJSONDataSource):
    """Add basemap (coast outline) to plot"""
    fig.patches('xs', 'ys', source=basemap, fill_color='lightgreen', line_alpha=1, fill_alpha=1)


def plot_protected_areas(fig: figure, protected_areas: GeoDataFrame, legend_offset=0):
    """Add marine protected areas to plot"""
    protected_source = GeoJSONDataSource(geojson=protected_areas.to_json())
    imma_view = CDSView(filter=BooleanFilter((protected_areas['ptype'] == 'imma').to_list()))
    mpa_view = CDSView(filter=BooleanFilter((protected_areas['ptype'] == 'mpa').to_list()))
    imma_patches = fig.patches('xs', 'ys', source=protected_source, fill_color='#ddd', line_color='#ddd',
                               line_alpha=1, fill_alpha=0.5, view=imma_view)
    mpa_patches = fig.patches('xs', 'ys', source=protected_source, fill_color='lightskyblue', line_color='lightskyblue',
                              line_alpha=1, fill_alpha=0.5, view=mpa_view)

    pa_legend = Legend(items=[('MPA', [mpa_patches]), ('IMMA', [imma_patches])],
                       location=(115, 2 + legend_offset), orientation='vertical')
    fig.add_layout(pa_legend)


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
    plot_basemap(fig, GeoJSONDataSource(geojson=basemap.to_json(default=str)))

    # Add annotations
    north_arrow(fig)
    scale_bar(fig, convert_from_deg=whale_df.crs.equals(4326))
    add_logo(fig, 'assets/logo.png', bounds, x_pos=0.95, y_pos=0.98, anchor='top_right')
    if timestamp:
        date_annotation(fig, timestamp)

    zoom_to_bounds(fig, bounds)

    return fig


def plot_encounters(fig: figure, vessel_pts: GeoDataFrame, max_dist=20000, timestamp=None, legend_offset=0):
    """Add encounter scatter/heatmap to figure"""
    mask = (~vessel_pts['encounter_dist'].isna()) & (vessel_pts['encounter_dist'] < max_dist)
    if timestamp:
        mask &= vessel_pts['timestamp'] <= timestamp
        mask &= vessel_pts['timestamp'] > (timestamp - timedelta(days=14))

    cmap = LinearColorMapper(Inferno256, low=max_dist, high=0)

    if mask.sum() > 0:
        vessel_data = (
            vessel_pts[mask]  # Only points with encounters
            .sort_values('encounter_dist')  # Sort by distance so that closest points are plotted last
            .iloc[::-1]  # Reverse order
        )
        vessel_source = GeoJSONDataSource(geojson=vessel_data.to_json(default=str))

        fig.scatter('x', 'y', source=vessel_source, color={'field': 'encounter_dist', 'transform': cmap},
                    fill_alpha=0.4, size=10, line_color=None)

    # Add colorbar
    color_bar = ColorBar(color_mapper=cmap, title='Encounter distance (m)')
    fig.add_layout(color_bar, 'right')

    # Add legend
    legend_dummies = {
        'Encounter': fig.scatter([0], [0], color='yellow', fill_alpha=0.8, size=10, line_color=None)
    }
    encounter_legend = Legend(items=[('Encounter', [legend_dummies['Encounter']])], location=(115, 62 + legend_offset), orientation='vertical')
    fig.add_layout(encounter_legend)


def encounters_map(whale_df: GeoDataFrame, vessel_df: GeoDataFrame, vessel_pts: GeoDataFrame,
                   protected_areas: GeoDataFrame, basemap: GeoDataFrame, bounds, max_dist=20000):
    """Produce a plot showing encounters between vessels & whales"""
    # Base figure containing basemap, MPAs, vessel & whale traces
    plot_width, plot_height = _fig_size(bounds)
    fig = figure(frame_width=plot_width, frame_height=plot_height, output_backend='webgl', toolbar_location=None)

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
        whale_cmap = whale_colormap(whale_pts_df)
        whale_data = whale_pts_df[whale_pts_df['timestamp'] == timestamp]
        whale_source = GeoJSONDataSource(geojson=whale_data.to_json(default=str))
        fig.scatter('x', 'y', source=whale_source, color={'field': 'name', 'transform': whale_cmap}, fill_alpha=1,
                    size=10)

    if (vessel_pts_df['timestamp'] == timestamp).sum() > 0:
        vessel_cmap, _ = vessel_colormap()
        vessel_data = vessel_pts_df[vessel_pts_df['timestamp'] == timestamp]
        vessel_source = GeoJSONDataSource(geojson=vessel_data.to_json(default=str))
        fig.scatter('x', 'y', source=vessel_source, color={'field': 'type', 'transform': vessel_cmap}, fill_alpha=1,
                    size=5)


def plot_vessels_fade(fig, vessels_seg_df, timestamp: datetime, legend_offset=0):
    """Plot vessel traces, fading out after a cutoff"""
    mask = vessels_seg_df['timestamp'] <= timestamp
    mask &= vessels_seg_df['timestamp'] > (timestamp - timedelta(days=14))

    cmap, vessel_colors = vessel_colormap()

    if mask.sum() > 0:
        vessels_data = vessels_seg_df[mask]
        vessels_data['fade'] = vessels_data['timestamp'].apply(_fade, plot_ts=timestamp, cutoff=14 * 24 * 3600,
                                                               max_alpha=0.5)

        src = GeoJSONDataSource(geojson=vessels_data.drop(columns=['timestamp']).to_json())
        fig.multi_line('xs', 'ys', source=src, color={'field': 'type', 'transform': cmap},
                       line_width=1, line_alpha='fade')

    # Include all vessel types in legend
    legend_dummies = {
        label: fig.line([0, 0], [0, 0], color=vessel_colors[label], line_width=2, line_alpha=0.5)
        for label in vessel_colors.keys()
    }

    vessel_legend = Legend(items=[(vtype, [legend_dummies[vtype]]) for vtype in vessel_colors.keys()],
                           location=(2, 2 + legend_offset), orientation='vertical')
    fig.add_layout(vessel_legend)


def plot_partial_vessel_traces(fig, vessels_pts_df, timestamp: datetime = None, legend_offset=0):
    """Plot vessel traces up to a given timestamp"""
    if timestamp is not None:
        mask = vessels_pts_df['timestamp'] <= timestamp
        mask &= vessels_pts_df['timestamp'] > (timestamp - timedelta(days=14))
    else:
        mask = slice(None)

    _, vessel_colors = vessel_colormap()

    if timestamp is None or mask.sum() > 0:
        vessel_data = vessels_pts_df[mask]

        for callsign, group in vessel_data.groupby('callsign'):
            vtype = group.iloc[0]['type']
            fig.line(group.geometry.x, group.geometry.y, color=vessel_colors[vtype], line_width=1, line_alpha=0.5)

    # Include all vessel types in legend
    legend_dummies = {
        label: fig.line([0, 0], [0, 0], color=vessel_colors[label], line_width=2, line_alpha=0.5)
        for label in vessel_colors.keys()
    }

    vessel_legend = Legend(items=[(vtype, [legend_dummies[vtype]]) for vtype in vessel_colors.keys()],
                           location=(2, 2 + legend_offset), orientation='vertical')
    fig.add_layout(vessel_legend)


def animation_frame(whales_df: GeoDataFrame, vessels_pts_df: GeoDataFrame, protected_areas: GeoDataFrame,
                    basemap_src: GeoJSONDataSource, bounds: list[float], timestamp: datetime = None,
                    encounters: GeoDataFrame = None):
    """
    Produce a plot showing the current location of vessels and whales
    If no timestamp is passed, plots all data given as a static map
    """
    plot_width, plot_height = _fig_size(bounds)
    fig = figure(frame_width=plot_width, frame_height=plot_height, toolbar_location=None, match_aspect=True)

    # Leave space for 2 lines of whale legend
    legend_offset = 60

    # Add layers
    with timer('plot_protected_areas'):
        plot_protected_areas(fig, protected_areas, legend_offset=legend_offset)
    with timer('plot_partial_vessel_traces'):
        plot_partial_vessel_traces(fig, vessels_pts_df, timestamp, legend_offset=legend_offset)
    with timer('plot_whale_pts'):
        plot_whale_pts(fig, whales_df, timestamp)
    if encounters is not None:
        with timer('plot_encounters'):
            plot_encounters(fig, encounters, max_dist=20000, timestamp=timestamp, legend_offset=legend_offset)
    with timer('plot_basemap'):
        plot_basemap(fig, basemap_src)
    if timestamp is not None:
        with timer('plot_location'):
            plot_location(fig, vessels_pts_df, whales_df, timestamp)

    fig.xaxis.axis_label = 'Lon'
    fig.yaxis.axis_label = 'Lat'

    # Add annotations
    north_arrow(fig, bounds, x_pos=0.03, y_pos=0.98, arrow_size=0.03)
    scale_bar(fig, convert_from_deg=whales_df.crs.equals(4326))
    logo_height = 60 / plot_height
    add_logo(fig, 'assets/logo.png', bounds, x_pos=0.98, y_pos=0.98, height=logo_height, anchor='top_right')
    add_logo(fig, 'assets/UoA_logo.png', bounds, x_pos=0.06, y_pos=0.98, height=0.8*logo_height, anchor='top_left')
    if timestamp is not None:
        date_annotation(fig, timestamp, bounds, x_pos=0.5, y_pos=0.02)

    zoom_to_bounds(fig, bounds)

    return fig


def animation_frame_fade(whales_seg_df: GeoDataFrame, vessels_seg_df: GeoDataFrame,
                         whales_pts_df: GeoDataFrame, vessels_pts_df: GeoDataFrame,
                         protected_areas: GeoDataFrame, basemap_src: GeoJSONDataSource,
                         bounds: list[float], timestamp: datetime, encounters: GeoDataFrame = None):
    """Produce a plot showing the current location of vessels and whales"""
    plot_width, plot_height = _fig_size(bounds)
    fig = figure(frame_width=plot_width, frame_height=plot_height, toolbar_location=None, match_aspect=True)

    # Leave room for one line of legend
    legend_offset = 30

    # Add layers
    with timer('plot_protected_areas'):
        plot_protected_areas(fig, protected_areas, legend_offset=legend_offset)
    with timer('plot_partial_vessel_traces'):
        plot_vessels_fade(fig, vessels_seg_df, timestamp, legend_offset=legend_offset)
    with timer('plot_whale_pts'):
        plot_whales_fade(fig, whales_seg_df, timestamp)
    if encounters is not None:
        with timer('plot_encounters'):
            plot_encounters(fig, encounters, max_dist=20000, timestamp=timestamp, legend_offset=legend_offset)
    with timer('plot_basemap'):
        plot_basemap(fig, basemap_src)
    with timer('plot_location'):
        plot_location(fig, vessels_pts_df, whales_pts_df, timestamp)

    fig.xaxis.axis_label = 'Lon'
    fig.yaxis.axis_label = 'Lat'

    # Add annotations
    north_arrow(fig, bounds, x_pos=0.03, y_pos=0.98, arrow_size=0.03)
    scale_bar(fig, convert_from_deg=whales_pts_df.crs.equals(4326))
    logo_height = 60 / plot_height
    add_logo(fig, 'assets/logo.png', bounds, x_pos=0.98, y_pos=0.98, height=logo_height, anchor='top_right')
    add_logo(fig, 'assets/UoA_logo.png', bounds, x_pos=0.06, y_pos=0.98, height=0.8*logo_height, anchor='top_left')
    date_annotation(fig, timestamp, bounds, x_pos=0.5, y_pos=0.02)

    zoom_to_bounds(fig, bounds)

    return fig
