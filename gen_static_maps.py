import os
from os import path

import pandas as pd
import geopandas as gpd

from bokeh.io import export_png
from bokeh.io.webdriver import webdriver_control
from bokeh.models import GeoJSONDataSource
from tqdm import tqdm

from data_processing import load_data
from plotting import heatmap
from plotting.util import fix_dateline
from utils import timer


if __name__ == '__main__':
    # Generates static maps for all bounds, time ranges, with or without encounters

    # Load base data
    whales, vessels, protected_areas, basemap, bounds = load_data.load_all()

    # Load vessel points
    vessel_data_files = [path.join('data', 'vessels', f'{vessel_type}_points_coast.gpkg') for vessel_type in
                         ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]
    vessel_data_sets = [gpd.read_file(vdf) for vdf in tqdm(vessel_data_files, desc="Loading vessel data")]

    vessel_points = pd.concat(vessel_data_sets)

    # Load whale points
    whales_interp = gpd.read_file('data/whales/whales_coast.gpkg')

    # Load encounters
    vessel_encounters = gpd.read_file('data/encounters/vessel_encounters.gpkg')
    # whale_encounters = gpd.read_file('data/encounters/whale_encounters.gpkg')

    # Set up timestamps
    ranges = [('2020-07-01', '2023-06-01'),
              ('2020-07-01', '2021-06-01'),
              ('2021-07-01', '2022-06-01'),
              ('2022-07-01', '2023-06-01')]
    range_labels = ['full', '2020', '2021', '2022']

    timestamps = {
        label: pd.date_range(start, end, freq='30min')
        for label, (start, end) in zip(range_labels, ranges)
    }

    whale_mask = {
        label: (whales_interp.timestamp >= start) & (whales_interp.timestamp < end)
        for label, (start, end) in zip(range_labels, ranges)
    }

    vessel_mask = {
        label: (vessel_points.timestamp >= start) & (vessel_points.timestamp < end)
        for label, (start, end) in zip(range_labels, ranges)
    }

    encounters_mask = {
        label: (vessel_encounters.timestamp >= start) & (vessel_encounters.timestamp < end)
        for label, (start, end) in zip(range_labels, ranges)
    }

    # Line segment dfs (for plotting fading lines)
    # whales_segs = points_to_segments(whales_interp, 'id')
    # vessel_segs = points_to_segments(vessel_points, 'callsign')

    # Convert all gdfs to lat/long
    whales_interp = whales_interp.to_crs(4326)
    vessel_points = vessel_points.to_crs(4326)
    protected_areas = protected_areas.to_crs(4326)
    # whales_segs = whales_segs.to_crs(4326)
    # vessel_segs = vessel_segs.to_crs(4326)
    basemap = basemap.to_crs(4326)
    bounds_full = vessel_points.geometry.total_bounds
    bounds_ant = [178, -50.5, 179.5, -47]
    bounds_auck = [165.25, -51.5, 167.25, -49.5]
    bounds_camp = [168.25, -53, 169.75, -52]

    whales_interp.geometry = whales_interp.geometry.apply(fix_dateline)
    vessel_points.geometry = vessel_points.geometry.apply(fix_dateline)
    protected_areas.geometry = protected_areas.geometry.apply(fix_dateline)
    # whales_segs.geometry = whales_segs.geometry.apply(fix_dateline)
    # vessel_segs.geometry = vessel_segs.geometry.apply(fix_dateline)
    basemap.geometry = basemap.geometry.apply(fix_dateline)

    vessel_encounters = vessel_encounters.to_crs(4326)
    # whale_encounters = whale_encounters.to_crs(4326)

    vessel_encounters.geometry = vessel_encounters.geometry.apply(fix_dateline)
    # whale_encounters.geometry = whale_encounters.geometry.apply(fix_dateline)

    basemap_src = GeoJSONDataSource(geojson=basemap.to_json(default=str))

    bds = {
        'full': bounds_full,
        'auck': bounds_auck,
        'camp': bounds_camp,
        'anti': bounds_ant,
    }

    # Generate frames
    folder = '/pvol/static_maps'
    os.makedirs(folder, exist_ok=True)

    done = 0
    for bname in ['full', 'auck', 'camp', 'anti']:
        for label, tss in timestamps.items():
            for use_encounters in [True, False]:
                if use_encounters:
                    fname = os.path.join(folder, f'map_{bname}_{label}_enc.png')
                else:
                    fname = os.path.join(folder, f'map_{bname}_{label}.png')

                if path.isfile(fname):
                    continue

                print(bname, label, use_encounters)
                with timer('Full frame'):
                    with timer('Generate frame'):
                        fig = heatmap.animation_frame(
                            whales_interp[whale_mask[label]], vessel_points[vessel_mask[label]],
                            protected_areas, basemap_src, bds[bname],
                            vessel_encounters[encounters_mask[label]] if use_encounters else None)
                    with timer('Export frame'):
                        export_png(fig, filename=fname)

                    # Close selenium drivers to prevent memory bloat
                    webdriver_control.cleanup()
