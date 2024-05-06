import os
from os import path

import geopandas as gpd
import pandas as pd
from bokeh.io import export_png
from bokeh.io.webdriver import webdriver_control
from bokeh.models import GeoJSONDataSource

from plotting import heatmap
from utils import timer

if __name__ == '__main__':
    # Generates static maps for all bounds, time ranges, with or without encounters
    # Load pre-processed data
    start = '2020-07-01'
    end = '2023-06-30'
    whales_interp = gpd.read_parquet(f'data/intermediate/whale_pts_final_{start}_{end}.parquet')
    vessel_points = gpd.read_parquet(f'data/intermediate/vessel_pts_final_{start}_{end}.parquet')
    protected_areas = gpd.read_parquet('data/intermediate/protected_final.parquet')
    whales_segs = gpd.read_parquet(f'data/intermediate/whale_seg_final_{start}_{end}.parquet')
    vessel_segs = gpd.read_parquet(f'data/intermediate/vessel_seg_final_{start}_{end}.parquet')
    basemap = gpd.read_parquet('data/intermediate/basemap_final.parquet')
    vessel_encounters = gpd.read_parquet(f'data/intermediate/vessel_encounters_final_{start}_{end}.parquet')

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

    basemap_src = GeoJSONDataSource(geojson=basemap.to_json(default=str))

    # Map bounds
    bounds_full = vessel_points.geometry.total_bounds
    bounds_ant = [178, -50.5, 179.5, -47]
    bounds_auck = [165.25, -51.5, 167.25, -49.5]
    bounds_camp = [168.25, -53, 169.75, -52]

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
                            encounters=vessel_encounters[encounters_mask[label]] if use_encounters else None)
                    with timer('Export frame'):
                        export_png(fig, filename=fname)

                    # Close selenium drivers to prevent memory bloat
                    webdriver_control.cleanup()
