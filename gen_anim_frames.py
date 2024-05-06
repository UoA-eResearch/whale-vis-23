import os
from os import path
import pandas as pd
import geopandas as gpd

from bokeh.io import export_png
from bokeh.io.webdriver import webdriver_control
from bokeh.models import GeoJSONDataSource

from plotting import heatmap
from utils import timer
import argparse


if __name__ == '__main__':
    # Get bounds name from argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('bounds', type=str, choices=['full', 'auck', 'camp', 'anti'])
    parser.add_argument('--encounters', action='store_true', help='Plot encounters')
    args = parser.parse_args()
    bname = args.bounds
    use_encounters = args.encounters

    start = '2020-07-01'
    end = '2023-06-30'
    timestamps = pd.date_range(start, end, freq='30min')

    whales_interp = gpd.read_file(f'data/intermediate/whale_pts_final_{start}_{end}.parquet')
    vessel_points = gpd.read_file(f'data/intermediate/vessel_pts_final_{start}_{end}.parquet')
    protected_areas = gpd.read_file('data/intermediate/protected_final.parquet')
    whales_segs = gpd.read_file(f'data/intermediate/whale_seg_final_{start}_{end}.parquet')
    vessel_segs = gpd.read_file(f'data/intermediate/vessel_seg_final_{start}_{end}.parquet')
    basemap = gpd.read_file('data/intermediate/basemap_final.parquet')


    if use_encounters:
        vessel_encounters = gpd.read_file(f'data/intermediate/vessel_encounters_final_{start}_{end}.parquet')
        # whale_encounters = gpd.read_file(f'data/intermediate/whale_encounters_final_{start}_{end}.parquet')
    else:
        vessel_encounters = None
        # whale_encounters = None

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
    }[bname]

    # Generate frames
    if not use_encounters:
        folder = '/pvol/frames'
    else:
        folder = '/pvol/frames_enc'
    os.makedirs(folder, exist_ok=True)

    done = 0
    for i, ts in enumerate(timestamps):
        fname = os.path.join(folder, f'frame_{bname}_{i:04d}.png')
        if path.isfile(fname):
            continue
        print(bname, i)
        with timer('Full frame'):
            with timer('Generate frame'):
                fig = heatmap.animation_frame_fade(
                    whales_segs, vessel_segs,
                    whales_interp, vessel_points,
                    protected_areas, basemap_src, bds, ts, vessel_encounters)
            with timer('Export frame'):
                export_png(fig, filename=fname)

            # Close selenium drivers to prevent memory bloat
            webdriver_control.cleanup()

        # Exit after generating some frames (selenium slows down after many frames)
        done += 1
        if done > 600:
            break
