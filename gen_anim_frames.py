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
from plotting.util import fix_dateline, points_to_segments
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

    # Load base data
    whales, vessels, protected_areas, basemap, bounds = load_data.load_all()

    # Load vessel points
    vessel_data_files = [path.join('data', 'vessels', f'{vessel_type}_points_coast.gpkg') for vessel_type in ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]
    vessel_data_sets = [gpd.read_file(vdf) for vdf in tqdm(vessel_data_files, desc="Loading vessel data")]

    vessel_points = pd.concat(vessel_data_sets)

    # Load whale points
    whales_interp = gpd.read_file('data/whales/whales_coast.gpkg')

    # Load encounters
    if use_encounters:
        vessel_encounters = gpd.read_file('data/encounters/vessel_encounters.gpkg')
        # whale_encounters = gpd.read_file('data/encounters/whale_encounters.gpkg')
    else:
        vessel_encounters = None

    start = '2020-07-01'
    end = '2023-06-30'
    timestamps = pd.date_range(start, end, freq='30min')
    whale_mask = (whales_interp.timestamp >= start) & (whales_interp.timestamp < end)
    vessel_mask = (vessel_points.timestamp >= start) & (vessel_points.timestamp < end)

    # Line segment dfs (for plotting fading lines)
    whales_segs = points_to_segments(whales_interp[whale_mask], 'id')
    vessel_segs = points_to_segments(vessel_points[vessel_mask], 'callsign')

    # Convert all gdfs to lat/long
    whales_interp = whales_interp.to_crs(4326)
    vessel_points = vessel_points.to_crs(4326)
    protected_areas = protected_areas.to_crs(4326)
    whales_segs = whales_segs.to_crs(4326)
    vessel_segs = vessel_segs.to_crs(4326)
    basemap = basemap.to_crs(4326)
    bounds_full = vessel_points.geometry.total_bounds
    bounds_ant = [178, -50.5, 179.5, -47]
    bounds_auck = [165.25, -51.5, 167.25, -49.5]
    bounds_camp = [168.25, -53, 169.75, -52]

    whales_interp.geometry = whales_interp.geometry.apply(fix_dateline)
    vessel_points.geometry = vessel_points.geometry.apply(fix_dateline)
    protected_areas.geometry = protected_areas.geometry.apply(fix_dateline)
    whales_segs.geometry = whales_segs.geometry.apply(fix_dateline)
    vessel_segs.geometry = vessel_segs.geometry.apply(fix_dateline)
    basemap.geometry = basemap.geometry.apply(fix_dateline)

    if use_encounters:
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
    }[bname]

    # Generate frames
    folder = 'frames'
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
                    whales_interp[whale_mask], vessel_points[vessel_mask],
                    protected_areas, basemap_src, bds, ts, vessel_encounters)
            with timer('Export frame'):
                export_png(fig, filename=fname)

            # Close selenium drivers to prevent memory bloat
            webdriver_control.cleanup()

        # Exit after generating some frames (selenium slows down after many frames)
        done += 1
        if done > 600:
            break
