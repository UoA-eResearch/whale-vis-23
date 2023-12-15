import os
from os import path
import pandas as pd
import geopandas as gpd

from bokeh.io import export_png, webdriver
from bokeh.models import GeoJSONDataSource
from tqdm import tqdm

from data_processing import load_data
from plotting import heatmap
from plotting.util import fix_dateline
from utils import timer
import argparse


if __name__ == '__main__':
    # Get bounds name from argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('bounds', type=str, choices=['full', 'auck', 'camp', 'anti'])
    args = parser.parse_args()
    bname = args.bounds

    # Load base data
    whales, vessels, protected_areas, basemap, bounds = load_data.load_all()

    # Load vessel points
    vessel_data_files = [path.join('data', 'vessels', f'{vessel_type}_points_coast.gpkg') for vessel_type in ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]
    vessel_data_sets = [gpd.read_file(vdf) for vdf in tqdm(vessel_data_files, desc="Loading vessel data")]

    vessel_points = pd.concat(vessel_data_sets)

    # Load whale points
    whales_interp = gpd.read_file('data/whales/whales_coast.gpkg')

    # Set up timestamps
    # ranges = [('2020-07-01', '2023-06-01'),
    #           ('2020-07-01', '2021-06-01'),
    #           ('2021-07-01', '2022-06-01'),
    #           ('2022-07-01', '2023-06-01')]
    # range_labels = ['full', '2020', '2021', '2022']
    #
    # timestamps = {
    #     label: pd.date_range(start, end, freq='30min')
    #     for label, (start, end) in zip(range_labels, ranges)
    # }
    #
    # whale_mask = {
    #     label: (whales_interp.timestamp >= start) & (whales_interp.timestamp < end)
    #     for label, (start, end) in zip(range_labels, ranges)
    # }
    #
    # vessel_mask = {
    #     label: (vessel_points.timestamp >= start) & (vessel_points.timestamp < end)
    #     for label, (start, end) in zip(range_labels, ranges)
    # }

    start = '2020-07-01'
    end = '2023-06-01'
    timestamps = pd.date_range(start, end, freq='30min')
    whale_mask = (whales_interp.timestamp >= start) & (whales_interp.timestamp < end)
    vessel_mask = (vessel_points.timestamp >= start) & (vessel_points.timestamp < end)

    # Convert all gdfs to lat/long
    whales_interp = whales_interp.to_crs(4326)
    vessel_points = vessel_points.to_crs(4326)
    protected_areas = protected_areas.to_crs(4326)
    basemap = basemap.to_crs(4326)
    bounds_full = vessel_points.geometry.total_bounds
    bounds_ant = [178, -50.5, 179.5, -47]
    bounds_auck = [165.25, -51.5, 167.25, -49.5]
    bounds_camp = [168.25, -53, 169.75, -52]

    whales_interp.geometry = whales_interp.geometry.apply(fix_dateline)
    vessel_points.geometry = vessel_points.geometry.apply(fix_dateline)
    protected_areas.geometry = protected_areas.geometry.apply(fix_dateline)
    basemap.geometry = basemap.geometry.apply(fix_dateline)

    basemap_src = GeoJSONDataSource(geojson=basemap.to_json(default=str))

    bds = {
        'full': bounds_full,
        'auck': bounds_auck,
        'camp': bounds_camp,
        'anti': bounds_ant,
    }[bname]

    # Generate frames
    os.makedirs('frames', exist_ok=True)
    for i, ts in enumerate(timestamps):
        print(bname, i)
        with timer('Full frame'):
            fname = f'frames/frame_{bname}_{i:04d}.png'
            if path.isfile(fname):
                continue

            with timer('Generate frame'):
                fig = heatmap.animation_frame(whales_interp[whale_mask], vessel_points[vessel_mask],
                                              protected_areas, basemap_src, bds, ts)
            with timer('Export frame'):
                export_png(fig, filename=fname)

            # Close selenium tab to prevent memory bloat
            with timer('Cleanup'):
                driver = webdriver.webdriver_control.get()
                print(driver)
                driver.close()
