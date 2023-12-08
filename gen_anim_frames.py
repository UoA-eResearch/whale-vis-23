import os
from os import path
import pandas as pd
import geopandas as gpd

from bokeh.io import export_png

from data_processing import load_data
from plotting import heatmap
from plotting.util import fix_dateline

if __name__ == '__main__':
    # Load base data
    whales, vessels, protected_areas, basemap, bounds = load_data.load_all()

    # Load vessel points
    vessel_data_files = [path.join('data', 'vessels', f'{vessel_type}_points_coast.gpkg') for vessel_type in ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]
    vessel_data_sets = [load_data.load_vessel_points(vdf, 2193) for vdf in vessel_data_files]

    vessel_points = pd.concat(vessel_data_sets)

    # Load whale points
    whales_interp = gpd.read_file('data/whales/whales_coast.gpkg')

    # Set up timestamps
    start = '2022-08-01'
    end = '2022-09-01'
    timestamps = pd.date_range(start, end, freq='30min')

    whale_mask = (whales_interp.timestamp >= start) & (whales_interp.timestamp < end)
    vessel_mask = (vessel_points.timestamp >= start) & (vessel_points.timestamp < end)

    # Convert all gdfs to lat/long
    whales_interp = whales_interp.to_crs(4326)
    vessel_points = vessel_points.to_crs(4326)
    protected_areas = protected_areas.to_crs(4326)
    basemap = basemap.to_crs(4326)
    bounds = vessel_points.geometry.total_bounds

    whales_interp.geometry = whales_interp.geometry.apply(fix_dateline)
    vessel_points.geometry = vessel_points.geometry.apply(fix_dateline)
    protected_areas.geometry = protected_areas.geometry.apply(fix_dateline)
    basemap.geometry = basemap.geometry.apply(fix_dateline)

    # Generate frames
    os.makedirs('frames', exist_ok=True)
    for i, ts in enumerate(timestamps):
        print(i)
        fig = heatmap.animation_frame(whales_interp[whale_mask], vessel_points[vessel_mask],
                                      protected_areas, basemap, bounds, ts)
        export_png(fig, filename=f'frames/frame_{i:04d}.png')
