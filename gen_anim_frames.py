import os
from os import path
import pandas as pd

from bokeh.io import export_png

from data_processing import load_data
from plotting import heatmap

if __name__ == '__main__':
    # Load base data
    whales, vessels, protected_areas, basemap, bounds = load_data.load_all()

    # Load vessel points
    vessel_data_files = [path.join('data', 'vessels', f'{vessel_type}_points.gpkg') for vessel_type in ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]
    vessel_data_sets = [load_data.load_vessel_points(vdf, 2193) for vdf in vessel_data_files]

    vessel_points = pd.concat(vessel_data_sets)

    # Load whale points
    whales_interp = load_data.load_whales('data/whales/df_all_3.csv', vessel_points.geometry.total_bounds, 2193, 10)

    # Set up timestamps
    start = '2022-08-01'
    end = '2022-09-01'
    timestamps = pd.date_range(start, end, freq='30min')

    whale_mask = (whales_interp.timestamp >= start) & (whales_interp.timestamp < end)
    vessel_mask = (vessel_points.timestamp >= start) & (vessel_points.timestamp < end)

    # Generate frames
    os.makedirs('frames', exist_ok=True)
    for i, ts in enumerate(timestamps):
        print(i)
        fig = heatmap.animation_frame(whales_interp[whale_mask], vessel_points[vessel_mask],
                                      protected_areas, basemap, bounds, ts)
        export_png(fig, filename=f'frames/frame_{i:04d}.png')
