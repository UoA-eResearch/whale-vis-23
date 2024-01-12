import os
from functools import partial
from os import path
import geopandas as gpd
import pandas as pd
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from pathlib import Path

from data_processing import encounters
from data_processing.load_data import load_vessel_points, load_basemap, reduce_poly_res, load_vessel_traces, load_whales
from data_processing.snaptocoast import internal_points_to_coast


def clean_vessel_data(file_name):
    """Procedure for cleaning and simplifying the input file 'data/vessels/AIS_3_29_filtered.gpkg'"""
    gdf = gpd.read_file(file_name)

    # Drop unwanted columns and rename for clarity
    initial_columns = ['AIS_timest', 'Vessel_nam', 'Call_sign',
                       # 'Length__me', 'Width__met', 'Orientatio', 'Speed__kno',  # May be interesting later
                       'new_type', 'geometry']
    rename_columns = {'AIS_timest': 'timestamp',
                      # 'Length__me': 'length', 'Width__met': 'width',
                      # 'Orientation': 'heading', 'Speed__kno': 'speed',  # May be interesting later
                      'Vessel_nam': 'name', 'Call_sign': 'callsign', 'new_type': 'type'}
    gdf = gdf[initial_columns].rename(columns=rename_columns)

    # Some ships don't have call-signs, fill these in with truncated names
    gdf.loc[gdf['callsign'].isna(), 'callsign'] = gdf.loc[gdf['callsign'].isna(), 'name'].str[:7]
    gdf.loc[gdf['callsign'] == ' ', 'callsign'] = gdf.loc[gdf['callsign'] == ' ', 'name'].str[:7]
    gdf.drop(columns=['name'], inplace=True)

    # Sort by callsign and timestamp
    gdf = (
        gdf.sort_values(['callsign', 'timestamp'])
        .reset_index(drop=True)
    )

    # Split by type and save the result
    for vessel_type, group in gdf.groupby('type'):
        group.reset_index(drop=True).to_file(f'data/vessels/{vessel_type}_points.gpkg', driver='GPKG')


def snap_points_to_coast(input_gdf, output_file):
    """Load all vessel data, snap to coast, and save pre-processed data"""
    assert not path.isfile(output_file), f'Output file {output_file} already exists'
    # Load coastline
    _, bounds = load_vessel_traces('data/vessels/fishing_all.gpkg', crs=2193)
    basemap = load_basemap('data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.gpkg',
                           crs=2193, bounds=bounds).explode(index_parts=True)
    basemap = reduce_poly_res(basemap, 10).reset_index(drop=True)
    # Drop smaller polygons (these won't be visible, and this significantly speeds up processing time)
    basemap = basemap[basemap.area > 30000]

    # Snap points to coast
    fn = partial(internal_points_to_coast, coasts=basemap)
    # points_df = points_df.reset_index(drop=True)
    input_gdf['geometry'] = process_map(fn, input_gdf['geometry'], chunksize=10000)

    input_gdf.to_file(output_file, driver='GPKG')


def gen_encounters(source_df, target_df, output_file):
    """Generate a dataframe of encounters between vessels and whales"""
    encounters_df = source_df.copy()
    encounters_df['encounter_dist'] = encounters.dist_to_whales(source_df, target_df)

    encounters_df.to_file(output_file, driver='GPKG')


if __name__ == '__main__':
    folder = Path('data') / 'vessels'

    vessel_data_files = [(folder / f'{vessel_type}_points.gpkg', folder / f'{vessel_type}_points_coast.gpkg') for
                         vessel_type in ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]

    # Snap vessel data to coast
    for input_file, output_file in vessel_data_files:
        if path.isfile(output_file):
            continue
        print(f'Processing {input_file}')
        points_df = load_vessel_points(input_file, 2193)
        snap_points_to_coast(points_df, output_file)

    # Snap whale data to coast
    whale_coast_file = 'data/whales/whales_coast.gpkg'
    if not path.isfile(whale_coast_file):
        print('Processing whales file')
        _, bounds = load_vessel_traces('data/vessels/fishing_all.gpkg', crs=2193)
        whales_interp = load_whales('data/whales/df_all_3.csv', bounds, 2193, 10)
        snap_points_to_coast(whales_interp, whale_coast_file)

    # Generate encounters
    print('Generating encounters')
    vessel_data = [gpd.read_file(coast_file) for _, coast_file in tqdm(vessel_data_files, desc='Loading vessel data')]
    vessel_points = pd.concat(vessel_data)
    whale_points = gpd.read_file(whale_coast_file)

    os.makedirs('data/encounters', exist_ok=True)
    if not path.isfile('data/encounters/vessel_encounters.gpkg'):
        gen_encounters(vessel_points, whale_points, 'data/encounters/vessel_encounters.gpkg')
    if not path.isfile('data/encounters/whale_encounters.gpkg'):
        gen_encounters(whale_points, vessel_points, 'data/encounters/whale_encounters.gpkg')
