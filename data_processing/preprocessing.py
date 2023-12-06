from functools import partial
from os import path
import geopandas as gpd
from tqdm.contrib.concurrent import process_map
from pathlib import Path

from data_processing.load_data import load_vessel_points, load_basemap, reduce_poly_res, load_vessel_traces
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


def snap_points_to_coast(input_file, output_file):
    """Load all vessel data, snap to coast, and save pre-processed data"""
    assert path.isfile(input_file), f'Input file {input_file} does not exist'
    assert not path.isfile(output_file), f'Output file {output_file} already exists'

    # Vessel points, interpolated to 10 minutes, for calculating encounters
    points_df = load_vessel_points(input_file, 2193)

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
    points_df['geometry'] = process_map(fn, points_df['geometry'], chunksize=10000)

    points_df.to_file(output_file, driver='GPKG')


if __name__ == '__main__':
    folder = Path('data') / 'vessels'

    vessel_data_files = [(folder / f'{vessel_type}_points.gpkg', folder / f'{vessel_type}_points_coast.gpkg') for
                         vessel_type in ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']]

    for input_file, output_file in vessel_data_files:
        print(f'Processing {input_file}')
        snap_points_to_coast(input_file, output_file)
