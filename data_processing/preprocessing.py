from functools import partial
from os import path
import pandas as pd
import geopandas as gpd
import swifter

from data_processing.load_data import load_vessel_points, load_basemap
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


def snap_vessels_to_coast():
    """Load all vessel data, snap to coast, and save pre-processed data"""
    # Vessel points, interpolated to 10 minutes, for calculating encounters
    vessel_data_files = [path.join('data', 'vessels', f'{vessel_type}_points.gpkg') for vessel_type in
                         ['Fishing', 'Passenger', 'Cargo', 'Tanker', 'Other']]
    vessel_data_sets = [load_vessel_points(vdf, 2193) for vdf in vessel_data_files]

    vessel_points = pd.concat(vessel_data_sets)

    # Snap points to coast
    basemap = load_basemap('data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.gpkg',
                                     crs=2193)

    fn = partial(internal_points_to_coast, coasts=basemap)
    vessel_points = vessel_points['geometry'].swifter.apply(fn)

    vessel_points.to_file(path.join('data', 'vessels', 'vessel_points_coast.gpkg'), driver='GPKG')
