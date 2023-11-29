import numpy as np
import pandas as pd
import geopandas as gpd
import topojson as tp
from joblib import Memory
from shapely import LineString

from data_processing.interpolation import interpolate_trace
from data_processing.snaptocoast import internal_points_to_coast

memory = Memory('cache_data/', verbose=1)


@memory.cache()
def load_whales(file_name, bounds, crs, interpolate_mins=None):
    """Load whale tracking points, optionally interpolate to finer time resolution"""
    wdf = pd.read_csv(file_name, parse_dates=['date'])

    # Drop unwanted columns
    wdf = wdf[['id', 'name', 'date', 'lat', 'lon']]

    if interpolate_mins is not None:
        # Resample data to finer time resolution and interpolate points
        # TODO: check this is valid in given crs
        wdf = (
            wdf.set_index('date')
            .sort_values('id')
            .groupby('id')
            .resample(f'{interpolate_mins}min')  # resample to finer time resolution
            .interpolate()                       # interpolate lat/long
            .ffill()                             # fill non-numeric cols (name, id)
            .reset_index(level='id', drop=True)  # drop id from index
            .reset_index()                       # move date index back to column
        )

    # Convert to geodataframe
    wgdf = gpd.GeoDataFrame(wdf, geometry=gpd.points_from_xy(wdf['lon'], wdf['lat']), crs=4326)
    wgdf.drop(columns=['lat', 'lon'], inplace=True)
    wgdf.rename(columns={'date': 'timestamp'}, inplace=True)

    # Reproject
    wgdf = wgdf.to_crs(crs)

    # Clip to bounds
    wgdf = wgdf.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]

    return wgdf


def load_whale_lines(filename, crs):
    """Loads whale tracking points and converts them to linestrings"""
    whale_points = load_whales(filename, [0, 0, np.inf, np.inf], crs)

    results = []
    for id, group in whale_points.groupby('id'):
        res = {'id': id, 'name': group.iloc[0]['name'],
               'geometry': LineString(group['geometry'].tolist())}
        results.append(res)

    return gpd.GeoDataFrame(results)


def _clean_vessel_data(file_name):
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


@memory.cache()
def load_vessel_traces(file_name, crs):
    vessels = (
        gpd.read_file(file_name)
        .set_crs(4326, allow_override=True)  # Some vessels have incorrect CRS set
        .to_crs(crs)                         # Reproject
        .drop(columns=['track_seg', 'begin', 'end', 'layer', 'path'])
    )

    return vessels, vessels.geometry.total_bounds


@memory.cache()
def load_vessel_points(filename, crs):
    gdf = gpd.read_file(filename)

    vessel_type = gdf.loc[0, 'type']

    # Interpolate data individually for each vessel
    results = {}
    for callsign, group in gdf.groupby('callsign'):
        res = interpolate_trace(group)
        res['callsign'] = callsign
        res['type'] = vessel_type
        results[callsign] = res

    # Combine into a single dataframe
    gdf = pd.concat(results, ignore_index=True)
    gdf['timestamp'] = gdf['timestamp'].dt.tz_localize(None)
    return gdf.to_crs(crs)


@memory.cache()
def load_protected_areas(bounds, crs):
    imma = (
        gpd.read_file('data/imma_hr.gpkg')  # Manually edited imma geometry to match NZ coastline
        [['Title', 'geometry']]
        .rename(columns={'Title': 'name'})
        .explode()
        .to_crs(crs)
    )

    mpa = (
        gpd.read_file('data/mpa3851.gpkg')
        [['Name', 'geometry']]
        .rename(columns={'Name': 'name'})
        .explode()
        .to_crs(crs)
    )

    imma['ptype'] = 'imma'
    mpa['ptype'] = 'mpa'

    protected_areas = (
        gpd.GeoDataFrame(pd.concat([imma, mpa], ignore_index=True))
        .clip(bounds)   # Clip to bounding box
        .reset_index(drop=True)
    )

    return protected_areas


@memory.cache()
def load_basemap(file_name, bounds, crs):
    basemap = (
        gpd.read_file(file_name)
        .to_crs(crs)    # Reproject
        .clip(bounds)   # Clip to bounding box
        .drop(columns=['name', 'macronated', 'grp_macron', 'TARGET_FID', 'grp_ascii', 'grp_name', 'name_ascii'])
        .reset_index(drop=True)
    )

    # Drop small polygons
    basemap = basemap[basemap.area > 1500]

    return basemap


def reducy_poly_res(gdf, tolerance):
    """Reduce polygon resolution to reduce file size/plotting time"""
    topo = tp.Topology(gdf, prequantize=False)
    return topo.toposimplify(tolerance).to_gdf()


@memory.cache()
def load_all(crs=2193):
    vessels, bounds = load_vessel_traces('data/vessels/fishing_all.gpkg', crs=crs)

    whales = load_whales('data/whales/df_all_3.csv', bounds, crs=crs)

    protected_areas = load_protected_areas(bounds, crs=crs)

    # Coastlines from linz https://data.linz.govt.nz/layer/51153-nz-coastlines-and-islands-polygons-topo-150k/
    basemap = load_basemap('data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.gpkg', bounds, crs=crs)

    # Simplify vessel lines
    vessels['geometry'] = list(map(lambda x: x.simplify(100), vessels['geometry']))

    # Simplify baselayer topologies
    basemap = reducy_poly_res(basemap, 10)
    protected_areas = reducy_poly_res(protected_areas, 10)

    # Snap trace points within land to nearest coast
    whales['geometry'] = whales['geometry'].progress_apply(internal_points_to_coast, coasts=basemap)

    return whales, vessels, protected_areas, basemap, bounds
