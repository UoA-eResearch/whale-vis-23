import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
import topojson as tp
from shapely import LineString

from data_processing.interpolation import interpolate_trace, split_on_gaps


def load_whales(file_name, bounds, crs, interpolate_mins=None):
    """Load whale tracking points, optionally interpolate to finer time resolution"""
    wdf = pd.read_csv(file_name, parse_dates=['date'])

    # Drop unwanted columns
    wdf = wdf[['id', 'name', 'date', 'lat', 'lon']]
    wdf.rename(columns={'date': 'timestamp'}, inplace=True)

    # Convert to geodataframe
    wgdf = gpd.GeoDataFrame(wdf, geometry=gpd.points_from_xy(wdf['lon'], wdf['lat']), crs=4326)
    wgdf.drop(columns=['lat', 'lon'], inplace=True)

    if interpolate_mins is not None:
        # Resample data to finer time resolution and interpolate points
        # TODO: check this is valid in given crs
        wgdf = wgdf.sort_values(['id', 'timestamp'])
        results = {}
        for id, group in wgdf.groupby('id'):
            for grp, sub_group in split_on_gaps(group, 36*3600):
                if len(sub_group) > 0:
                    res = interpolate_trace(sub_group, interval=60*interpolate_mins)
                    res['id'] = f'{id}_{grp}'
                    res['name'] = sub_group.iloc[0]['name']

                    results[f'{id}_{grp}'] = res

        wgdf = pd.concat(results, ignore_index=True)

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

    return gpd.GeoDataFrame(results, crs=crs)


def load_vessel_traces(file_name, crs):
    vessels = (
        gpd.read_file(file_name)
        .set_crs(4326, allow_override=True)  # Some vessels have incorrect CRS set
        .to_crs(crs)                         # Reproject
        .drop(columns=['track_seg', 'begin', 'end', 'layer', 'path'])
    )

    return vessels, vessels.geometry.total_bounds


def load_vessel_points(filename, crs):
    gdf = gpd.read_file(filename)

    vessel_type = gdf.loc[0, 'type']

    # Interpolate data individually for each vessel
    results = {}
    for callsign, group in gdf.groupby('callsign'):
        for grp, sub_group in split_on_gaps(group):
            if len(sub_group) > 4:
                res = interpolate_trace(sub_group)
                res['callsign'] = f'{callsign}_{grp}'
                res['type'] = vessel_type
                results[f'{callsign}_{grp}'] = res

    # Combine into a single dataframe
    gdf = pd.concat(results, ignore_index=True)
    gdf['timestamp'] = gdf['timestamp'].dt.tz_localize(None)
    return gdf.to_crs(crs)


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


def load_basemap(file_name, crs, bounds=None):
    basemap = (
        gpd.read_file(file_name)
        .to_crs(crs)    # Reproject
        .drop(columns=['name', 'macronated', 'grp_macron', 'TARGET_FID', 'grp_ascii', 'grp_name', 'name_ascii'])
        .reset_index(drop=True)
    )

    if bounds is not None:
        # basemap = basemap.clip(bounds)   # Clip to bounding box
        # Rather than clipping to bounds, exclude polygons outside bounds, but keep intersecting polygons unaltered
        bbox = shapely.geometry.box(*bounds)
        basemap = basemap[basemap.intersects(bbox)]

    # Drop small polygons
    basemap = basemap[basemap.area > 1500]

    return basemap.explode(index_parts=False).reset_index()


def reduce_poly_res(gdf, tolerance):
    """Reduce polygon resolution to reduce file size/plotting time"""
    topo = tp.Topology(gdf, prequantize=False)
    return topo.toposimplify(tolerance).to_gdf()


def load_all(crs=2193):
    vessels, bounds = load_vessel_traces('data/vessels/fishing_all.gpkg', crs=crs)

    whales = load_whales('data/whales/df_all_3.csv', bounds, crs=crs)

    protected_areas = load_protected_areas(bounds, crs=crs)

    # Coastlines from linz https://data.linz.govt.nz/layer/51153-nz-coastlines-and-islands-polygons-topo-150k/
    basemap = load_basemap('data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.gpkg',
                           crs=crs, bounds=bounds)

    # Simplify vessel lines
    vessels['geometry'] = list(map(lambda x: x.simplify(100), vessels['geometry']))

    # Simplify baselayer topologies
    basemap = reduce_poly_res(basemap, 10)
    protected_areas = reduce_poly_res(protected_areas, 10)

    return whales, vessels, protected_areas, basemap, bounds
