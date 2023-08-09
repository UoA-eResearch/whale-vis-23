import pandas as pd
import geopandas as gpd
import topojson as tp

from data_processing.interpolation import interpolate_trace


def load_whales(file_name, bounds, crs, interpolate_mins=None):
    # TODO: clip whale paths to vessel extents
    # TODO: join whale paths by ID
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

    # Reproject
    wgdf = wgdf.to_crs(crs)

    # Clip to bounds
    wgdf = wgdf.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]

    return wgdf


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


def load_vessel_traces(file_name, crs):
    vessels = gpd.read_file(file_name)

    # Some vessels have incorrect CRS set
    vessels = vessels.set_crs(4326, allow_override=True)

    # Reproject
    vessels = vessels.to_crs(crs)

    return vessels, vessels.geometry.total_bounds


def load_vessel_points(filename, crs):
    gdf = gpd.read_file(filename)

    vessel_type = gdf.loc[0, 'type']

    # Interpolate data individually for each vessel
    results = {}
    for callsign, group in gdf.groupby('callsign'):
        res = interpolate_trace(group)
        res['callsign'] = callsign
        res['type'] = vessel_type
        results['callsign'] = res

    # Combine into a single dataframe
    gdf = pd.concat(results, ignore_index=True)
    return gdf.to_crs(crs)


def load_protected_areas(crs):
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

    return gpd.GeoDataFrame(pd.concat([imma, mpa], ignore_index=True))


def load_basemap(file_name, crs):
    basemap = gpd.read_file(file_name)

    # Reproject
    basemap = basemap.to_crs(crs)

    return basemap


def reducy_poly_res(gdf, tolerance):
    """Reduce polygon resolution to reduce file size/plotting time"""
    topo = tp.Topology(gdf, prequantize=False)
    return topo.toposimplify(tolerance).to_gdf()


def load_all(crs=2193):
    vessels, bounds = load_vessel_traces('data/vessels/fishing_all.gpkg', crs=crs)

    whales = load_whales('data/whales/df_all_3.csv', bounds, crs=crs)

    protected_areas = load_protected_areas(crs=crs)

    # Coastlines from linz https://data.linz.govt.nz/layer/51153-nz-coastlines-and-islands-polygons-topo-150k/
    basemap = load_basemap('data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.shp', crs=crs)

    # Simplify baselayer topologies
    basemap = reducy_poly_res(basemap, 10)
    protected_areas = reducy_poly_res(protected_areas, 10)

    return whales, vessels, protected_areas, basemap, bounds
