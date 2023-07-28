import pandas as pd
import geopandas as gpd


def load_whales(file_name, bounds, crs):
    # TODO: clip whale paths to vessel extents
    # TODO: join whale paths by ID
    wdf = pd.read_csv(file_name)

    # Drop unwanted columns
    wdf = wdf[['id', 'name', 'date', 'lat', 'lon']]

    # Convert to geodataframe
    wgdf = gpd.GeoDataFrame(wdf, geometry=gpd.points_from_xy(wdf['lon'], wdf['lat']), crs=4326)

    # Reproject
    wgdf = wgdf.to_crs(crs)

    # Clip to bounds
    wgdf = wgdf.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]

    return wgdf


def load_vessels(file_name, crs):
    vessels = gpd.read_file(file_name)

    # Some vessels have incorrect CRS set
    vessels = vessels.set_crs(4326, allow_override=True)

    # Reproject
    vessels = vessels.to_crs(crs)

    return vessels, vessels.geometry.total_bounds


def load_basemap(file_name, crs):
    basemap = gpd.read_file(file_name)

    # Drop non-land areas
    basemap = basemap[basemap['TA2022_V1_00'] != '999']

    # Reproject
    basemap = basemap.to_crs(crs)

    return basemap


def load_all(crs=2193):
    vessels, bounds = load_vessels('data/vessels/fishing_all.gpkg', crs=crs)

    whales = load_whales('data/whales/df_all_3.csv', bounds, crs=crs)

    basemap = load_basemap('data/territorial-authority-2022-generalised.gpkg', crs=crs)

    # TODO: Load MPAs

    return whales, vessels, basemap, bounds
