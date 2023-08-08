import pandas as pd
import geopandas as gpd
import topojson as tp


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
    vessels, bounds = load_vessels('data/vessels/fishing_all.gpkg', crs=crs)

    whales = load_whales('data/whales/df_all_3.csv', bounds, crs=crs)

    protected_areas = load_protected_areas(crs=crs)

    # Coastlines from linz https://data.linz.govt.nz/layer/51153-nz-coastlines-and-islands-polygons-topo-150k/
    basemap = load_basemap('data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.shp', crs=crs)

    # Simplify baselayer topologies
    basemap = reducy_poly_res(basemap, 10)
    protected_areas = reducy_poly_res(protected_areas, 10)

    return whales, vessels, protected_areas, basemap, bounds
