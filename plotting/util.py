import pandas as pd
from shapely import LineString, Point, Polygon, MultiPolygon
from tqdm import tqdm


def points_to_segments(gdf, grouper):
    """Converts a GeoDataFrame containing points to one where each row is a line segment between two points"""
    assert gdf.crs == 2193, "Expecting projected CRS"
    # Sort by timestamp
    gdf = gdf.sort_values('timestamp')

    results = []
    for _, group in tqdm(desc='Converting vessel points to segments', iterable=gdf.groupby(grouper)):
        res = group.iloc[1:].copy()
        res['geometry'] = list(map(LineString, zip(group['geometry'][:-1],
                                                   group['geometry'].shift(-1)[:-1])))

        # Drop any rows where the segment jumps an unreasonable distance
        res['_len'] = res['geometry'].length
        res = res.drop(res[res['_len'] > 10000].index)
        res.drop(columns=['_len'], inplace=True)

        results.append(res)

    out = pd.concat(results, ignore_index=True)

    return out


def fix_dateline(geom):
    """Apply to a gdf column to fix points that cross the dateline"""
    if isinstance(geom, LineString):
        return LineString([(x + 360, y) if x < 0 else (x, y) for x, y in geom.coords])
    elif isinstance(geom, Point):
        return Point((geom.x + 360, geom.y)) if geom.x < 0 else geom
    elif isinstance(geom, Polygon):
        return Polygon([(x + 360, y) if x < 0 else (x, y) for x, y in geom.exterior.coords])
    elif isinstance(geom, MultiPolygon):
        return MultiPolygon([fix_dateline(p) for p in geom.geoms])
    else:
        raise ValueError(f'Unexpected geometry type: {type(geom)}')
