import pandas as pd
from shapely import LineString, Point, Polygon
from tqdm import tqdm


def points_to_segments(gdf, grouper):
    """Converts a GeoDataFrame containing points to one where each row is a line segment between two points"""
    # Sort by timestamp
    gdf = gdf.sort_values('timestamp')

    results = []
    for _, group in tqdm(desc='Converting vessel points to segments', iterable=gdf.groupby(grouper)):
        res = group.iloc[:-1].copy()
        res['geometry'] = list(map(LineString, zip(group['geometry'][:-1],
                                                   group['geometry'].shift(-1)[:-1])))
        results.append(res)

    out = pd.concat(results, ignore_index=True)

    # Convert timestamp to seconds since epoch (for easier comparison later)
    out['timestamp'] = out['timestamp'].astype('int64') // 10 ** 9

    return out


def fix_dateline(geom):
    """Apply to a gdf column to fix points that cross the dateline"""
    if isinstance(geom, LineString):
        return LineString([(x + 360, y) if x < 0 else (x, y) for x, y in geom.coords])
    elif isinstance(geom, Point):
        return Point((geom.x + 360, geom.y)) if geom.x < 0 else geom
    elif isinstance(geom, Polygon):
        return Polygon([(x + 360, y) if x < 0 else (x, y) for x, y in geom.exterior.coords])
    else:
        raise ValueError(f'Unexpected geometry type: {type(geom)}')