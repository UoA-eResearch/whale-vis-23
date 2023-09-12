import pandas as pd
from shapely import LineString
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
