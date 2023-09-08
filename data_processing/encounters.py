import pandas as pd
import tqdm

tqdm.tqdm.pandas()


def _nearest_whale(row, whales):
    # matches = whales[whales['timestamp'] == row['timestamp']]
    matches = whales.loc[row['timestamp']]
    if isinstance(matches, pd.Series):
        return matches['geometry'].distance(row['geometry'])
    dists = matches.distance(row['geometry'])
    # match = matches.loc[dists.idxmin(), 'id']
    return dists.min()


def dist_to_whales(vessels, whales):
    assert vessels.crs == whales.crs
    assert 'timestamp' in vessels.columns
    assert 'timestamp' in whales.columns

    # Mask vessels points that have a whale present
    mask = vessels['timestamp'].isin(whales['timestamp'])

    # Create NaN column for distance to nearest whale
    out = pd.Series(index=vessels.index, dtype=float)

    # Map whales by timestamp for more efficient inner loop
    whales_copy = whales.copy()
    whales_copy['timestamp'] = whales_copy['timestamp'].astype('int64') // 10**9
    vessels_copy = vessels.copy()
    vessels_copy['timestamp'] = vessels_copy['timestamp'].astype('int64') // 10**9
    whales_copy = whales_copy.set_index('timestamp')

    # Calculate distance to nearest whale
    out[mask] = vessels_copy[mask].progress_apply(_nearest_whale, whales=whales_copy, axis=1)

    return out
