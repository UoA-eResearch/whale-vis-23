import pandas as pd
import tqdm

tqdm.tqdm.pandas()


def _nearest_whale(row, whales):
    matches = whales[whales['timestamp'] == row['timestamp']]
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

    # Calculate distance to nearest whale
    out[mask] = vessels[mask].progress_apply(_nearest_whale, whales=whales, axis=1)

    return out
