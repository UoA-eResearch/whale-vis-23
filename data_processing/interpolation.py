import datetime

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.interpolate import interp1d


def ceil_time(dt, delta):
    """Round up datetime to nearest delta"""
    return dt + (datetime.datetime.min - dt) % delta


def interpolate_trace(df, interval=600):
    """Resample trace to be at regular intervals (e.g.: 10 minutes)"""
    # Extract x/y coordinates, convert timestamps to seconds since start
    input = df.copy().reset_index()
    input['x'] = input['geometry'].x
    input['y'] = input['geometry'].y
    input['seconds'] = (input['timestamp'] - input.loc[0, 'timestamp']).dt.total_seconds()

    interp = interp1d(input['seconds'], input[['x', 'y']].to_numpy().T, fill_value=np.nan)

    # Create arrays of seconds since start, and datetimes for new sampling interval
    # (start time is rounded up to nearest interval, stop time is rounded down)
    start_time = ceil_time(input.loc[0, 'timestamp'].tz_localize(None), datetime.timedelta(seconds=interval))
    start_seconds = (start_time - input.loc[0, 'timestamp'].tz_localize(None)).seconds

    out_seconds = np.arange(start_seconds, input['seconds'].max(), interval, dtype=float)
    out_dt = [input.loc[0, 'timestamp'] + datetime.timedelta(seconds=s) for s in out_seconds]

    # Build output dataframe
    out = pd.DataFrame(interp(out_seconds).T).rename(columns={0: 'x', 1: 'y'})
    out['timestamp'] = out_dt

    out_gdf = gpd.GeoDataFrame(out, geometry=gpd.points_from_xy(out['x'], out['y']), crs=input.crs)
    out_gdf.drop(columns=['x', 'y'], inplace=True)

    return out_gdf
