import json

import geopandas as gpd
import numpy as np
import pandas as pd

# Full time range
start = snakemake.wildcards.start
end = snakemake.wildcards.end
timestamps = pd.date_range(start, end, freq='30min')

# Filter every nth timestamp if no whales present
interval = int(snakemake.wildcards.interval)
no_vessel_interval = int(snakemake.params.no_vessel_interval)

if no_vessel_interval is None or no_vessel_interval < interval:
    # If no_vessel_interval is not set or is less than interval, vessel free animation should not be slower than
    # that with vessels
    no_vessel_interval = interval

# Show all frames with a whale present
whales_interp = gpd.read_parquet(snakemake.input.whales)
has_whale = timestamps.map(lambda ts: any(whales_interp['timestamp'] == ts))

# Show some frames if a vessel is present
vessels_interp = gpd.read_parquet(snakemake.input.vessels)
has_vessel = timestamps.map(lambda ts: any(vessels_interp['timestamp'] == ts))

in_interval = np.arange(len(timestamps)) % interval == 0
in_interval = in_interval & has_vessel

# Show even fewer frames if no vessel is present
in_no_vessel_interval = np.arange(len(timestamps)) % no_vessel_interval == 0

timestamps = timestamps[has_whale | in_interval | in_no_vessel_interval]

# Convert to dataframe and save
pd.DataFrame(index=timestamps).to_parquet(snakemake.output[0])
