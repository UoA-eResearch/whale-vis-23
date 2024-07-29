import json

import geopandas as gpd
import numpy as np
import pandas as pd

# Full time range
start = snakemake.wildcards.start
end = snakemake.wildcards.end
timestamps = pd.date_range(start, end, freq='30min')

# Filter every other timestamp if no whales present
interval = int(snakemake.wildcards.interval)
whales_interp = gpd.read_parquet(snakemake.input[0])
has_whale = timestamps.map(lambda ts: any(whales_interp['timestamp'] == ts))
in_interval = np.arange(len(timestamps)) % interval == 0
timestamps = timestamps[has_whale | in_interval]

# Convert to dataframe and save
pd.DataFrame(index=timestamps).to_parquet(snakemake.output[0])
