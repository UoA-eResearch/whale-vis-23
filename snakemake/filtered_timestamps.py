import json

import geopandas as gpd
import numpy as np
import pandas as pd

# Full time range
start = snakemake.wildcards.start
end = snakemake.wildcards.end
timestamps = pd.date_range(start, end, freq='30min')

# Load bounds and filter whales
bounds_id = snakemake.wildcards.bounds
with open(snakemake.input.bounds) as f:
    all_bounds = json.load(f)

assert bounds_id in all_bounds, f'Bounds {bounds_id} not found in {all_bounds.keys()}'

whales_interp = gpd.read_parquet(snakemake.input.whales)
select_bounds = all_bounds[bounds_id]
whales_interp = whales_interp.cx[select_bounds[0]:select_bounds[2], select_bounds[1]:select_bounds[3]]

# Filter every other timestamp if no whales present
interval = int(snakemake.wildcards.interval)
has_whale = timestamps.map(lambda ts: any(whales_interp['timestamp'] == ts))
in_interval = np.arange(len(timestamps)) % interval == 0
timestamps = timestamps[has_whale | in_interval]

# Convert to dataframe and save
pd.DataFrame(index=timestamps).to_parquet(snakemake.output[0])
