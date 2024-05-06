import geopandas as gpd
import pandas as pd

from plotting.util import fix_dateline

# Load data
if len(snakemake.input) == 1:
    data = gpd.read_file(snakemake.input[0])
else:
    data_sets = [gpd.read_file(fname) for fname in snakemake.input]
    data = pd.concat(data_sets)

# Mask to date range
start = snakemake.wildcards.start
end = snakemake.wildcards.end
mask = (data.timestamp >= start) & (data.timestamp < end)

# Convert crs and fix dateline
output = data[mask].to_crs(4326)
output.geometry = output.geometry.apply(fix_dateline)

# Save result
output.to_file(snakemake.output[0])