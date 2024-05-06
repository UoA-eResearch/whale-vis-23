import geopandas as gpd
import pandas as pd

from plotting.util import points_to_segments

# Load data
if len(snakemake.input) == 1:
    data = gpd.read_file(snakemake.input[0])
else:
    data_sets = [gpd.read_file(fname) for fname in snakemake.input]
    data = pd.concat(data_sets)

# Mask to time range
start = snakemake.wildcards.start
end = snakemake.wildcards.end

mask = (data.timestamp >= start) & (data.timestamp < end)

grouper = snakemake.params.grouper

# Convert from lines to segments
data_segments = points_to_segments(data[mask], grouper)
data_segments.to_file(snakemake.output[0])
