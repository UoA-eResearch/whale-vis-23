from plotting.util import points_to_segments
from util import load_files

# Load data
data = load_files(snakemake.input)

# Mask to time range
start = snakemake.wildcards.start
end = snakemake.wildcards.end

mask = (data.timestamp >= start) & (data.timestamp < end)

grouper = snakemake.params.grouper

# Convert from lines to segments
data_segments = points_to_segments(data[mask], grouper)
data_segments.to_file(snakemake.output[0])
