from plotting.util import fix_dateline
from util import load_files

# Load data
data = load_files(snakemake.input)

# Mask to date range
start = snakemake.wildcards.start
end = snakemake.wildcards.end
mask = (data.timestamp >= start) & (data.timestamp < end)

# Convert crs and fix dateline
output = data[mask].to_crs(4326)
output.geometry = output.geometry.apply(fix_dateline)

# Save result
output.to_file(snakemake.output[0])