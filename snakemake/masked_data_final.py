import json

from plotting.util import fix_dateline
from util import load_files

# Load data
data = load_files(snakemake.input.data)

# Mask to date range
start = snakemake.wildcards.start
end = snakemake.wildcards.end
mask = (data.timestamp >= start) & (data.timestamp < end)

# Convert crs and fix dateline
output = data[mask].to_crs(4326)
output.geometry = output.geometry.apply(fix_dateline)

# Mask to bounds
bounds_id = snakemake.wildcards.bounds
with open(snakemake.input.bounds) as f:
    all_bounds = json.load(f)

assert bounds_id in all_bounds, f'Bounds {bounds_id} not found in {all_bounds.keys()}'

select_bounds = all_bounds[bounds_id]
output = output.cx[select_bounds[0]:select_bounds[2], select_bounds[1]:select_bounds[3]]

# Save result
output.to_file(snakemake.output[0])