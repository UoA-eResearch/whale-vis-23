import json
import pandas as pd
from datetime import datetime

# Given filtered timestamps and date ranges, calculate the frame numbers from which to generate videos
timestamps = pd.read_parquet(snakemake.input[0]).index

# Take yearly intervals from start date
start = datetime.strptime(snakemake.wildcards.start, '%Y-%m-%d')
end = datetime.strptime(snakemake.wildcards.end, '%Y-%m-%d')
years = range(start.year, end.year + 1)

# Calculate frame numbers
output = {}
for year in years:
    start_year = datetime(year, start.month, start.day)
    end_year = datetime(year + 1, start.month, start.day)
    start_frame = int(sum(timestamps < start_year))
    end_frame = int(sum(timestamps < end_year))
    if start_frame != end_frame:  # Don't save empty years
        output[str(year)] = [start_frame, end_frame]

# Save frame numbers
with open(snakemake.output[0], 'w') as f:
    json.dump(output, f, indent=4)
