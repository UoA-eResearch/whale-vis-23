from util import load_files

from data_processing.preprocessing import gen_encounters

# Load data
vessel_points = load_files(snakemake.input.vessels)
whale_points = load_files(snakemake.input.whales)

# Mask by date
start = snakemake.wildcards.start
end = snakemake.wildcards.end

vessel_points = vessel_points[(vessel_points.timestamp >= start) & (vessel_points.timestamp < end)]
whale_points = whale_points[(whale_points.timestamp >= start) & (whale_points.timestamp < end)]

# Generate encounters from {source} to {target}
if snakemake.wildcards.source == 'vessel':
    gen_encounters(vessel_points, whale_points, snakemake.output[0])
else:
    gen_encounters(whale_points, vessel_points, snakemake.output[0])
