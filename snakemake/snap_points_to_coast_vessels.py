from data_processing.load_data import load_vessel_points, load_vessel_traces, load_whales
from data_processing.preprocessing import snap_points_to_coast

# Load data
points_df = load_vessel_points(snakemake.input[0], 2193)

# Snap to coast and save
snap_points_to_coast(points_df, snakemake.output[0])