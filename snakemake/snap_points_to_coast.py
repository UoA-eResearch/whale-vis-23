from data_processing.load_data import load_vessel_points, load_vessel_traces, load_whales
from data_processing.preprocessing import snap_points_to_coast

# Load data (whale or vessel)
# TODO: generate interpolated whale/vessel data and load, rather than interpolating here
if snakemake.params.whale:
    _, bounds = load_vessel_traces('data/vessels/fishing_all.gpkg', crs=2193)
    points_df = load_whales(snakemake.input[0], bounds, 2193, 10)
else:
    points_df = load_vessel_points(snakemake.input[0], 2193)

# Snap to coast and save
snap_points_to_coast(points_df, snakemake.output[0])