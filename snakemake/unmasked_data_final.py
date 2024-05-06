from data_processing import load_data
from plotting.util import fix_dateline

# Load data
_, _, protected_areas, basemap, bounds = load_data.load_all()

# Convert crs and fix dateline
protected_areas = protected_areas.to_crs(4326)
basemap = basemap.to_crs(4326)

protected_areas.geometry = protected_areas.geometry.apply(fix_dateline)
basemap.geometry = basemap.geometry.apply(fix_dateline)

# Save result
protected_areas.to_file(snakemake.output.protected_areas)
basemap.to_file(snakemake.output.basemap)