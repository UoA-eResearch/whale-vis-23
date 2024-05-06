import geopandas as gpd


gpd.read_file(snakemake.input[0]).to_parquet(snakemake.output[0], engine='pyarrow')
