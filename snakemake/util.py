import geopandas as gpd
import pandas as pd


def load_files(input_files: list[str]) -> gpd.GeoDataFrame:
    """Load either a single file or multiple files into a single GeoDataFrame."""
    if isinstance(input_files, str):
        data = gpd.read_file(input_files)
    elif len(input_files) == 1:
        data = gpd.read_file(input_files[0])
    else:
        data_sets = [gpd.read_file(fname) for fname in input_files]
        data = pd.concat(data_sets)

    return data
