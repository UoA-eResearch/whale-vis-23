VTYPES = ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']

def all_files(start, end):
    return [
        f'data/intermediate/vessel_pts_final_{start}_{end}.parquet',
        f'data/intermediate/whale_pts_final_{start}_{end}.parquet',
        f'data/intermediate/vessel_seg_final_{start}_{end}.parquet',
        f'data/intermediate/whale_seg_final_{start}_{end}.parquet',
        'data/intermediate/protected_final.parquet',
        'data/intermediate/basemap_final.parquet',
    ]

rule all:
    # Full files required by frame generation
    input:
        *all_files('2020-07-01', '2023-06-30')

rule test:
    # Small output files used by testing code
    input:
        *all_files('2020-08-01', '2020-08-30')

# Vessel data pre-processing: initial cleaning, interpolation, snap to coast, and segment
rule clean_vessel_data:
    input:
        'data/vessels/AIS_3_29_filtered.gpkg'
    output:
        expand('data/vessels/{vtype}_points.gpkg', vtype=VTYPES)
    script:
        'snakemake/clean_vessel_data.py'

rule vessel_points_to_coast:
    input:
        'data/vessels/{vtype}_points.gpkg'
    output:
        'data/vessels/{vtype}_points_coast.gpkg'
    params:
        whale=False
    script:
        'snakemake/snap_points_to_coast.py'

rule vessel_segments:
    input:
        expand('data/vessels/{vtype}_points_coast.gpkg', vtype=VTYPES)
    output:
        temp('data/intermediate/vessel_segments_{start}_{end}.gpkg')
    params:
        grouper='callsign'
    script:
        'snakemake/points_to_segments.py'

# Whale data pre-processing: snap to coast and segment
rule whale_points_to_coast:
    input:
        'data/whales/df_all_3.csv',
        'data/vessels/fishing_all.gpkg'  # Required for bounds
    output:
        'data/whales/whales_coast.gpkg'
    params:
        whale=True
    script:
        'snakemake/snap_points_to_coast.py'

rule whale_segments:
    input:
        'data/whales/whales_coast.gpkg'
    output:
        temp('data/intermediate/whale_segments_{start}_{end}.gpkg')
    params:
        grouper='id'
    script:
        'snakemake/points_to_segments.py'

# Final data pre-processing: mask dates, convert crs, fix dateline, convert to parquet format
rule vessel_pts_final:
    input:
        expand('data/vessels/{vtype}_points_coast.gpkg', vtype=VTYPES),
    output:
        temp('data/intermediate/vessel_pts_final_{start}_{end}.gpkg',)
    script:
        'snakemake/masked_data_final.py'

rule whale_pts_final:
    input:
        'data/whales/whales_coast.gpkg',
    output:
        temp('data/intermediate/whale_pts_final_{start}_{end}.gpkg',)
    script:
        'snakemake/masked_data_final.py'

rule vessel_seg_final:
    input:
        'data/intermediate/vessel_segments_{start}_{end}.gpkg',
    output:
        temp('data/intermediate/vessel_seg_final_{start}_{end}.gpkg')
    script:
        'snakemake/masked_data_final.py'

rule whale_seg_final:
    input:
        'data/intermediate/whale_segments_{start}_{end}.gpkg'
    output:
        temp('data/intermediate/whale_seg_final_{start}_{end}.gpkg')
    script:
        'snakemake/masked_data_final.py'

rule unmasked_data_to_crs:
    input:
        imma = 'data/imma_hr.gpkg',
        mpa = 'data/mpa3851.gpkg',
        basemap = 'data/linz_coastlines/nz-coastlines-and-islands-polygons-topo-150k.gpkg',
    output:
        protected_areas = temp('data/intermediate/protected_final.gpkg'),
        basemap = temp('data/intermediate/basemap_final.gpkg')
    script:
        'snakemake/unmasked_data_final.py'

rule to_parquet:
    input:
        'data/intermediate/{file}.gpkg',
    output:
        'data/intermediate/{file}.parquet',
    script:
        'snakemake/to_parquet.py'