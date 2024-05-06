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

rule vessel_segments:
    input:
        expand('data/vessels/{vtype}_points_coast.gpkg', vtype=VTYPES)
    output:
        temp('data/intermediate/vessel_segments_{start}_{end}.gpkg')
    params:
        grouper='callsign'
    script:
        'snakemake/points_to_segments.py'

rule whale_segments:
    input:
        'data/whales/whales_coast.gpkg'
    output:
        temp('data/intermediate/whale_segments_{start}_{end}.gpkg')
    params:
        grouper='id'
    script:
        'snakemake/points_to_segments.py'

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