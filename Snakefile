VTYPES = ['Fishing', 'Other', 'Cargo', 'Passenger', 'Tanker']
BOUNDS = ['full', 'auck', 'camp', 'anti']


def all_files(start, end, interval=3):
    timestamps = [f'data/timestamps/filtered_timestamps_{bound}_{start}_{end}_{interval}.parquet' for bound in BOUNDS]
    frame_nos = [f'data/timestamps/frame_numbers_{bound}_{start}_{end}_{interval}.json' for bound in BOUNDS]
    vessel_pts_final = [f'data/intermediate/vessel_pts_final_{bound}_{start}_{end}.parquet' for bound in BOUNDS]
    whale_pts_final = [f'data/intermediate/whale_pts_final_{bound}_{start}_{end}.parquet' for bound in BOUNDS]
    vessel_seg_final = [f'data/intermediate/vessel_seg_final_{bound}_{start}_{end}.parquet' for bound in BOUNDS]
    whale_seg_final = [f'data/intermediate/whale_seg_final_{bound}_{start}_{end}.parquet' for bound in BOUNDS]
    return [
        f'data/intermediate/vessel_pts_final_{start}_{end}.parquet',
        f'data/intermediate/whale_pts_final_{start}_{end}.parquet',
        f'data/intermediate/vessel_seg_final_{start}_{end}.parquet',
        f'data/intermediate/whale_seg_final_{start}_{end}.parquet',
        f'data/intermediate/vessel_encounters_final_{start}_{end}.parquet',
        f'data/intermediate/whale_encounters_final_{start}_{end}.parquet',
        'data/intermediate/protected_final.parquet',
        'data/intermediate/basemap_final.parquet',
        *timestamps,
        *frame_nos,
        *vessel_pts_final,
        *whale_pts_final,
        *vessel_seg_final,
        *whale_seg_final,
        f'frame_gen_{start}_{end}_{interval}.sh',
        f'ffmpeg_commands_{start}_{end}_{interval}.sh',
        f'ffmpeg_commands_{start}_{end}_{interval}_enc.sh'
    ]

wildcard_constraints:
    # date format: YYYY-MM-DD
    start = '[0-9]{4}-[0-9]{2}-[0-9]{2}',
    end = '[0-9]{4}-[0-9]{2}-[0-9]{2}',
    bounds = '[a-z]*',
    interval = '[0-9]*'

rule all:
    # Full files required by frame generation
    input:
        *all_files('2020-07-01', '2023-06-30')

rule test:
    # Small output files used by testing code
    input:
        *all_files('2020-08-01', '2020-08-30')

# Timestamps
rule filtered_timestamps:
    input:
        whales='data/intermediate/whale_pts_final_{bounds}_{start}_{end}.parquet',
        vessels='data/intermediate/vessel_pts_final_{bounds}_{start}_{end}.parquet'
    output:
        'data/timestamps/filtered_timestamps_{bounds}_{start}_{end}_{interval}.parquet'
    params:
        no_vessel_interval=20  # Show every nth frame if no vessels present. Set to None/0 for same speed as no whales
    script:
        'snakemake/filtered_timestamps.py'

# Given filtered timestamps and date ranges, calculate the frame numbers from which to generate videos
rule frame_numbers:
    input:
        'data/timestamps/filtered_timestamps_{bounds}_{start}_{end}_{interval}.parquet',
    output:
        'data/timestamps/frame_numbers_{bounds}_{start}_{end}_{interval}.json'
    script:
        'snakemake/frame_numbers.py'

def all_frame_numbers(bounds):
    return {
        bds: f'data/timestamps/frame_numbers_{bds}_{{start}}_{{end}}_{{interval}}.json'
        for bds in BOUNDS
    }

rule frame_gen_sh:
    input:
        **all_frame_numbers(BOUNDS)
    output:
        'frame_gen_{start}_{end}_{interval}.sh'
    params:
        frame_dir='/pvol/frames'
    script:
        'snakemake/frame_gen_sh.py'

# Generate ffmpeg commands for each video
rule ffmpeg_sh:
    input:
        **all_frame_numbers(BOUNDS)
    output:
        'ffmpeg_commands_{start}_{end}_{interval}.sh'
    params:
        frame_dir='/pvol/frames',
        output_dir='output/'
    script:
        'snakemake/ffmpeg_sh.py'

rule ffmpeg_sh_enc:
    input:
        **all_frame_numbers(BOUNDS)
    output:
        'ffmpeg_commands_{start}_{end}_{interval}_enc.sh'
    params:
        frame_dir='/pvol/frames_enc',
        output_dir='output_enc/'
    script:
        'snakemake/ffmpeg_sh.py'

# Vessel data pre-processing: initial cleaning, interpolation, snap to coast, and segment
rule clean_vessel_data:
    input:
        'data/vessels/AIS_3_29_filtered.gpkg'
    output:
        expand('data/vessels/{vtype}_points.gpkg', vtype=VTYPES)
    script:
        'snakemake/clean_vessel_data.py'

# rule vessel_points_to_coast:
#     input:
#         'data/vessels/{vtype}_points.gpkg'
#     output:
#         'data/vessels/{vtype}_points_coast.gpkg'
#     script:
#         'snakemake/snap_points_to_coast_vessels.py'

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
# rule whale_points_to_coast:
#     input:
#         'data/whales/df_all_3.csv',
#         'data/vessels/fishing_all.gpkg'  # Required for bounds
#     output:
#         'data/whales/whales_coast.gpkg'
#     script:
#         'snakemake/snap_points_to_coast_whales.py'

rule whale_segments:
    input:
        'data/whales/whales_coast.gpkg'
    output:
        temp('data/intermediate/whale_segments_{start}_{end}.gpkg')
    params:
        grouper='id'
    script:
        'snakemake/points_to_segments.py'

# Encounters
rule encounters:
    input:
        vessels=expand('data/vessels/{vtype}_points_coast.gpkg', vtype=VTYPES),
        whales='data/whales/whales_coast.gpkg',
    output:
        temp('data/intermediate/{source}_encounters_{start}_{end}.gpkg')
    script:
        'snakemake/encounters.py'

# Final data pre-processing: mask dates, convert crs, fix dateline, convert to parquet format
rule vessel_pts_final:
    input:
        data=expand('data/vessels/{vtype}_points_coast.gpkg', vtype=VTYPES),
        bounds='data/bounds.json'
    output:
        temp('data/intermediate/vessel_pts_final_{bounds}_{start}_{end}.gpkg',)
    script:
        'snakemake/masked_data_final.py'

rule whale_pts_final:
    input:
        data='data/whales/whales_coast.gpkg',
        bounds='data/bounds.json'
    output:
        temp('data/intermediate/whale_pts_final_{bounds}_{start}_{end}.gpkg',)
    script:
        'snakemake/masked_data_final.py'

rule vessel_seg_final:
    input:
        data='data/intermediate/vessel_segments_{start}_{end}.gpkg',
        bounds='data/bounds.json'
    output:
        temp('data/intermediate/vessel_seg_final_{bounds}_{start}_{end}.gpkg')
    script:
        'snakemake/masked_data_final.py'

rule whale_seg_final:
    input:
        data='data/intermediate/whale_segments_{start}_{end}.gpkg',
        bounds='data/bounds.json'
    output:
        temp('data/intermediate/whale_seg_final_{bounds}_{start}_{end}.gpkg')
    script:
        'snakemake/masked_data_final.py'

rule encounters_final:
    input:
        'data/intermediate/{source}_encounters_{start}_{end}.gpkg',
    output:
        temp('data/intermediate/{source}_encounters_final_{start}_{end}.gpkg')
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