# Whale & Vessel Visualisation

## Snakemake
Only the un-processed input files are included with this repository. Pre-processing is required to generate all the input files used in the generation of maps/videos.
The first step is to run:
```bash
python data_processing/preprocessing.py
```

These files can be generated using the Snakemake pipeline provided in the `snakemake` directory, using the command:
```bash
snakemake all --cores {n}
```
where `{n}` is the number of cores available on your machine.

## Generating video files
Videos are generated one frame at a time, by the script `gen_anim_frames.py`. This takes as arguments the bounds of the area to be visualised, and whether to visualise encounters.
Other parameters, such as the date range, frequency, and output directory are set in the script itself.

Individual plots are generated using the bokeh library, which requires a selenium webdriver in order to save output figures.

The script `gen_anim_frames.py` will only generate 600 frames at a time, as there appears to be a slight memory leak in the use of selenium. It therefore needs to be run multiple times in order to generate all required frames. This can be done by using the shell scripts `run.sh` and `run_enc.sh`.

Once all frames are generated, ffmpeg can be used to assemble these into videos. The notebook `video_output.ipynb` will generate shell scripts with the necessary ffmpeg commands for creating each individual output video.

## Generating static maps
Static maps are generated using the script `gen_static_maps.py`. This takes no arguments, with date ranges and output directories set in the script itself. Running this script once will generate all expected maps in the output directory.