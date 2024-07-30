import json

# Generate a shell script will call run.sh/gen_anim_frames in order to generate all animation frames
shell_script = ['#!/usr/bin/env sh', '']

frame_dir = snakemake.params.frame_dir
shell_script.append(f'mkdir -p {frame_dir}')
shell_script.append(f'mkdir -p {frame_dir}_enc')
shell_script.append('mkdir -p ./logs')
shell_script.append('')

for bounds_name, fname in snakemake.input.items():
    with open(fname) as f:
        bounds = json.load(f)

    # Calculate total number of frames
    start_frame = min([frame_nums[0] for frame_nums in bounds.values()])
    end_frame = max([frame_nums[1] for frame_nums in bounds.values()])
    frame_no = end_frame - start_frame + 1
    # gen_anim_frames generates 600 frames at a time
    #   however, batches occasionally crash due to issues with selenium/chromedriver
    #   so, run an excessive (2x) number of batches to ensure all frames are generated
    #   batches with no work to do will still take ~5 seconds to load data,
    #   but this is negligible relative to the ~50 minutes for 600 frames
    batch_no = (frame_no // 600) * 2

    shell_script.append(f'# {bounds_name}')
    shell_script.append(f'sh run.sh {bounds_name} {batch_no} > ./logs/{bounds_name}.log &')  # Run in background - 1 batch per bounds
    shell_script.append(f'sh run.sh {bounds_name} {batch_no} --encounters > ./logs/{bounds_name}_enc.log &')
    shell_script.append('')

with open(snakemake.output[0], 'w', newline='') as f:
    f.write('\n'.join(shell_script))
