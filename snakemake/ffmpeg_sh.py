import json

from plotting import video

# Generate a shell script that will call ffmpeg to generate each video (per year, per bounds)
shell_script = ['#!/usr/bin/env sh', '']
frame_dir = snakemake.params.frame_dir

for bounds_name, fname in snakemake.input.items():
    with open(fname) as f:
        bounds = json.load(f)

    shell_script.append(f'# {bounds_name}')
    # ffmpeg command for each yearly video
    for year, frame_nums in bounds.items():
        command = video.ffmpeg(frame_dir, f'frame_{bounds_name}_%04d.png',
                               f'output/video_{bounds_name}_{year}.mp4',
                               start_frame=frame_nums[0], end_frame=frame_nums[1])
        shell_script.append(command)
        shell_script.append('')

    # ffmpeg command for full video
    start_frame = min([frame_nums[0] for frame_nums in bounds.values()])
    end_frame = max([frame_nums[1] for frame_nums in bounds.values()])

    command = video.ffmpeg(frame_dir, f'frame_{bounds_name}_%04d.png',
                           f'output/video_{bounds_name}_full.mp4',
                            start_frame=start_frame, end_frame=end_frame)

    shell_script.append(command)
    shell_script.append('')

with open(snakemake.output[0], 'w', newline='') as f:
    f.write('\n'.join(shell_script))
