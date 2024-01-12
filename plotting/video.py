from os import path
import os
import subprocess


def ffmpeg(frames_dir, frame_format, output_name, start_frame=0, end_frame=None, framerate=30, skip_frames=0,
           call_command=False):
    """
    Use ffmpeg to convert png files in frames/ to mp4 video
    """
    if end_frame is not None:
        # Use only frames between start_frame and end_frame
        frames_command = ['-frames:v', str(end_frame - start_frame + 1)]
    else:
        frames_command = []

    # Pad image size to even dimensions
    pad_command = 'pad=ceil(iw/2)*2:ceil(ih/2)*2'
    # pad_command = ["-vf", '"pad=ceil(iw/2)*2:ceil(ih/2)*2"']

    if skip_frames > 1:
        # Use only every nth frame
        skip_command = ["-vf", f"\"{pad_command},select='not(mod(n,{skip_frames}))',setpts=N/{framerate}/TB\""]
    else:
        skip_command = ["-vf", f"\"{pad_command}\""]

    # Call ffmpeg
    command = ['ffmpeg', '-framerate', str(framerate), '-start_number', str(start_frame),
               '-i', f'{frames_dir}/{frame_format}',
               *frames_command,
               '-c:v', 'libx264', '-profile:v', 'high', '-crf', '20', '-pix_fmt', 'yuv420p',
               *skip_command,
               f'{output_name}',
               '-y']

    print(' '.join(command))

    if call_command:
        print()
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            line = proc.stderr.readline()
            if not line:
                break
            print(line.decode('utf-8').strip())

    return ' '.join(command)