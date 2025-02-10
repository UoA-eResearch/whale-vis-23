#!/usr/bin/env sh
ffmpeg -framerate 30 -start_number 17520 -i /pvol/frames/frame_full_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/full_2021.mp4 -y

ffmpeg -framerate 30 -start_number 0 -i /pvol/frames/frame_auck_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/auck_2020.mp4 -y

ffmpeg -framerate 30 -start_number 17520 -i /pvol/frames/frame_auck_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/auck_2021.mp4 -y

ffmpeg -framerate 30 -start_number 0 -i /pvol/frames/frame_camp_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/camp_2020.mp4 -y

ffmpeg -framerate 30 -start_number 17520 -i /pvol/frames/frame_camp_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/camp_2021.mp4 -y

ffmpeg -framerate 30 -start_number 35040 -i /pvol/frames/frame_auck_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/auck_2022.mp4 -y

ffmpeg -framerate 30 -start_number 35040 -i /pvol/frames/frame_camp_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/camp_2022.mp4 -y

ffmpeg -framerate 30 -start_number 35040 -i /pvol/frames/frame_full_%04d.png -frames:v 17520 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output/full_2022.mp4 -y

