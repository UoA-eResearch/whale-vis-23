import glob,shutil,tqdm

for f in tqdm.tqdm(glob.glob('frames/*.png')):
    shutil.move(f, '/pvol/frames/')

