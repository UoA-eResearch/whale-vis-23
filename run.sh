#!/usr/bin/env sh
# run the script as: ./run.sh bname
# where bname is auck, camp, full, anti

# Re-run the gen_anim_frames script multiple times, as it gets slower over time
echo "$1"
for i in `seq 1 $2`
do
    # run the command multiple times (to refresh memory)
    if [ "$3" = "--encounters" ]; then
        python -m gen_anim_frames $1 --encounters
    else
        python -m gen_anim_frames $1
    fi
done