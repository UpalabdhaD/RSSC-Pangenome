#! usr/bin/bash

SCRIPT="<PROJECT_ROOT>/scripts/ANI_clustermap.py"

pw_ani_path=$(realpath $1)
threshold=$2
cluster_member_out=$(realpath $3)
filtered_network=$(realpath $4)
metadata=$(realpath $5)
merged_out=$(realpath $6)

python $SCRIPT \
    "$pw_ani_path" \
    "$threshold" \
    "$cluster_member_out" \
    "$filtered_network" \
    "$metadata" \
    "$merged" \
    "$merged_out"