#!/bin/bash

# List of filenames (without extensions)
FILENAMES=("phgg" "dipg" "atrt" "medullo")

# Loop through each filename in the list
for file_prefix in "${FILENAMES[@]}"
do
  echo "Processing file: $file_prefix"

  # Construct filenames for disruptive and nondisruptive cases
  FILE_NAME="${file_prefix}_top10mse"

  # Run the two commands in parallel
  taskset -c 0-16 python scripts/SuPreMo.py /pollard/data/projects/shzhang/cbtn/results/svs_tumortype_abc_FINAL/$FILE_NAME.txt \
    --file $FILE_NAME \
    --dir /pollard/data/projects/shzhang/cbtn/results/svs_tumortype_abc_FINAL/ \
    --get_Akita_scores \
    --get_maps \
    --get_tracks \
    --genome hg38 \
    --fa /pollard/home/shzhang/akita_variant_scoring/data/hg38.fa \
    --shift_by 0 #&

#   taskset -c 9-16 python scripts/SuPreMo.py /pollard/data/projects/shzhang/cbtn/svs_for_abc/$NONDISRUPTIVE_FILE \
#     --file $NONDISRUPTIVE_FILE \
#     --dir /pollard/data/projects/shzhang/cbtn/svs_for_abc/ \
#     --get_Akita_scores \
#     --get_maps \
#     --get_tracks \
#     --genome hg38 \
#     --fa /pollard/home/shzhang/akita_variant_scoring/data/hg38.fa \
#     --shift_by 0 &

#   # Wait for both commands to finish before processing the next file
#   wait
  echo "Completed processing for $file_prefix"
done

echo "All files processed."
