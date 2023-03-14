#!/bin/bash
#$ -o sgeout/
#$ -e sgeout/
#$ -cwd
#$ -q high@bronte,high@cerbero,high@chimera,high@gerione,high@nemeo,high@ortro

# Load modules
module load python3/3.8.5

pred_path=$1
# Define path to output directory
out_dir=$2
# Define debug level
debug_lvl=debug

# Define path to list of lists
# in_lists="${out_dir}/db_files/lists.dat"
in_lists="./lists.txt"
# Define path to input list
curr_list=$(sed "${SGE_TASK_ID}q;d" "${in_lists}")
# Loop through each input file
while read in_path; do
  # Debug
  echo "Current input path is ${in_path}"
  # Run Python script (substitute input with output)
  python3 ./parse_to_binary.py -i "${pred_path}" -db -o "${out_dir}" -ll "${debug_lvl}"
# Read current line
done < "${curr_list}"
