#!/bin/bash

#$ -o sgeout/
#$ -e sgeout/
#$ -cwd
#$ -q high@bronte,high@cerbero,high@chimera,high@gerione,high@nemeo,high@ortro

# Load modules
module load python3/3.8.5

db_path=$1
# Define path to output directory
out_dir=$2
# Define path to list
in_lists=$3
# Define path to input list
curr_file=$(sed "${SGE_TASK_ID}q;d" "$3")
# Define debug level
debug_lvl=info

python3 ./bin/parse_to_binary.py -i "${curr_file}" -db "${db_path}" -o "${out_dir}" "${debug_lvl}"
