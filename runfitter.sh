#!/bin/bash
JOB_OUTPUT_FILE="output-$(basename $0)-$(date '+%s').txt"
source /data/snoplus/home/kate/install/rat-aging-const/env_rat-aging-const.sh
srun --comment="apr03" constaging_test apr03 >logs/$JOB_OUTPUT_FILE 2>&1 &
echo "OUTPUT GOING TO: $JOB_OUTPUT_FILE"

