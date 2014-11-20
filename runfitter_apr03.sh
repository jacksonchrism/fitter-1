#!/bin/bash
JOB_OUTPUT_FILE="output-$(basename $0)-$(date '+%s').txt"
source /data/snoplus/home/kate/install/rat-aging-linear/env_rat-aging-linear.sh
srun --comment="linear apr03" linearaging_test apr03 >logs/$JOB_OUTPUT_FILE 2>&1 &
echo "OUTPUT GOING TO: $JOB_OUTPUT_FILE"

