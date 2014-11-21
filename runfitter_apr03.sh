#!/bin/bash
JOB_OUTPUT_FILE="output-$(basename $0)-$(date '+%s').txt"
source /data/snoplus/home/jackson/SNO+/snoing/snoing/install/env_rat-aging-linear.sh
srun --comment="linear apr03" linearaging apr03 >logs/$JOB_OUTPUT_FILE 2>&1 &
echo "OUTPUT GOING TO: $JOB_OUTPUT_FILE"

