#!/bin/bash
JOB_OUTPUT_FILE="output-$(basename $0)-$(date '+%s').txt"
source /data/snoplus/home/kate/install/rat-tracking-aging/env_rat-tracking-aging.sh
srun --comment="10xmay02 fit 2par" 10xmay02 >logs/$JOB_OUTPUT_FILE 2>&1 &
echo "OUTPUT GOING TO: $JOB_OUTPUT_FILE"

