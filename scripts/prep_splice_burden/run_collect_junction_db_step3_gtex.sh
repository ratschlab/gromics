#!/bin/bash

set -e


THRESH=0

basedir=/cluster/work/grlab/projects/GTEx/rna/results
junction_map=${basedir}/gtex_junctions.t${THRESH}.junction_map.pickle
logname=${junction_map%.junction_map.pickle}.lsf.log

mem=20000
threads=1

echo "python $(pwd)/collect_junction_db_step3_gather.py $junction_map ${basedir}/splice_burden/projection" | bsub -M ${mem} -J gathr_gtex -We 48:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -n $threads -o $logname

