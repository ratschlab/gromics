#!/bin/bash

set -e

threads=1
mem=20000
pmem=$(($mem / $threads))

THRESH=0

basedir=/cluster/work/grlab/projects/GTEx/rna/results/
junction_map=${basedir}/gtex_junctions.t${THRESH}.junction_map.pickle
outdir=${basedir}/splice_burden/projection
mkdir -p $outdir

for fname in $(ls -1 ${basedir}/alignments/*.conf_2.filt.hdf5)
do
    fbase=$(basename $fname)
    outfname=${outdir}/${fbase%.hdf5}.t${THRESH}.projected.hdf5
    logname=${outdir}/${fbase%.hdf5}.t${THRESH}.projected.log
    if [ ! -f ${outfname} ]
    then
        echo "python $(pwd)/collect_junction_db_step2_project_counts.py $junction_map $fname ${outfname}" | bsub -M ${mem} -J compr_gtex -We 2:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -n $threads -o $logname
    fi
done
