#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/GTEx/rna/results/

THRESH=0

python collect_junction_db_step1_coordinates.py "${basedir}/alignments/*.conf_2.filt.hdf5" ${THRESH} ${basedir}/gtex_junctions  
