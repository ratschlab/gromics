import sys
import os
import numpy as np
import glob
import h5py
import re
import pickle

if len(sys.argv) < 4:
    sys.stderr.write('Usage: %s <pattern_of_sparse_alignment_hdfs> <threshold> <outbasename>\n' % sys.argv[0])
    sys.exit(1)
hdf5_pattern = sys.argv[1]
thresh = int(sys.argv[2])
outfname_junctions = sys.argv[3] + '.t%s' % sys.argv[2] + '.junction_map.pickle'

flist = glob.glob(hdf5_pattern)
### collect all possible junction coordinates
junctions = dict()
junction_list = dict()
all_junctions = 0
all_junctions_last = 1
for i,fname in enumerate(flist):
    sys.stderr.write('processing %i/%i - curr junctions found: %i (+ %.2f%%)\n' % (i + 1, len(flist), all_junctions, (float(all_junctions) / all_junctions_last * 100) - 100))
    all_junctions_last = max(1, all_junctions)
    with h5py.File(fname, 'r') as IN:
        for key in IN:
            if key.endswith('_introns_m'):
                strand = 2
            elif key.endswith('_introns_p'):
                strand = 1
            else:
                continue

            chrm = key.split('_')[0]
            if not (chrm, strand) in junction_list:
                junction_list[(chrm, strand)] = set()
                junctions[(chrm, strand)] = dict()
            kk_set = set([(_[0], _[1]) for _ in IN[key][:] if _[2] >= thresh])
            kk_new = kk_set - set(junctions[(chrm, strand)])
            for kk in kk_new:
                junctions[(chrm, strand)][kk] = len(junction_list[(chrm, strand)])
                junction_list[(chrm, strand)].add(kk)
                all_junctions += 1

### rewrite indices to be sorted by position rather than by incoming order
sys.stderr.write('re-ordering coordinate dict\n')
for key in junctions:
    if len(junctions[key]) == 0:
        continue
    pos = np.array([[k[0], k[1]] for k in junctions[key]])
    sidx = np.argsort(pos[:, 1])
    ssidx = np.argsort(pos[sidx, 0], kind='mergesort')
    sidx = sidx[ssidx]
    pos = pos[sidx, :]
    for i,p in enumerate(pos):
        junctions[key][(p[0], p[1])] = i

sys.stderr.write('writing results to %s\n' % outfname_junctions)
pickle.dump(junctions, open(outfname_junctions, 'wb'), -1)

