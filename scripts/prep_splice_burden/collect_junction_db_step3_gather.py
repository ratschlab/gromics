import sys
import os
import pickle
import h5py
import numpy as np
import re
import glob

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <coordinate_db> <projected_dir>\n' % sys.argv[0])
    sys.exit(1)

junc_pickle = sys.argv[1]
projected_dir = sys.argv[2]
thresh = os.path.basename(junc_pickle).split('.')[-3]

print('junction coordinate db is: ' + junc_pickle)
outfname = re.sub(r'.junction_map.pickle$', '', junc_pickle) + '.hdf5'
print('writing junction count db to: ' + outfname)

### load coordinate DB
junctions = pickle.load(open(junc_pickle, 'rb'))

### collect files containing projected counts
flist = glob.glob(os.path.join(projected_dir, '*.%s.projected.hdf5' % thresh))

### gather single result files into joint DB
OUT = h5py.File(outfname, 'w')
samples = []
for i,fname in enumerate(flist):
    sys.stderr.write('processing %i/%i\n' % (i + 1, len(flist)))
    samples.append(os.path.basename(fname).split('.')[0])
    with h5py.File(fname, 'r') as IN:
        for key in IN:
            if IN[key].shape[0] == 0:
                continue
            print('.. key ' + key)
            if not key in OUT:
                OUT.create_dataset(name=key, data=IN[key][:], compression='gzip', maxshape=(IN[key].shape[0], None)) 
                ### add positions
                if key.endswith('_m'):
                    dbkey = (key.split('_')[1], 2)
                elif key.endswith('_p'):
                    dbkey = (key.split('_')[1], 1)
                pos = np.array([[_[0], _[1], junctions[dbkey][_]]for _ in junctions[dbkey]])
                sidx = np.argsort(pos[:, 2])
                pos = pos[sidx, :]
                OUT.create_dataset(name=re.sub(r'^junctions_', 'pos_', key), data=np.array(pos)[:, :2], compression='gzip')
            else:
                tmp = OUT[key].shape
                OUT[key].resize((tmp[0], tmp[1] + 1))
                OUT[key][:, tmp[1]] = IN[key][:, 0]
        ### generate entries for all contigs not present in the file
        for dbkey in junctions:
            if len(junctions[dbkey]) == 0:
                continue
            if dbkey[1] == 2:
                key = 'junctions_%s_m' % dbkey[0]
            else:
                key = 'junctions_%s_p' % dbkey[0]
            if not key in IN:
                print('...adding key %s' % key)
                if not key in OUT:
                    OUT.create_dataset(name=key, data=np.zeros((len(junctions[dbkey]), 1), dtype='int'), compression='gzip', maxshape=(len(junctions[dbkey]), None))
                    pos = np.array([[_[0], _[1], junctions[dbkey][_]]for _ in junctions[dbkey]])
                    sidx = np.argsort(pos[:, 2])
                    pos = pos[sidx, :]
                    OUT.create_dataset(name=re.sub(r'^junctions_', 'pos_', key), data=np.array(pos)[:, :2], compression='gzip')
                else:
                    tmp = OUT[key].shape
                    OUT[key].resize((tmp[0], tmp[1] + 1))
                    OUT[key][:, tmp[1]] = np.zeros((len(junctions[dbkey]),), dtype='int')
                    
OUT.create_dataset(name='samples', data=np.array(samples, dtype='str').view(np.chararray).encode('utf-8'), compression='gzip')
OUT.close()

