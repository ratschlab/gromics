import sys
import numpy as np
import os
import h5py
import pickle
import re

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <annotation_gtf>\n' % sys.argv[0])
    sys.exit(1)
infile = sys.argv[1]

CONF = 2

def get_tags_gtf(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags

### collect junction information from GTF
trans2gene = dict()
transcripts = []
chrms = []
exons = []
strands = []
for l, line in enumerate(open(infile, 'r')):
    if l > 0 and l % 100000 == 0:
        sys.stdout.write('.')
        if l % 1000000 == 0:
            sys.stdout.write('%i\n' % l)
        sys.stdout.flush()
    if line[0] == '#':
        continue
    sl = line.strip().split('\t')
    if sl[2].lower() != 'exon':
        continue
    tags = get_tags_gtf(sl[8])
    trans2gene[tags['transcript_id']] = tags['gene_id']
    chrms.append(re.sub(r'chr', '', sl[0]))
    transcripts.append(tags['transcript_id'])
    exons.append([int(sl[3]) - 1, int(sl[4])]) 
    strands.append(sl[6])
transcripts = np.array(transcripts)
chrms = np.array(chrms)
exons = np.array(exons)
strands = np.array(strands)

sidx = np.argsort(transcripts)
transcripts = transcripts[sidx]
exons = exons[sidx, :]
strands = strands[sidx]
chrms = chrms[sidx]

junctions = []
jstrands = []
jchrs = []
jgenes = []
_,fidx = np.unique(transcripts, return_index=True)
lidx = np.r_[fidx[1:], transcripts.shape[0]]
for i in range(fidx.shape[0]):
    tidx = np.arange(fidx[i], lidx[i])
    eidx = np.argsort(exons[tidx, 0])
    if eidx.shape[0] > 1:
        junctions.append(np.reshape(exons[tidx[eidx], :].ravel()[1:-1], (eidx.shape[0] - 1, 2)))
        jchrs.extend(chrms[tidx[:-1]])
        jstrands.extend(strands[tidx[:-1]])
        jgenes.extend([trans2gene[x] for x in transcripts[tidx[:-1]]])
junctions = np.vstack(junctions)
jgenes = np.array(jgenes)
jchrs = np.array(jchrs)
jstrands = np.array(jstrands)

del exons, strands, transcripts

jid = np.array(['.'.join(x) for x in np.c_[jchrs.astype('str'), junctions.astype('str')]])
_, uidx = np.unique(jid, return_index=True)
junctions = junctions[uidx, :]
jstrands = jstrands[uidx]
jchrs = jchrs[uidx]
jgenes = jgenes[uidx]

### sort everything by coordinates
sidx = np.argsort(junctions[:, 0])
junctions = junctions[sidx, :]
jstrands = jstrands[sidx]
jchrs = jchrs[sidx]
jgenes = jgenes[sidx]
sidx = np.argsort(jchrs, kind='mergesort')
junctions = junctions[sidx, :]
jstrands = jstrands[sidx]
jchrs = jchrs[sidx]
jgenes = jgenes[sidx]

### prepare output file
OUT = h5py.File(re.sub(r'.gtf$', '', infile) + '.junctions.hdf5', 'w') 
OUT.create_dataset(name='gene_names', data=jgenes.view(np.chararray).encode('utf-8'), compression='gzip')
OUT.create_dataset(name='strand', data=jstrands.view(np.chararray).encode('utf-8'), compression='gzip')
OUT.create_dataset(name='pos', data=junctions, compression='gzip')
OUT.create_dataset(name='chrms', data=jchrs.view(np.chararray).encode('utf-8'), compression='gzip')
OUT.close()
