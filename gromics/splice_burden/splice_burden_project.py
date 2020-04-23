import sys
import os
import scipy as sp
import h5py
import pickle
import re
from spladder.classes import gene as cgene
from spladder.classes import splicegraph as csplicegraph
from spladder.classes import segmentgraph as csegmentgraph
sys.modules['modules.classes.gene'] = cgene
sys.modules['modules.classes.splicegraph'] = csplicegraph
sys.modules['modules.classes.segmentgraph'] = csegmentgraph


def build_index(IN): #f_idx, l_idx, pos):

    print('building index')
    index = dict()
    for key in IN:
        if not key.startswith('pos_'):
            continue
        if key.endswith('_m'):
            dbkey = (key.split('_')[1], 2)
        else:
            dbkey = (key.split('_')[1], 1)
        index[dbkey] = dict()
        for i, p in enumerate(IN[key][:]):
            index[dbkey][(p[0], p[1])] = (key, i)

    return index


def main():

    if len(sys.argv) < 5:
        sys.stderr.write('Usage: %s <junction_db.hdf5> <spladder.pickle> <spladder.count.hdf5> <outfile>\n' % sys.argv[0])
        sys.exit(1)
    junction_db = sys.argv[1]
    spladder_pickle = sys.argv[2]
    spladder_count = sys.argv[3]
    outfname = sys.argv[4]

    STRANDS = ['.', '+', '-']
    CONF = 2

    ### given an edge, return the cooridnates of the edge (chrm, strand, start, stop)
    def compute_position(e):
        
        gidx = gene_ids_edges_spladder[e] 
        eidx = edge_index_spladder[e]

        ### get index pair of segments
        if genes[gidx].segmentgraph.segments.size == 0:
            genes[gidx].segmentgraph = csegmentgraph.Segmentgraph(genes[gidx])
        a,b = sp.unravel_index(eidx, genes[gidx].segmentgraph.seg_edges.shape)
        
        ### get INTRON coordinates
        start = genes[gidx].segmentgraph.segments[1, a]
        stop = genes[gidx].segmentgraph.segments[0, b]

        return (genes[gidx].chr, genes[gidx].strand, start, stop)


    ### define a mapping between edges in the foreground set and the outgroup
    def compute_projection(e):

        ### get coordinates
        chrm, strand, start, stop = compute_position(e)

        ### check in outgroup introns
        try:
            return junction_index[(chrm, STRANDS.index(strand))][(start, stop)]
        except KeyError:
            try:
                return junction_index[(chrm, 0)][(start, stop)]
            except KeyError:
                return -1

    ### preprocessing to harmonize genes used in the analysis
    print('loading sample count data from %s' % spladder_count)
    IN_TC = h5py.File(spladder_count, 'r')
    gene_ids_edges_spladder = IN_TC['gene_ids_edges'][:, 0]
    assert sp.all(gene_ids_edges_spladder == sp.sort(gene_ids_edges_spladder))
    edge_index_spladder = IN_TC['edge_idx'][:].astype('int')

    ### load spladder event data
    print('loading sample genes from %s' % spladder_pickle)
    genes, features = pickle.load(open(spladder_pickle, 'rb'), encoding='latin1')

    ### load outgroup junction db file
    print('loading outgroup junction db from %s' % junction_db)
    IN_GT = h5py.File(junction_db, 'r')

    ### index junctions for faster search
    junction_index = build_index(IN_GT) #f_idx, l_idx, pos_gt)

    ### generate sample junctions file
    print('generating sample junction file')
    OUT = h5py.File(outfname, 'w')
    OUT.create_dataset(name='gene_ids', data=gene_ids_edges_spladder, compression='gzip')
    OUT.create_dataset(name='gene_names', data=IN_TC['gene_names'][:], compression='gzip')
    OUT.create_dataset(name='strains', data=IN_TC['strains'][:], compression='gzip')
    OUT.create_dataset(name='strand', data=sp.array([genes[_].strand for _ in gene_ids_edges_spladder]).view(sp.chararray).encode('utf-8'), compression='gzip')
    OUT.create_dataset(name='chrms', data=sp.array([genes[_].chr for _ in gene_ids_edges_spladder]).view(sp.chararray).encode('utf-8'), compression='gzip')

    ### compute projection
    print('computing coordinate projection')
    project = []
    positions = []
    for i in range(edge_index_spladder.shape[0]):
        if i > 0 and i % 1000 == 0:
            sys.stderr.write('.')
            if i % 10000 == 0:
                sys.stderr.write('%i\n' % i)
            sys.stderr.flush()
        project.append(compute_projection(i))
        positions.append(compute_position(i)[2:])
    project_dict = dict()
    for i,p in enumerate(project):
        if p == -1:
            continue
        try:
            project_dict[p[0]][p[1]].append(i)
        except KeyError:
            try:
                project_dict[p[0]][p[1]] = [i]
            except KeyError:
                project_dict[p[0]] = dict()
                project_dict[p[0]][p[1]] = [i]

    ### generate empty edge object
    print('creating empty output set')
    out_edges = sp.zeros((IN_TC['edges'].shape[0], IN_GT['strains'].shape[0]), dtype='int')
    ### fill in values iterating over chromosomes
    for key in IN_GT.keys():
        if not key.startswith('junctions_'):
            continue
        sys.stderr.write('processing %s\n' % key)
        pkey = re.sub(r'^junctions_', 'pos_', key)
        if pkey not in project_dict:
            sys.stderr.write('\nWARNING: %s not found - skipping\n' % pkey)
            continue
        cs = 10*IN_GT[key].chunks[0]
        for i in range(0, IN_GT[key].shape[0], cs):
            print('%i/%i' % (i, IN_GT[key].shape[0]))
            tmp = IN_GT[key][i:i+cs, :]
            for j in range(tmp.shape[0]):
                if not i+j in project_dict[pkey]:
                    continue
                for k in project_dict[pkey][i+j]:
                    out_edges[k, :] = tmp[j, :]

    OUT.create_dataset(name='pos', data=sp.array(positions), compression='gzip')
    OUT.create_dataset(name='edges_outgroup', data=out_edges, compression='gzip')
    OUT.create_dataset(name='strains_outgroup', data=IN_GT['strains'][:], compression='gzip')
    OUT.close()
    IN_TC.close()
    IN_GT.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
