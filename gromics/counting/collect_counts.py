import sys
import os
import scipy as sp
import glob 
import pdb
import gzip
import h5py
import re

from ..utils import hdf5

def parse_options(argv):

    """Parses options from the command line """

    from argparse import ArgumentParser

    parser = ArgumentParser(prog='collect_counts',
                            description='This script collects expression counts in TSV format and aggregates them in a single HDF5 file. Please note \
                                         that there are several ways to provide the list of input files (see options below). You must specifiy the input \
                                         files using one of these options. If you specify multiple options, they are parsed in the following precedence: \
                                         -i, -f, -p.\n')
    parser.add_argument('-i', '--input', dest='infiles_fnames', metavar='STR', nargs='+', help='list of expression count files in TSV format', default='-')
    parser.add_argument('-f', '--filelist', dest='infiles_flist', metavar='STR', help='text file listing expression count files', default='-')
    parser.add_argument('-p', '--pattern', dest='infiles_fpattern', metavar='STR', help='search pattern describing list of expression count files', default='-')
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='STR', help='name of output file (will be hdf5)', default='-', required=True)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='switch on verbose output [off]', default=False)
    
    return parser.parse_args(argv[1:])


def main():

    options = parse_options(sys.argv)

    if options.infiles_fnames != '-':
        files = options.infiles_fnames
    elif options.infiles_flist != '-':
        files = sp.loadtxt(options.infiles_flist, dtype='str', delimiter='\t')
    elif options.infiles_fpattern != '-':
        files = glob.glob(options.infiles_fpattern)
    else:
        sys.stderr.write('ERROR: You need to provide at least one form of input via -i, -f, or -p\n')
        return 1

    OUT = h5py.File(options.outfile, 'w')

    counts = []
    header = ['feature']
    labels = []
    for f, fname in enumerate(files):
        if options.verbose:
            print('(%i / %i) Loading %s' % (f + 1, len(files), fname), file=sys.stderr)
        data = sp.loadtxt(fname, dtype='str', delimiter='\t')

        #sid = re.sub(r'.tsv$', '', fname.split('/')[-1])
        assert data[0, 0] == 'gene_id', 'ERROR: data has no header!'
        sid = data[0, 1]

        if f == 0:
            OUT.create_dataset('sids', data=sp.array([sid]).view(sp.chararray).encode('utf-8'), dtype='|S128', chunks=True, compression='gzip', maxshape=(None,))
            OUT.create_dataset('gids', data=data[1:, 0].view(sp.chararray).encode('utf-8'))
            OUT.create_dataset('counts', data=data[1:, 1][:, sp.newaxis].astype('int'), chunks=True, compression='gzip', maxshape=(data.shape[0], None))
        else:
            assert(sp.all(OUT['gids'][:].view(sp.chararray).decode('utf-8') == data[1:, 0]))
            hdf5.append(OUT, sp.array([sid], dtype='|S128'), 'sids')
            hdf5.append(OUT, data[1:, 1].astype('int'), 'counts')
        del data
    OUT.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
