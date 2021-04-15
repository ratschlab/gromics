import sys
import os
import numpy as np
import glob 
import pdb
import gzip
import h5py
import re

def parse_options(argv):

    """Parses options from the command line """

    from argparse import ArgumentParser

    parser = ArgumentParser(prog='count_hdf2tsv',
                            description='This script takes an aggregated expression count file in hdf5 format and converts it into tsv.')
    parser.add_argument('-i', '--infile', dest='infile', metavar='STR', help='input file name', required=True)
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='STR', help='output file name', required=True)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='switch on verbose output [off]', default=False)
    
    return parser.parse_args(argv[1:])


def main():

    options = parse_options(sys.argv)

    if options.verbose:
        print('Writing to %s' % options.outfile)
        print('Reading data from %s' % options.infile)

    IN = h5py.File(options.infile, 'r')
    cs = IN['counts'].chunks[0] 

    with open(options.outfile, 'w') as OUT:
        genes = IN['gids'][:].view(np.chararray).decode('utf-8')
        sids = IN['sids'][:].view(np.chararray).decode('utf-8')

        header = np.r_[['gene_id'], sids]
        OUT.write('\t'.join(header) + '\n') 
        for c in range(0, genes.shape[0], cs):
            curr_count = IN['counts'][c:min(c+cs, genes.shape[0]), :].astype('str')
            _out_str = []
            for cc in range(c, min(c+cs, genes.shape[0])):
                _out_str.append('\t'.join(np.r_[[genes[cc]], curr_count[cc-c, :]]))
            OUT.write('\n'.join(_out_str) + '\n')
    IN.close()
            
    return 0

if __name__ == "__main__":
    sys.exit(main())
