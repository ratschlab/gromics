import sys
import scipy as sp
import h5py
import os
import re

def parse_options(argv):

    """Parses options from the command line """

    from argparse import ArgumentParser

    parser = ArgumentParser(prog='compute_lib_size')
    parser.add_argument('-i', '--input', dest='infile', metavar='STR', help='expression counts in hdf5 format', default='-', required=True)
    parser.add_argument('-a', '--annotation', dest='annotation', metavar='STR', help='annotation file (needed for coding / autosome filtering) []', default='-')
    parser.add_argument('--coding', dest='coding', action='store_true', help='only use coding genes for normalization', default=False)
    parser.add_argument('--autosomes', dest='autosomes', metavar='LIST', nargs='+', help='list of autosomes to exclude from normalization []', default=[])
    
    return parser.parse_args(argv[1:])


def get_tags_gff3(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags


def get_tags_gtf(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags


def main():

    options = parse_options(sys.argv)

    ### get the data
    print('loading expression data from ' + options.infile)
    IN = h5py.File(options.infile, 'r')
    expression = IN['counts'][:]
    genes = IN['gids'][:].view(sp.chararray).decode('utf-8')
    strains = IN['sids'][:].view(sp.chararray).decode('utf-8')
    IN.close()

    ### get list of protein coding genes
    coding = []
    if options.coding or len(options.autosomes) > 0:
        for line in open(options.annotation, 'r'):
            if line[0] in ['#']:
                continue
            sl = line.strip().split('\t')
            if sl[2].lower() != 'gene':
                continue
            if options.annotation.lower().endswith('.gtf'):
                tags = get_tags_gtf(sl[8])
                gene = tags['gene_id']
                is_coding = tags['gene_type'] == 'protein_coding'
            elif options.annotation.lower().endswith('.gff'):
                tags = get_tags_gff3(sl[8])
                try:
                    gene = tags['ID'] 
                    is_coding = tags['gene_type'] == 'protein_coding'
                except KeyError:
                    continue
            if not options.coding or is_coding:
                coding.append([sl[0], gene])
        coding = sp.array(coding)        

        ### filter for autosomes
        if len(options.autosomes) > 0:
            k_idx = sp.where(~sp.in1d(coding[:, 0], options.autosomes))[0]
            coding = coding[k_idx, :]
        coding = coding[:, 1]

        ### filter expression
        k_idx = sp.where(sp.in1d(genes, coding))[0]
        genes = genes[k_idx]
        expression = expression[k_idx, :]

    ### compute normalizers
    libsize_uq = sp.array([sp.percentile(x, 75) for x in expression.T])
    libsize_tc = expression.sum(axis=0)

    s_idx = sp.argsort(libsize_uq)[::-1]

    out = open(re.sub('.hdf5$', '', options.infile) + '.libsize.tsv', 'w')
    print('sample\tlibsize_75percent\tlibsize_total_count', file=out)
    for i in s_idx:
        print('\t'.join([strains[i], str(libsize_uq[i]), str(libsize_tc[i])]), file=out)
    out.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
