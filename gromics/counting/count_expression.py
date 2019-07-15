import sys
import pdb
import pysam
import time
import re
import scipy as sp
import h5py
import pickle
import os

def parse_options(argv):
    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='anno', metavar='FILE', help='annotation file in GTF/GFF3 format', default='-')
    required.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='outfile to store counts in tab delimited format [stdin]', default='-')
    required.add_option('-A', '--alignment', dest='alignment', metavar='FILE', help='alignment in sam or bam format [stdin - sam]', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-H', '--header', dest='header', metavar='STRING', help='header to be used for file in output, default is filename (without extension) [filename]', default='-')
    optional.add_option('-F', '--fields', dest='fields', metavar='STRING', help='annotation fields [exon], comma separated', default='exon')
    optional.add_option('-f', '--filters', dest='filters', metavar='STRING', help='file containing filter maps in hdf5 [-]', default='-')
    optional.add_option('-n', '--filternames', dest='filternames', metavar='STRING', help='list of filter names to use, comma separated, names must be present in t filter hdf5 [names in hdf5 in lex order]', default='-')
    optional.add_option('-t', '--filtertypes', dest='filtertypes', metavar='STRING', help='list of filter types to use, comma separated, either one or same number as filternames, possible types: any, start, all [any]', default='-')
    optional.add_option('-c', '--filtercombs', dest='filtercombs', metavar='STRING', help='list of filter-index combinations: 0,2,4:0,1:... (index relative to filter name list) [one filter in hdf5 at a time]', default='-')
    optional.add_option('-m', '--mask_gene_overlap', dest='mask_gene_overlap', action='store_true', help='mask genomic positions that are annotated with different genes [off]', default=False)
    optional.add_option('-M', '--mask_alternative_overlap', dest='mask_alternative_overlap', action='store_true', help='mask genomic positions that are annotated with both intronic and exonic positions [off]', default=False)
    optional.add_option('-b', '--bam_force', dest='bam_force', action='store_true', help='force BAM as input even if file ending is different from .bam - does not work for STDIN', default=False)
    optional.add_option('-B', '--best_only', dest='best_only', action='store_true', help='count only the best alignment per read [off]', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args(argv)
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def compress_g(g, idx2gene):
    """Find reduced g"""

    g = sorted(g, key = lambda x: len(idx2gene[x]))[::-1]
    g_ = [g[0]]
    seen = idx2gene[g[0]]

    for gg in g[1:]:
        if not all([i in seen for i in idx2gene[gg]]):
            g_.append(gg)
            seen += idx2gene[gg]
    
    return sp.array(g_)

def valid_after_filter(filtermap, filtertype, positions):
    """Description"""

    if filtertype == 'all':
        return not sp.all(filtermap[:, positions])
    elif filtertype == 'start':
        return not filtermap[:, positions[0]]
    elif filtertype == 'any':
        return not sp.any(filtermap[:, positions])
    else:
        return False

def get_filter_settings(options):
    """Parse filter settings from command line options.""" 
    
    if options.filternames != '-':
        filter_names = options.filternames.split(',')
    else:
        hdf_in = h5py.File(options.filters, 'r')
        filter_names = sorted(hdf_in.keys())
        hdf_in.close()

    if options.filtercombs != '-':
        filter_combs = []
        for fc in options.filtercombs.split(':'):
            filter_combs.append(fc.split(','))
            filter_combs[-1] = [int(x) for x in filter_combs[-1]]
    else:
        filter_combs = [[x] for x in range(len(filter_names))]

    if options.filtertypes == '-':
        filter_types = ['any'] * len(filter_names)
    else:
        ft = options.filtertypes.split(',')
        if len(ft) == 1:
            filter_types = [ft[0]] * len(filter_names)
        else:
            assert(len(ft) == len(filter_names))
            filter_types = ft
    
    return (filter_names, filter_combs, filter_types)
    
def compress_counts(count_list, genes):
    """Takes a list of gene IDs and compresses them to a list of tuples"""

    a = 0
    g = 0
    compressed_list = []

    print(" [compressing gene list] ", end=' ', file=sys.stderr)

    while g < len(genes):
        while g < len(genes) and (a == len(count_list) or genes[g] < count_list[a]):
            g += 1
        if g < len(genes):
            b = a
            while a < len(count_list) and genes[g] == count_list[a]:
                a += 1
            compressed_list.append([genes[g], a - b])
            g += 1
    return compressed_list
    
def condense_compressed_counts(compressed_counts):

    t0 = time.time()
    for idx in range(len(compressed_counts)):
        compressed_counts[idx] = sorted(compressed_counts[idx], key = lambda x: x[0])
        for i in range(1, len(compressed_counts[idx])):
            if compressed_counts[idx][i-1][0] == compressed_counts[idx][i][0]:
                compressed_counts[idx][i][1] += compressed_counts[idx][i-1][1]
                compressed_counts[idx][i-1][1] = -1
        compressed_counts[idx] = [x for x in compressed_counts[idx] if x[1] >= 0]
        t1 = time.time() - t0
        print("... done. took %.2f secs" % t1, file=sys.stderr)

    return compressed_counts


def main():
    """Main Program Procedure"""

    options = parse_options(sys.argv)

    options.anno_hdf5 = options.anno
    if options.mask_gene_overlap:
        options.anno_hdf5 += '.mask_go'
    if options.mask_alternative_overlap:
        options.anno_hdf5 += '.mask_ao'

    time_total = time.time()

    ### get filters
    filters = []
    if options.filters != '-':
        ### get filter names
        (filter_names, filter_combs, filter_types) = get_filter_settings(options)

        ### subset to filter names that occur in filter_combs
        filter_combs_flat = list(set([j for sublist in filter_combs for j in sublist]))
        filter_names = [filter_names[i] for i in range(len(filter_names)) if i in filter_combs_flat]
        filter_types = [filter_types[i] for i in range(len(filter_types)) if i in filter_combs_flat]
        filter_combs = [[filter_combs_flat.index(x) for x in j] for j in filter_combs]

        hdf_in = h5py.File(options.filters, 'r')
        for fn in filter_names:
            filters.append(dict())
            for c in hdf_in[fn]:
                filters[-1][c] = hdf_in[fn][c][:]
        hdf_in.close()
    else:
        filter_names = []
        filter_combs = []
        filter_types = []

    anno = dict()
    ### iterate over alignment file(s)
    for fname in options.alignment.split(','):
        options.is_bam = False
        ### open file stream
        if fname == '-':
            infile = sys.stdin
        elif (len(fname) > 3 and fname[-3:] == 'bam') or options.bam_force:
            infile = pysam.Samfile(fname, 'rb')
            options.is_bam = True
        else:
            infile = open(fname, 'r')

        if options.verbose:
            if options.alignment == '-':
                print("Reading alignment from stdin\n", file=sys.stderr)
            else:
                print("Reading alignment from %s\n" % options.alignment, file=sys.stderr)

        if len(anno) == 0:

            ### parse annotation into memory or read from hdf5
            if options.verbose:
                t0 = time.time()
                print('Loading annotation from %s ...' % (options.anno_hdf5 + '.pickle'), file=sys.stderr)
            (idx2gene, gene2idx) = pickle.load(open(options.anno_hdf5 + '.pickle', 'rb'))

            hdf_in = h5py.File(options.anno_hdf5 + '.exons.hdf5', 'r')
            for c in hdf_in:
                anno[c] = hdf_in[c][:]
                if options.verbose:
                    t1 = time.time() - t0
                    print("... %s took %i secs" % (c, t1), file=sys.stderr)
                    t0 = time.time()
            hdf_in.close()

        ### count reads
        counter = 1
        t0 = time.time()
        tmp_count = [[] for i in range(1 + len(filter_combs))]
        compressed_counts = [[] for i in range(1 + len(filter_combs))]
        genes = sorted(idx2gene.keys())
        for line in infile:
            if counter % 10000 == 0:
                print('.', end=' ', file=sys.stderr)
                if counter % 100000 == 0:
                    if len(tmp_count[0]) > 5000000:
                        for idx in range(len(tmp_count)):
                            compressed_counts[idx].extend(compress_counts(sorted(tmp_count[idx]), genes))
                        tmp_count = [[] for i in range(1 + len(filter_combs))]
                    t1 = time.time() - t0
                    print('%i (last 100000 took %.2f secs)' % (counter, t1), file=sys.stderr) 
                    t0 = time.time()

            counter += 1

            if options.is_bam:

                if line.is_unmapped:
                    continue

                if options.best_only and line.is_secondary:
                    continue

                chrm = infile.getrname(line.tid)
                pos = line.pos #- 1 is already 0 based in pysam now
                broken = False

                #read_pos = line.positions --> alternative to code below
                read_pos = []
                for o in line.cigar:
                    #if o[0] in [0, 2]: ### new behavior -> ignore deletions for counting
                    if o[0] == 0:
                        read_pos.extend(list(range(pos, pos + o[1])))

                    if not o[0] in [1, 4, 5]:
                        pos += o[1]

                try:
                    g = sp.unique(anno[chrm][read_pos])    
                except IndexError:
                    try:
                        read_pos = read_pos[(read_pos >= 0) & (read_pos < anno[chrm].shape[0])]
                        g = sp.unique(anno[chrm][read_pos])    
                    except:
                        continue

                g = g[g > 1]
                if g.shape[0] == 0:
                    continue

                ### resolve overlapping genes if we haven't masked them
                if not options.mask_gene_overlap and g.shape[0] > 1:
                    g = compress_g(g, idx2gene)
                tmp_count[0].extend(g)
                #print line.qname

                ### get validity for each filter
                if len(filter_names) > 0:
                    is_valid = sp.ones((len(filter_names), ), dtype = 'bool')
                    for idx, fn in enumerate(filter_names):
                        try:
                            is_valid[idx] = valid_after_filter(filters[idx][chrm], filter_types[idx], read_pos)
                        except KeyError:
                            continue
                    ### generate filter combination counts
                    for idx, comb in enumerate(filter_combs):
                        if sp.all(is_valid[comb]):
                            tmp_count[idx + 1].extend(g)
            else:
                sl = line.strip().split('\t')
                if len(sl) < 9:
                    print("ERROR: invalid SAM line\n%s" % line, file=sys.stderr)
                    sys.exit(1)

                (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
                size = [int(i) for i in size]

                #chrm = sl[2].replace('chr', '')
                chrm = sl[2]
                pos = int(sl[3]) - 1
                broken = False

                ## is unmapped ?
                if (int(sl[1]) & 4) == 4:
                    continue

                ## is secondary ?
                if options.best_only and (int(sl[1]) & 256 == 256):
                    continue

                for o in range(len(op)):
                    if op[o] in ['M', 'D']:
                        for p in range(size[o]):
                            try:
                                g = anno[chrm][pos + p]
                                if g > 1:
                                    tmp_count.append(g)
                                    break
                            except KeyError:
                                continue
                            except IndexError:
                                if chrm in ['chrM', 'M', 'chrM_rCRS']:
                                    continue
                                else:
                                    print('ERROR: %i exceeds length of %s' % (pos + p, chrm), file=sys.stderr)
                    if broken:
                        break
                    if not op[o] in ['H', 'I']:
                        pos += size[o]
                           
        ### close file stream
        if not fname == '-':
            infile.close()

        ### compress remaining counts
        for idx in range(len(tmp_count)):
            compressed_counts[idx].extend(compress_counts(sorted(tmp_count[idx]), genes))
        del tmp_count

        ### condense count lists
        print("Sorting and condensing compressed list ...", file=sys.stderr)
        compressed_counts = condense_compressed_counts(compressed_counts)

        ### resolve gene combinations
        for idx in range(len(compressed_counts)):
            extend_list = []
            for a in range(len(compressed_counts[idx]) -1, -1, -1):
                if len(idx2gene[compressed_counts[idx][a][0]]) > 1:
                    for g in idx2gene[compressed_counts[idx][a][0]]:
                        extend_list.append([gene2idx[tuple([g])], compressed_counts[idx][a][1]])
                    del compressed_counts[idx][a]
            compressed_counts[idx].extend(extend_list)
        compressed_counts = condense_compressed_counts(compressed_counts)

        ### remove gene IDs that encode combinations
        genes = [genes[i] for i in range(len(genes)) if len(idx2gene[genes[i]]) < 2]

        ### report counts to outfile
        if options.verbose:
            print("Summarizing gene counts ...", file=sys.stderr)

        for idx in range(len(compressed_counts)):
            if idx > 0 and (idx - 1) < len(filter_combs):
                comb_tag = '_'.join([filter_names[i] for i in filter_combs[idx - 1]])
            else:
                comb_tag = ''

            outfile = open(options.outfile + comb_tag, 'w')
            a = 0
            g = 0
            ### writing header
            if options.header == '-':
                print('\t'.join(['gene_id', fname.rsplit('.', 1)[0]]), file=outfile)
            else:
                print('\t'.join(['gene_id', options.header]), file=outfile)

            ### seek to first position that mapped to gene (0 means not gene found)
            while g < len(genes):
                while g < len(genes) and (a == len(compressed_counts[idx]) or genes[g] < compressed_counts[idx][a][0]):
                    print('%s\t0' % idx2gene[genes[g]][0], file=outfile)
                    if options.verbose and g % 100 == 0:
                        print("%.2f / 100 percent \r" % (float(g) / len(genes) * 100), end=' ', file=sys.stderr)
                    g += 1
                while a < len(compressed_counts[idx]) and g < len(genes) and genes[g] == compressed_counts[idx][a][0]:
                    print('%s\t%i' % (idx2gene[genes[g]][0], compressed_counts[idx][a][1]), file=outfile)
                    a += 1
                    g += 1
                    if options.verbose and g % 100 == 0:
                        print("%.2f / 100 percent \r" % (float(g) / len(genes) * 100), end=' ', file=sys.stderr)
            outfile.close()

        if options.verbose:
            t1 = time.time() - time_total
            print("\n... done - total run took %i secs." % t1, file=sys.stderr)

    return 0

if __name__ == "__main__":
    sys.exit(main())
