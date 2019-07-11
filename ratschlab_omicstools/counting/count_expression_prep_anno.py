import scipy as sp
import sys
import time
import h5py
import pickle

def parse_options(argv):
    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='anno', metavar='FILE', help='annotation file in GTF/GFF3 format', default='-')
    required.add_option('-g', '--genome', dest='genome', metavar='FILE', help='reference genome in FASTA format', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-m', '--mask_gene_overlap', dest='mask_gene_overlap', action='store_true', help='mask genomic positions that are annotated with different genes [off]', default=False)
    optional.add_option('-M', '--mask_alternative_overlap', dest='mask_alternative_overlap', action='store_true', help='mask genomic positions that are annotated with both intronic and exonic positions [off]', default=False)
    optional.add_option('-F', '--fields', dest='fields', metavar='STRING', help='annotation fields [exon], comma separated', default='exon')
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args(argv)
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def get_tags_gff(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.split(';'):
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

def parse_anno(options, format='gff'):
    """This function reads the gff3 input file and returns the information in an
       internal data structure"""

    anno = dict()
    idx2gene = dict()
    gene2idx = dict()

    if options.verbose:
        print("Parsing annotation from %s ..." % options.anno, file=sys.stderr)
    
    ### initial run to get the transcript to gene mapping
    if format == 'gff':
        if options.verbose:
            print("... init structure", file=sys.stderr)

        trans2gene = dict() ### dict with: keys = transcript IDs, values = gene IDs
        for line in open(options.anno, 'r'):
            if line[0] == '#':
                continue
            sl = line.strip().split('\t')
            if sl[2] in ['mRNA', 'transcript', 'mrna', 'miRNA', 'tRNA', 'snRNA', 'snoRNA', 'ncRNA', 'mRNA_TE_gene', 'rRNA', 'pseudogenic_transcript', 'transposon_fragment']:
                tags = get_tags_gff(sl[8])
                trans2gene[tags['ID']] = tags['Parent']

    ### get contig lengths
    seq = 0
    curr_id = ''
    contigs = dict()
    for line in open(options.genome, 'r'):
        if line[0] == '>':
            if seq > 0:
                contigs[curr_id] = seq
                seq = 0
            curr_id = line.strip().split()[0][1:]
            continue
        seq += len(line.strip())
    if seq > 0:
        contigs[curr_id] = seq

    ### init genome structure
    for c in contigs:
        if options.verbose:
            print('reserving memory for contig %s of len %s' % (c, contigs[c]), file=sys.stderr)
        anno[c] = sp.zeros((contigs[c] + 1,), dtype = 'int32')

    ### init list of  considered GFF fields
    fields = options.fields.split(',')

    ### generate a list of exons with attached gene/transcript information
    ### one list per chromsome
    counter = 1
    gene_counter = 2 ### 0 is default for no coverage and 1 is mask for overlap

    exons = dict() # contains the exon list per transcript, only need this for mask_alternative_overlap

    t0 = time.time()
    for line in open(options.anno, 'r'):
        if options.verbose and counter % 10000 == 0:
            print('.', end=' ', file=sys.stderr)
            if counter % 100000 == 0:
                t1 = time.time() - t0
                print("%i - took %.2f secs" % (counter, t1), file=sys.stderr)
                t0 = time.time()
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        if not sl[2] in fields:
            continue

        if sl[2] != 'exon':
            print('Currently only >exon< is supported', file=sys.stderr)
            sys.exit(1)

        if format == 'gff':
            tags = get_tags_gff(sl[8])
            trans_id = tags['Parent']
            gene_id = trans2gene[trans_id]
        elif format == 'gtf':
            tags = get_tags_gtf(sl[8])
            trans_id = tags['transcript_id']
            gene_id = tags['gene_id']
        else:
            sys.stderr.write('ERROR: Unknown format: %s \n' % format)
            sys.exit(1)

        if tuple([gene_id]) not in gene2idx:
            gene2idx[tuple([gene_id])] = gene_counter
            idx2gene[gene_counter] = tuple([gene_id])
            gene_counter += 1

        ### store for each position of the transcriptome a tuple containing all overlapping gene IDs
        ### assume positions are 1 based and in closed intervals
        try:
            start = int(sl[3]) - 1
        except ValueError:
            start = 0
        try:
            stop = int(sl[4])
        except ValueError:
            stop = 1

        chrm = sl[0]
        if chrm == 'chrM_rCRS':
            chrm = 'chrM'

        if not chrm in exons:
            exons[chrm] = dict()

        if options.mask_alternative_overlap:
            try:
                exons[chrm][trans_id].append([start, stop])
            except KeyError:
                exons[chrm][trans_id] = [[start, stop]]

        ### check, if there is already a different gene ID present, form a combination ID
        if sp.any(anno[chrm][start:stop] > 0):
            for p in range(start, stop):
                if anno[chrm][p] == 0:
                    new_set = tuple([gene_id])
                else:
                    new_set = tuple(set(idx2gene[anno[chrm][p]]) | set([gene_id]))
                try:
                    anno[chrm][p] = gene2idx[new_set]
                except KeyError:
                    anno[chrm][p] = gene_counter
                    gene2idx[new_set] = gene_counter
                    idx2gene[gene_counter] = new_set
                    gene_counter += 1
        else:
            anno[chrm][start:stop] = sp.array([gene2idx[tuple([gene_id])]] * (stop - start), dtype = 'int32')
    if options.verbose:
        print("... done", file=sys.stderr)

    ### mask all positions in the genome, where we have more than one annotated gene
    if options.mask_gene_overlap:
        total_pos = 0
        total_masked = 0
        if options.verbose:
            print('\nMasking positions due to gene overlap:', file=sys.stderr)
        for c in anno:
            masked_pos = 0
            p_idx = sp.where(anno[c] > 1)[0]
            pos = p_idx.shape[0]
            #print >> sys.stderr, 'found %i positions' % p_idx.shape[0]
            for p in p_idx:
                if len(idx2gene[anno[c][p]]) > 1:
                    anno[c][p] = 1
                    masked_pos += 1
            total_pos += pos
            total_masked += masked_pos
            if options.verbose:
                print('\t%s: %i (%i) masked (total) - %.2f %%' % (c, masked_pos, pos, masked_pos / float(max(1, pos)) * 100), file=sys.stderr)
        if options.verbose:
            print("Total positions: %i\nMasked positions: %i (%.2f %%)" % (total_pos, total_masked, total_masked / float(max(1, total_pos)) * 100), file=sys.stderr)
            print("... done", file=sys.stderr)

    ### mask all positions in the genome, where exonic and intronic positions are annotated
    if options.mask_alternative_overlap:
        if options.verbose:
            print('\nMasking positions due to exon/intron overlap:', file=sys.stderr)
        for c in exons:
            masked_pos = 0
            for t in exons[c]:
                if len(exons[c][t]) < 2:
                    continue
                ### pre-process exon (sort by start)
                tmp = sp.array(exons[c][t], dtype='int')
                s_idx = sp.argsort(tmp[:, 0])
                tmp = tmp[s_idx, :]
                ### mask positions that are intronic and exonic
                for e in range(1, tmp.shape[0]):
                    p_idx = sp.where(anno[c][tmp[e - 1, 1]:tmp[e, 0]] > 1)[0]
                        
                    if p_idx.shape[0] > 0:
                        anno[c][p_idx + tmp[e - 1, 1]] = 1
                        masked_pos += p_idx.shape[0]
            total_masked += masked_pos
            if options.verbose:
                print('\t%s: %i pos masked' % (c, masked_pos), file=sys.stderr)
        if options.verbose:
            print('Masked positions: %i' % total_masked, file=sys.stderr)
            print("... done", file=sys.stderr)

    
    if options.verbose:
        print("Storing exon array in HDF5 %s ..." % (options.anno_hdf5 + '.exons.hdf5'), file=sys.stderr)

    ### store annotation in hdf5
    hdf_out = h5py.File(options.anno_hdf5 + '.exons.hdf5', 'w')
    for c in list(anno.keys()):
        hdf_out.create_dataset(name = c, data = anno[c])
    hdf_out.close()

    if options.verbose:
        print("... pickling gene ID map", file=sys.stderr)

    pickle.dump((idx2gene, gene2idx), open(options.anno_hdf5 + '.pickle', 'wb'))

    if options.verbose:
        print("... done", file=sys.stderr)

    return (anno, idx2gene, gene2idx)

def main():

    options = parse_options(sys.argv)
    options.verbose = True

    options.anno_hdf5 = options.anno
    if options.mask_gene_overlap:
        options.anno_hdf5 += '.mask_go'
    if options.mask_alternative_overlap:
        options.anno_hdf5 += '.mask_ao'

    if options.anno.lower().endswith('gff3') or options.anno.lower().endswith('gff'):
        ### read annotation from GFF3
        parse_anno(options, format='gff')
    elif options.anno.lower().endswith('gtf'):
        ### read annotation from GTF
        parse_anno(options, format='gtf')
    else:
        sys.stderr.write('ERROR: Unknown format for input annotation file %s\n' % options.anno) 
        sys.exit(1)

    return 0


if __name__ == "__main__":
    sys.exit(main())
