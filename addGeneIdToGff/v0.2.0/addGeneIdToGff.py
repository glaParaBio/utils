#!/usr/bin/env python3

import gffutils
import argparse
import sys

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description='Add gene ID attribute to GFF features (i.e. to exon, CDS etc)')

parser.add_argument('gff', type= str, help= 'Input GFF file [%(default)s]', default= '-', nargs= '?')

parser.add_argument('--gene-feature', '-g', type= str, help= 'List of feature types (column 3) identifying a gene %(default)s',
        default= ['gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'], nargs= '+')
parser.add_argument('--geneid', '-id', type= str, help= 'Name of gene ID attribute [%(default)s]', default= 'gene_id')
parser.add_argument('--overwrite-id', '-o', help= 'Overwrite the geneid attribute if it exists. Default is to exit with error', action= 'store_true')
parser.add_argument('--add-tss', '-tss', help= 'Add a record to mark the TSS of features of type [%(default)s]. Use "" to skip', default= 'mRNA')
parser.add_argument('--version', action='version', version='%(prog)s 0.2.0')

args = parser.parse_args()

if args.gff == '-':
    gff = sys.stdin.readlines()
else:
    gff = open(args.gff).readlines()

db = gffutils.create_db(''.join(gff), ':memory:', from_string= True, merge_strategy= 'create_unique')

for d in db.directives:
    print('##%s' % d)
 
tss_id = 1
for feature in db.all_features():
    parents = list(db.parents(feature))
    #if not args.overwrite_id and args.geneid in feature.attributes:
    #    print(feature.attributes)
    #    sys.stderr.write('Error: attribute key "%s" already exists\n' % args.geneid)
    #    sys.exit(1)
    parent_gene = []
    for p in parents:
        if p.featuretype in args.gene_feature:
            parent_gene.append(p)
    assert len(parent_gene) < 2
    if feature.featuretype in args.gene_feature:
        feature.attributes[args.geneid] = feature.attributes['ID']
    elif len(parent_gene) == 1:
        parent_gene = parent_gene[0]
        
        gene_id = parent_gene.attributes['ID']
        assert len(gene_id) == 1

        if not args.overwrite_id and args.geneid in feature.attributes and gene_id != feature.attributes[args.geneid]:
            sys.stderr.write('Error: attribute key "%s" already exists and it is different from the new one.\nUse --overwrite-id to overwrite\n' % args.geneid)
            sys.exit(1)
        feature.attributes[args.geneid] = gene_id
    if args.add_tss != '' and feature.featuretype == args.add_tss:
        if feature.strand == '+':
            start = feature.start
        elif feature.strand == '-':
            start = feature.end
        else:
            raise Exception('Invalid strand: %s' % feature.strand)
        if args.geneid not in feature.attributes:
            sys.stderr.write('Cannot assign TSS to feature:\n%s\nbecause it does not have attibute "%s" assigned.\nPerhaps you need to edit argument --geneid to ensure the parent of this feature is included\n' % (feature, args.geneid))
            sys.exit(1)
        tss = gffutils.Feature(seqid= feature.seqid, source= feature.source, featuretype= 'TSS', start= start, end= start, strand= feature.strand, attributes= {'ID': ['TSS_%s' % tss_id], 'Parent': feature['ID'], args.geneid: feature.attributes[args.geneid]})
        tss_id = tss_id + 1
        print(tss)
    print(feature)
