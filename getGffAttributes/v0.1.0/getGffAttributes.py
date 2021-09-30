#!/usr/bin/env python3

import argparse
import sys
import urllib.parse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description='Get attributes from GFF')

parser.add_argument('--gff', '-gff', type= str, help= 'Input GFF file [%(default)s]', default= '-', nargs= '?')
default = ['ID', 'Name', 'description']
parser.add_argument('--attributes', '-a', nargs= '+', help= 'Attribute keys [%s]' % ', '.join(default), default= default)
parser.add_argument('--type', '-t', help= 'Only parse records of type (3rd column) [%(default)s]', default= 'gene')
parser.add_argument('--rename', '-r', nargs= '+', help= 'Rename attributes in output header, use syntax `attribute:new_name`')
parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

args = parser.parse_args()

def gff_attribute(attribute_str, key, unquote= True):
    # ID=exon_PBANKA_1307600.1-E1;Parent=PBANKA_1307600.1;gene_id=PBANKA_1307600;biotype=exon
    xkey = key + '='
    lst = attribute_str.split(';')
    kv = [x for x in lst if x.startswith(xkey)] 
    assert len(kv) <= 1
    if len(kv) == 0:
        return ''
    kv = kv[0]
    value = kv[len(xkey):]
    if unquote:
        value= urllib.parse.unquote(value)
    return value

if args.gff == '-':
    fin = sys.stdin
else:
    fin = open(args.gff)

rename = {}
for x in args.rename:
    k,v = x.split(':')
    if k not in args.attributes:
        raise Exception('Cannot rename %s because "%s" is not a requested attribute' %(k, k))
    rename[k] = v

header = []
for k in args.attributes:
    if k in rename:
        header.append(rename[k])
    else:
        header.append(k)

print('\t'.join(header))

for line in fin:
    if line.startswith('#'):
        continue
    line = line.strip().split('\t')
    if line[2] != args.type:
        continue
    attrs = line[8]
    values = []
    for k in args.attributes:
         values.append(gff_attribute(attrs, k))
    print('\t'.join(values))
fin.close()        
