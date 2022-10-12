#!/usr/bin/env python3

from signal import signal, SIGPIPE, SIG_DFL
import argparse
import os
import pathlib
import pandas

signal(SIGPIPE, SIG_DFL) 

parser = argparse.ArgumentParser(description= "Concatenate the orthologues files from orthofinder and make them in 'long' format")
parser.add_argument("orthofinder_basedir", help= "Result directory produced by orthofinder; it should contain the 'Orthologues' subdirectory with files to concatenate")
parser.add_argument('--include-unassigned', '-u', action= 'store_true', help= "Include genes without orthologues")
parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

args = parser.parse_args()

tsv = sorted(pathlib.Path(os.path.join(args.orthofinder_basedir, 'Orthologues')).glob('**/*__v__*.tsv'))

print('\t'.join(['orthogroup', 'species1', 'pid1', 'species2', 'pid2']))

for x in tsv:
    assert x.name.count('__v__') == 1
    with open(x) as fin:
        for line in fin:
            line = line.strip().split('\t')
            if line[0] == 'Orthogroup':
                sp1, sp2 = line[1:]
                continue
            og = line[0]
            pid1 = line[1].split(', ') 
            pid2 = line[2].split(', ') 
            for p1 in pid1:
                for p2 in pid2:
                    print('\t'.join([og, sp1, p1, sp2, p2]))

if args.include_unassigned:
    tsv = os.path.join(args.orthofinder_basedir, 'Orthogroups/Orthogroups_UnassignedGenes.tsv')
    ung = pandas.read_csv(tsv, sep= '\t')
    ung = pandas.melt(ung, id_vars= ['Orthogroup'], var_name= 'species1', value_name= 'pid1')
    ung.rename(columns= {'Orthogroup': 'orthogroup'}, inplace=True)
    ung = ung[ung.pid1.notna()]
    for idx,row in ung.iterrows():
        print('\t'.join([row.orthogroup, row.species1, row.pid1, 'NA', 'NA']))
        print('\t'.join([row.orthogroup, 'NA', 'NA', row.species1, row.pid1]))
