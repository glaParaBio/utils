#!/usr/bin/env python3

import os
import argparse
import subprocess
import tempfile
import sys
import re
import math
import collections

parser = argparse.ArgumentParser(description= """
Calculate the log2 fold enrichment of target intervals using the read depth 
profile from Genrich.

If the target intervals are is the narrowPeak file from Genrich, the output has header:

    #chrom       
    start  
    end    
    peak_id     
    score   
    strand 
    auc          
    pvalue       
    fdr          
    summit 
    logfc:avg            # LogFC averaged across samples
    logfc:<input.bam.1>  # logFC of input sample 1, 2,... N
    logfc:<input.bam.2>  #
    logfc:<input.bam.N>  #

* MEMO: Execute Genrich -k option to produce the bedgraph profile of read depth

* Requires bedtools

* Tested with Genrich v0.6.0
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--bedgraph', '-b',
                    required= True,
                    metavar= 'FILE',
                    help='''\
Input bedgraph from Genrich''')

parser.add_argument('--narrowPeak', '-n',
                    default= '-',
                    metavar= 'FILE',
                    help='''\
Input peak file in bed format. Typically this is the narrowPeak
file from Genrich [%(default)s]''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.3.0')

# bedgraph file has columns:
# chr     start   end     experimental    control -log(p)
# narroPeak file has columns:
# chrom  start   end     peak_id score   strand  auc     pvalue  fdr     summit

class Bdgfile:
    def __init__(self, header_line):
        assert line.startswith('# experimental file: ')
        self.experiment = re.sub(';.*', '', re.sub('^# experimental file: ', '', line))
        assert len(self.experiment) > 0
        self.bdgfile = tempfile.NamedTemporaryFile(prefix= 'genrich_logfc.', suffix= '.tmp.bedgraph', delete= False, dir= os.getcwd()).name
        self.npos = 0.0
        self.mass = 0.0
        self.mass_ctrl = 0.0
        self.fout= open(self.bdgfile, 'w')
        sys.stderr.write('Processing experiment %s\n' % self.experiment) 

    def logfc(self, narrowPeak, delete_bdg= True):
        self.fout.close()

        cmd = r"""
        set -euf -o pipefail
        intersectBed -a {bdgfile} -b {narrowPeak} -sorted -wa -wb \
        | groupBy -g 6,7,8 -c 4,5 -o max,mean
        """.format(bdgfile= self.bdgfile, narrowPeak= narrowPeak)
        
        p = subprocess.Popen(cmd, shell= True, executable= 'bash', stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout,stderr= p.communicate()

        if delete_bdg:
            os.remove(self.bdgfile)
        if not p.returncode == 0:
            raise Exception(stderr.decode())

        avg_depth = self.mass / self.npos
        avg_depth_ctrl = self.mass_ctrl / self.npos
        peak_logfc = []
        for line in stdout.decode().strip().split('\n'):
            line = line.strip().split('\t')
            assert len(line) == 5
            depth = float(line[3])
            depth_ctrl = float(line[4])
            if depth_ctrl == 0:
                sys.stderr.write('Warning: Zero depth in control for %s\n' % '\t'.join(line))
                depth_ctrl = avg_depth_ctrl / 10 # This is quiet arbitrary
            if depth == 0:
                sys.stderr.write('Warning: Zero depth in treatment for %s\n' % '\t'.join(line))
                depth = avg_depth / 10 # This is quiet arbitrary

            logfc = math.log2(depth / avg_depth) - math.log2(depth_ctrl / avg_depth_ctrl)
            peak_logfc.append(logfc)
        return peak_logfc

def extend_narrowPeak(narrowPeak, peak_logfc):
    """Read narrowPeak file and add columns for logfc
    """
    np = open(narrowPeak).readlines()
    i = 0
    while True:
        if np[i].startswith('#'):
            sys.stdout.write(np[i])
            del(np[i])
            i += 1
        else:
            break

    ncols = len(np[0].strip().split('\t'))
    if ncols == 10:
        header = ['#chrom', 'start', 'end', 'peak_id', 'score', 'strand', 'auc', 'pvalue', 'fdr', 'summit']
    else:
        header = ['x' + str(i) for i in range(ncols)]
        header[0] = '#' + header[0]
    header.append('logfc:avg')

    for peakfile in peak_logfc:
        header.append('logfc:' + peakfile)
        assert len(np) == len(peak_logfc[peakfile])
    print('\t'.join(header))

    i = 0
    for i in range(len(np)):
        logfcs = []
        for peakfile in peak_logfc:
            logfcs.append(peak_logfc[peakfile][i])
        line = np[i].strip().split('\t')
        fold_changes = [2**x for x in logfcs]
        line.append('%.4f' % math.log2(sum(fold_changes) / len(fold_changes)))
        [line.append('%.4f' % x) for x in logfcs]
        print('\t'.join(line))

if __name__ == "__main__":
    args = parser.parse_args()
    # The bedgraph from Genrich is a concatenation of all the input experimental files. 
    # 
    #   # experimental file: genrich/WT-507-10h-3.ATAC.nsort.bam; control file: genrich/gDNA.bam
    #   chr     start   end     experimental    control -log(p)
    #
    narrowPeak = tempfile.NamedTemporaryFile(prefix= 'genrich_logfc.', suffix= '.tmp.narrowPeak', delete= False, dir= os.getcwd()).name
    fout = open(narrowPeak, 'w')
    if args.narrowPeak == '-':
        fin = sys.stdin
    else:
        fin = open(args.narrowPeak)
    for line in fin:
        fout.write(line)
    fout.close()
    fin.close()

    np = open(narrowPeak).readlines()
    np = [x for x in np if not x.strip().startswith('#')]
    if len(np) == 0:
        sys.stderr.write('Exiting: No intervals found in %s\n' % narrowPeak)
        sys.exit(0)

    with open(args.bedgraph) as bdg:
        first = True
        peak_logfc = collections.OrderedDict()
        for line in bdg:
            line = line.strip()
            if line.startswith('#'):
                if first:
                    bdg = Bdgfile(line)
                    first = False
                else:
                    peak_logfc[bdg.experiment] = bdg.logfc(narrowPeak)
                    bdg = Bdgfile(line)
                continue

            line= line.split('\t')
            if line[0] == 'chr':
                exp_idx = line.index('experimental')
                ctrl_idx = line.index('control')
            else:
                for i in range(int(line[1]), int(line[2])):
                    bdg.npos += 1
                    bdg.mass += float(line[exp_idx])
                    bdg.mass_ctrl += float(line[ctrl_idx])
                    out= [line[0], str(i), str(i+1), line[exp_idx], line[ctrl_idx]]
                    bdg.fout.write('\t'.join(out) + '\n')
        peak_logfc[bdg.experiment] = bdg.logfc(narrowPeak)

    extend_narrowPeak(narrowPeak, peak_logfc)
    os.remove(narrowPeak)
