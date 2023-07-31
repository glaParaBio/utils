<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Example](#example)
* [Setup](#setup)

<!-- vim-markdown-toc -->

# Description

Help from `./pafplot.R -h` (printout here may be outdated):

```
usage: ./pafplot.R [-h] [--input INPUT] --output OUTPUT [--ncol NCOL]
                   [--nrow NROW] [--min-qcovpct MIN_QCOVPCT] [--mapq MAPQ]
                   [--pagesize PAGESIZE [PAGESIZE ...]]
                   [--panel-spacing PANEL_SPACING] [-v]

DESCRIPTION 
Plot one assembly vs another contig by contig using the PAF format as input

USAGE 

* Prepare alignment (note that `paftools sam2paf` converts BAM to PAF) 

    minimap2 -x asm5 ref.fa assembly.fa > aln.paf 

* Plot 4 contigs per A4 page possibly spanning multiple pages 

    pafplot.R -i aln.paf -o out.pdf --nrow 4 

* Plot selected contigs and reduce pagesize to 21x16 cm

    grep -P "contig_34|scaffold_68" aln.paf | pafplot.R -o out.pdf --nrow 2 -s 21 16 

Version 0.1.0

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input in PAF format. Use - for reading from stdin [-]
  --output OUTPUT, -o OUTPUT
                        File for PDF output
  --ncol NCOL, -c NCOL  Number of columns per page [1]
  --nrow NROW, -r NROW  Number of rows per page. "auto" will fit all in one page [auto]
  --min-qcovpct MIN_QCOVPCT, -p MIN_QCOVPCT
                        Exclude alignments that cover less than this percent
                        of the contig [2%]
  --mapq MAPQ, -q MAPQ  Do not draw connecting lines when the mapping quality is below mapq [5]
  --pagesize PAGESIZE [PAGESIZE ...], -s PAGESIZE [PAGESIZE ...]
                        One or two arguments for WIDTH and HEIGHT in cm for the page.
                        If only one number is given it used for both width and height.
                        "A4" is shortcut for "21.0 29.7" and "A4landscape" for A4 in
                        landscape format - can be a prefix and case insensitive
                        e.g. "A4L" [A4]
  --panel-spacing PANEL_SPACING, -ps PANEL_SPACING
                        Vertical panel spacing in units of lines of text
  -v, --version         show program's version number and exit
```

# Example

See [example.pdf](test_data/example.pdf), which you should be able to reproduce with:

```
zcat test_data/aln.paf.gz | ./pafplot.R -i - -o example.pdf
```

# Setup

[requirements.txt](requirements.txt) lists the dependencies; package versions
other than the specified ones are likely to work fine. This is an R script, install
R and the required packages using your favourite method. With conda/mamba:

```
mamba install --file requirements.txt 
```
