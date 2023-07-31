#!/usr/bin/env Rscript

suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))

VERSION <- '0.1.0'

docstring<- sprintf('DESCRIPTION \\n\\
Plot one assembly vs another contig by contig using the PAF format as input\\n\\
\\n\\
USAGE \\n\\
\\n\\
* Prepare alignment (note that `paftools sam2paf` converts BAM to PAF) \\n\\
\\n\\
    minimap2 -x asm5 ref.fa assembly.fa > aln.paf \\n\\
\\n\\
* Plot 4 contigs per A4 page possibly spanning multiple pages \\n\\
\\n\\
    pafplot.R -i aln.paf -o out.pdf --nrow 4 \\n\\
\\n\\
* Plot selected contigs and reduce pagesize to 21x16 cm\\n\\
\\n\\
    grep -P "contig_34|scaffold_68" aln.paf | pafplot.R -o out.pdf --nrow 2 -s 21 16 \\n\\
\\n\\
Version %s', VERSION)

parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')

def <- '-'
parser$add_argument('--input', '-i', help= sprintf('Input in PAF format. Use - for reading from stdin [%s]', def), default= def)
parser$add_argument('--output', '-o', help= 'File for PDF output', required= TRUE)

def <- 1
parser$add_argument('--ncol', '-c', help= sprintf('Number of columns per page [%s]', def), type= 'integer', default= def)

def <- 'auto'
parser$add_argument('--nrow', '-r', help= sprintf('Number of rows per page. "auto" will fit all in one page [%s]', def), type= 'character', default= def)

def <- 2.0
parser$add_argument('--min-qcovpct', '-p', help= sprintf('\\
Exclude alignments that cover less than this percent\\n\\
of the contig [%s%%%%]', def), type= 'double', default= def)

def <- 5
parser$add_argument('--mapq', '-q', help= sprintf('\\
Do not draw connecting lines when the mapping quality is below mapq [%s]', def), type= 'double', default= def)

def <- 'A4'
parser$add_argument('--pagesize', '-s', help= sprintf('\\
One or two arguments for WIDTH and HEIGHT in cm for the page.\\n\\
If only one number is given it used for both width and height.\\n\\
"A4" is shortcut for "21.0 29.7" and "A4landscape" for A4 in\\n\\
landscape format - can be a prefix and case insensitive\\n\\
e.g. "A4L" [%s]', def), 
type= 'character', default= def, nargs= '+')

parser$add_argument('--panel-spacing', '-ps', help= sprintf('\\
Vertical panel spacing in units of lines of text'),
type= 'double')

# NB: argparse v1.1.1+ required for -v option to work.
parser$add_argument("-v", "--version", action= 'version', version= VERSION)

psize <-  function(size){
    # Return a vector of length 2 with page size in *inches* for width and
    # height
    size <- toupper(size)
    if(length(size) > 2){
        stop('Invalid argument for page size. Got length > 2')
    }
    if(length(size) == 1 & size[1] == 'A4'){
        return(c(21.0, 29.7)/2.54)
    }
    if(length(size) == 1 & grepl(sprintf('^%s', size[1]), 'A4LANDSCAPE')){
        # A4L also 
        return(c(29.7, 21.0)/2.54)
    }
    size <- as.numeric(size)
    if(NA %in% size){
        stop('Invalid argument for page size')
    }
    if(length(size) == 1){
        return(c(size, size)/2.54)
    } else {
        return(size/2.54)
    }
}

ylabeller <- function(labels) {
    outlabs <- rep(NA, length(labels))
    for(i in 1:length(labels)){
        lab <- labels[i]
        if(startsWith(lab, prefix)){
            outlabs[i] <- substring(lab, nchar(prefix) + 2)
        } else {
            outlabs[i] <- lab
        }
    }
    return(outlabs)
}

if(sys.nframe() == 0){
    # Script is being executed from the command line

    # ---------------- [Validate arguments] -------------
    options(scipen= 10)

    xargs <- parser$parse_args()

    pagesize <- psize(xargs$pagesize)

    col.names= c('qname', 'qlen', 'qstart', 'qend', 'strand', 'sname', 'slen', 'sstart', 'send', 'mapq')
    xcut <- 'cut -f 1-9,12'
    if(xargs$input == '-'){
        paf <- fread(cmd= sprintf('cat /dev/stdin | %s', xcut), col.names= col.names)
    } else {
        paf <- fread(cmd= sprintf('%s %s', xcut, xargs$input), col.names= col.names)
    }
    naln <- paf[, list(.N), by= qname]
    cat(sprintf('Found %s alignments in %s contigs\n', sum(naln$N), length(unique(paf$qname))))
        
    #paf <- paf[mapq >= xargs$mapq]
    #if(nrow(paf) == 0){
    #    cat('No alignment found\n')
    #    quit()
    #}
    paf[, qcovpct := 100 * (qend-qstart)/qlen]
    
    pafx <- paf[qcovpct > xargs$min_qcovpct]
    cat(sprintf('%s alignments excluded\n', sum(naln$N) - nrow(pafx)))
    if(nrow(pafx) == 0){
        cat('No alignment left after filtering for query coverage\n')
        quit()
    }
    pafx <- merge(pafx, naln, by= 'qname', sort= FALSE)
    pafx <- merge(pafx, pafx[, list(shown= .N), by= 'qname'], by= 'qname', sort= FALSE)

    bp_unit <- 'bp'
    if(max(pafx$slen) > 10000000 | max(pafx$qlen) > 10000000){
        paf[, c('qstart', 'qend', 'sstart', 'send', 'slen', 'qlen') :=
            list(qstart/1000000, qend/1000000, sstart/1000000, send/1000000, slen/1000000, qlen/1000000)]
        bp_unit <- 'Mb'
    }
    else if(max(pafx$slen) > 10000 | max(pafx$qlen) > 10000){
        pafx[, c('qstart', 'qend', 'sstart', 'send', 'slen', 'qlen') := list(qstart/1000, qend/1000, sstart/1000, send/1000, slen/1000, qlen/1000)]
        bp_unit <- 'kb'
    }
    srank <- pafx[, list(qcovpct= max(qcovpct)), by= list(qname, sname)][, list(sname, srank= as.character(rank(-qcovpct))), by= qname]
    pafx <- merge(pafx, srank, by= c('sname', 'qname'), sort= FALSE)
    
    pafx[, strip := sprintf('%s | aln = %s | shown = %s', qname, N, shown)]
    pafx[, strip := factor(strip, levels= unique(strip))]
    
    pafx[, ysname := as.numeric(factor(sname, levels= unique(sname))), by= qname]
    pafx[, ysname := ysname/max(.SD$ysname), by= qname]
    
    maxsname <- max(pafx[, list(N= length(unique(sname))), by= qname]$N)
    qyend <- 0.05
    if(1 / maxsname <= qyend){
        qyend <- 1 / (maxsname * 2) 
    }

    gg <- ggplot(data= pafx) +
        scale_x_continuous(expand= c(0.01, 0.01)) +
        scale_y_continuous(expand= c(0, 0), labels= NULL, limits= c(-qyend/2, 1 + qyend)) +
        coord_cartesian(clip= 'off') +
        geom_text(data= unique(pafx[, list(strip, sname, ysname)]), aes(x= 0, y= ysname, label= paste0(sname, '  ')), size= 2.5, hjust= 1)

    if(nrow(pafx[mapq >= xargs$mapq]) > 0){
        gg <- gg + geom_segment(data= pafx[mapq >= xargs$mapq], aes(x= qstart, xend= sstart, y= qyend/2, yend= ysname - qyend/2, colour= srank), linetype= 'dotted') +
                   geom_segment(data= pafx[mapq >= xargs$mapq], aes(x= qend, xend= send, y= qyend/2, yend= ysname - qyend/2, colour= srank), linetype= 'dotted')
    }

    if(xargs$nrow == "auto"){
        xnrow <- ceiling(length(unique(pafx$qname)) / xargs$ncol)
    } else {
        xnrow <- as.integer(xargs$nrow)
    }

    gg <- gg +
        # Draw full contigs 
        geom_rect(data= unique(pafx[, list(strip, qlen)]), aes(xmin= 0, xmax= qlen, ymin= -qyend/2, ymax= qyend/2), colour= NA, fill= 'grey50', alpha= 0.5) +
        geom_rect(data= unique(pafx[, list(strip, slen, ysname)]), aes(xmin= 0, xmax= slen, ymin= ysname - qyend/2, ymax= ysname + qyend/2), fill= 'grey50', colour= NA, alpha= 0.5) +
        
        geom_rect(aes(xmin= sstart, xmax= send, ymin= ysname - qyend/2 , ymax= ysname + qyend/2), fill= 'white', colour= NA) +
        geom_rect(aes(xmin= sstart, xmax= send, ymin= ysname - qyend/2, ymax= ysname + qyend/2, fill= srank), alpha= 0.5, colour= NA) +
        
        geom_rect(aes(xmin= qstart, xmax= qend, ymin= -qyend/2, ymax= qyend/2), fill= 'white', colour= NA) +
        geom_rect(aes(xmin= qstart, xmax= qend, ymin= -qyend/2, ymax= qyend/2, fill= srank), alpha= 0.5, colour= NA) +

        geom_text(aes(x= sstart + (send - sstart)/2, y= ysname, label= ifelse(strand == '+', '>', '<'))) +
        facet_wrap_paginate(~ strip, scales= 'free', ncol= xargs$ncol, nrow= xnrow, page= 1) +
        xlab(sprintf('Contig size %s', bp_unit)) +
        ylab('') +
        theme(
            panel.grid.major.y= element_blank(),
            panel.grid.major.x= element_blank(),
            panel.background= element_rect('white'), 
            axis.line = element_line(colour = "black"),
            axis.title.y= element_blank(),
            axis.line.y= element_blank(),
            axis.line.x= element_line(size= 0.25, colour= 'grey50'),
            axis.ticks.y= element_blank(),
            legend.position= 'none', 
            plot.margin= unit(c(0.5, 0.5, 0.5, max(nchar(pafx$sname)) / 3), units= 'char'),
        ) 
    
    if(! is.null(xargs$panel_spacing)){
        gg <- gg + theme(panel.spacing.y = unit(xargs$panel_spacing, 'lines'))
    }

    n <- n_pages(gg)
    
    pdf(xargs$output, w= pagesize[1], h= pagesize[2])
    for(i in 1:n){
        cat(sprintf('Printing page %s of %s\n', i, n))
        print(gg + facet_wrap_paginate(~ strip, scales= 'free', ncol= xargs$ncol, nrow= xnrow, page= i)) 
    }
    invisible(dev.off())
}

