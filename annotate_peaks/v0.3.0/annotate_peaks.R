#!/usr/bin/env Rscript

# New in v0.3.0: 
# * Accepts annotation in GTF format. We intentionally keep GFF as the only documented
#   one. 

suppressWarnings(library(data.table))
suppressWarnings(library(argparse))

parser <- ArgumentParser(description= 'Annotate peaks with genes in their vicinity')
parser$add_argument('--peaks', '-p', help= 'Peaks to annotate in bed format. It must have a header line of column names', required= TRUE, metavar= 'FILE')
parser$add_argument('--gff', '-gff', help= 'Annotation file in GFF format', required= TRUE, metavar= 'FILE')
parser$add_argument('--genome', '-g', help= 'Tab separated file of chromsome names and sizes. The fasta index (.fai file) of the reference genome is suitable)', required= TRUE, metavar= 'FILE')

def <- 'mRNA'
parser$add_argument('--feature-type', '-f', help= 'Feature type (column 3) in GFF identifying genes or transcript. Typically use "gene" or "mRNA"', default= def, metavar= sprintf('[%s]', def))

def <- 'ID'
parser$add_argument('--gene-key', '-k', help= 'Attribute key in the gff file of the feature to use for annotation. Typically "ID" or "Parent"', default= def, metavar= sprintf('[%s]', def))

def <- c('description')
parser$add_argument('--extra', '-x', help= sprintf('Extra attributes to add to the annotation [%s]', paste(def, collapse= ' ')), default= def, metavar= 'NAME', nargs= '*')

def <- 10000
parser$add_argument('--span', '-s', type= 'integer', help= 'Extend each peak by this many bases and assign genes whose TSS is in the intersection', default= def, metavar= sprintf('[%s]', def))

def <- 'summit'
parser$add_argument('--summit', '-sm', help= "Column name in the peak file giving the distance of the peak summit from the peak start. If set to '', use the peak mid-point", default= def,  metavar= sprintf('[%s]', def))

def <- -1
parser$add_argument('--max-rank', '-r', help= 'Keep at most max-rank features closest to each peak. No limit if < 0', default= def,  metavar= sprintf('[%s]', def))

parser$add_argument('--verbose', '-V', action= 'store_true', help= 'Verbose mode mostly for debugging')
parser$add_argument('--version', '-v', action= 'version', version= '0.3.0')

args <- parser$parse_args()

datestr <- function() {
    tm <- as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S")
    str <- strftime(tm , "%Y-%m-%dT%H:%M:%S")
    return(str)
}

get_gff_attribute <- function(gff_attr, key) {
    values <- rep(NA, length(gff_attr))
    for(i in 1:length(gff_attr)) {
        attr_list <- strsplit(gff_attr[i], ';')[[1]]
        attr_list <- sapply(attr_list, sub, pattern= '^ ', replacement= '', USE.NAMES= FALSE)
        keep <- grepl(sprintf('^%s( |=)', key), attr_list)
        if(sum(keep) > 1) {
            stop(sprintf('Attribute "%s" found more than once in\n%s', key, gff_attr[i]))
        }
        if(sum(keep) == 0) {
            value <- NA
        } else {
            value <- attr_list[keep]
            value <- sub(sprintf('^%s( |=)', key), '', value)
            value <- gsub('^"|"$', '', value)
        }
        values[i] <- value
    }
    return(values)
}

prepare_tss <- function(gff, feature_type, gene_key, extra, verbose) {
    xgff <- fread(cmd= sprintf('grep -v "^#" %s | awk -v FS="\t" \'$3 == "%s"\'', gff, feature_type), header= FALSE)
    if(nrow(xgff) == 0) {
        stop(sprintf('There are no records of type "%s" in file "%s"', feature_type, gff))
    }
    setnames(xgff, c('V1', 'V4', 'V5', 'V7', 'V9'), c('chrom', 'start', 'end', 'strand', 'gff_attr'))
    stopifnot(xgff$strand %in% c('+', '-'))

    xgff[, tss_start := ifelse(strand == '+', start - 1, end - 1)]
    
    xgff[, tss_id := get_gff_attribute(gff_attr, gene_key)]
    if(sum(is.na(xgff$tss_id)) > 0) {
        stop(sprintf('Some GFF/GTF records of type "%s" do not have attribute "%s"', feature_type, gene_key))    
    }

    if(gene_key %in% extra) {
        extra <- extra[extra != gene_key]
    }
    if(length(extra) == 0 | extra[1] == '') {
        extra <- c()
    }
    for(xt in extra) {
        xgff[, tag__ := get_gff_attribute(gff_attr, xt)]
        xgff [, tag__ := sapply(tag__, URLdecode)]
        setnames(xgff, 'tag__', xt)
    }
    outnames <- c('chrom', 'tss_start', 'tss_id', 'strand', extra)
    xgff <- xgff[, outnames, with= FALSE]
    tss <- rbind(
            xgff[strand == '+'][, list(tss_start= min(tss_start)), by= list(chrom, tss_id)],
            xgff[strand == '-'][, list(tss_start= max(tss_start)), by= list(chrom, tss_id)]
            )
    xgff <- merge(xgff, tss, by= c('chrom', 'tss_id', 'tss_start'), sort= FALSE)
    xgff$tss_end <- xgff$tss_start + 1
    xgff$unused <- '.'
    xgff <- xgff[, c('chrom', 'tss_start', 'tss_end', 'tss_id', 'unused', 'strand', extra), with= FALSE]

    # Prefix every column with a string that will make very unlikely to have
    # duplicate columns with the peak file (not great but good enough). 
    setnames(xgff, names(xgff), paste0('gff__', names(xgff))) 
    setnames(xgff, names(xgff)[1], paste0('#', names(xgff)[1])) # Comment out header line
    xgff <- unique(xgff[order(`#gff__chrom`, gff__tss_start, gff__tss_end)])
    
    if(verbose) {
        write(sprintf('Found %s features\n', nrow(xgff)), stderr())
    }
    stopifnot(length(xgff$gff__tss_id) == length(unique(xgff$gff__tss_id)))
    return(xgff)
}

closest <- function(peak_file, tss, summit, verbose= FALSE) {
    # Prepare peak file. Get the summit 
    peaks <- fread(peak_file)
    peaks[, id__ := 1:nrow(peaks)]

    if(summit != '') {
        offset <- peaks[[summit]]
    } else {
        offset <- round((peaks[,3] - peaks[,2])/2)[[1]]
        offset <- ifelse(offset == 0, 1, offset)
    }
    peak_summit <- peaks[, c(1,2)]
    setnames(peak_summit, names(peak_summit), c('#chrom', 'start'))
    peak_summit[, ref_start := start + offset]
    peak_summit[, ref_end := ref_start + 1]
    peak_summit[, start := NULL]
    peak_summit[, id__ := 1:nrow(peaks)]

    peak_file2 <- tempfile(pattern= paste0(basename(peak_file), '.', datestr(), '.'), tmpdir= '.', fileext= '.tmp.bed')
    write.table(peak_summit, peak_file2, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- tempfile(pattern= paste0(basename(peak_file), '.', datestr(), '.'), tmpdir= '.', fileext= '.tmp.bed')
    write.table(tss, tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <-  tempfile(pattern= paste0(basename(peak_file), '.', datestr(), '.'), tmpdir= '.', fileext= '.tmp.closest.bed')

    cmd <- sprintf("set -e \
    set -o pipefail \
    closestBed -a %s -b %s -D b > %s", peak_file2, tss_file, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    x <- file.remove(c(tss_file, peak_file2))
    if(exit != 0) {
        x <- file.remove(out_intx)
        stop(sprintf('Executing %s', cmd))
    }
    intx <- fread(out_intx, header= FALSE, col.names= c(names(peak_summit), names(tss), 'tss_distance'))
    x <- file.remove(out_intx)
    intx[, `#chrom` := NULL]
    intx[, ref_start := NULL]
    intx[, ref_end := NULL]
    intx[, `#gff__chrom` := NULL]
    intx[, gff__tss_start := NULL]
    intx[, gff__tss_end := NULL]
    intx[, gff__unused := NULL]

    intx <- merge(peaks, intx, by= 'id__')
    intx[, 'id__' := NULL]
    return(intx) 
}

intersection <- function(peak_file, genome_file, tss, span, summit, verbose) {

    peaks <- fread(peak_file)
    peaks[, start__ := peaks[, 2]]
    peaks[, end__ := peaks[, 3]]
    peak_file2 <- tempfile(pattern= paste0(basename(peak_file), '.', datestr(), '.'), tmpdir= '.', fileext= '.tmp.bed')
    write.table(peaks, peak_file2, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- tempfile(pattern= paste0(basename(peak_file), '.', datestr(), '.'), tmpdir= '.', fileext= '.tmp.bed')
    write.table(tss, tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <-  tempfile(pattern= paste0(basename(peak_file), '.', datestr(), '.'), tmpdir= '.', fileext= '.tmp.tsv')

    cmd <- sprintf("set -e \
    set -o pipefail \
    slopBed -b %s -i %s -g %s \\
    | intersectBed -header -a - -wa -wb -b %s > %s", span, peak_file2, genome_file, tss_file, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    x <- file.remove(c(tss_file, peak_file2))
    if(exit != 0) {
        file.remove(out_intx)
        stop(sprintf('Executing %s', cmd))
    }
    
    hdr <- c(names(peaks), names(tss))
    stopifnot(length(hdr) == length(unique(hdr)))

    intx <- fread(out_intx, header= FALSE, col.names= hdr)
    file.remove(out_intx)
    intx[, 2] <- intx$start__
    intx[, 3] <- intx$end__
    intx[, start__ := NULL]
    intx[, end__ := NULL]
    peaks[, start__ := NULL]
    peaks[, end__ := NULL]

    if(summit != '') {
        if(!summit %in% names(intx)) {
            stop(sprintf('Summit column "%s" not found. Perhaps you need to adjust the argument to option --summit?', summit))
        }
        offset <- intx[[summit]]
    } else {
        offset <- round((intx[,3] - intx[,2])/2)[[1]]
        offset <- ifelse(offset == 0, 1, offset)
    }
    intx$peak__refpoint <- intx[,2] + offset 

    stopifnot(intx$gff__strand %in% c('+', '-'))

    intx[, tss_distance := ifelse(gff__strand == "+", peak__refpoint - gff__tss_start,
                                  gff__tss_start - peak__refpoint)]
    
    intx[, `#gff__chrom` := NULL]
    intx[, gff__tss_start := NULL]
    intx[, gff__tss_end := NULL]
    intx[, peak__refpoint := NULL]
    intx[, gff__unused := NULL]
    
    hdr <- names(intx)[1:3]
    setnames(intx, 1:3, c('chrom', 'start', 'end'))

    intx <- intx[order(chrom, start, end, abs(tss_distance))]
    setnames(intx, 1:3, hdr)
    
    if(verbose) {
        write(sprintf('Intersection returned %s features\n', nrow(intx)), stderr())
    }
    return(intx)
}

merge_intx <- function(region_intx, xclosest, max_rank) {
    stopifnot(identical(names(region_intx), names(xclosest)))
    hdr <- names(region_intx)[1:3]
    setnames(region_intx, 1:3, c('chrom', 'start', 'end'))
    setnames(xclosest, 1:3, c('chrom', 'start', 'end'))
    xclosest[, id__ := paste(chrom, start, end, sep= '_')]
    region_intx[, id__ := paste(chrom, start, end, sep= '_')]
    stopifnot(region_intx$id__ %in% xclosest$id__)
    miss <- xclosest[! id__ %in% region_intx$id__]
    intx <- rbind(region_intx, miss)

    dist_rank <- unique(intx[, list(id__, gff__tss_id, tss_distance)])
    dist_rank <- dist_rank[order(id__, abs(tss_distance))]
    dist_rank[, tss_distance_rank := 1:nrow(.SD), by= id__]

    intx_hdr <- names(intx)
    intx <- merge(intx, dist_rank, by= c('id__', 'gff__tss_id', 'tss_distance'))
    setcolorder(intx, intx_hdr)
    intx[, id__ := NULL]
    intx <- intx[order(chrom, start, end, tss_distance_rank)]
    setnames(intx, 1:length(hdr), hdr)
    setnames(intx, names(intx), sub('^gff__', '', names(intx)))

    if(max_rank > 0) {
        intx <- intx[tss_distance_rank <= max_rank]
    }

    return(intx)
}

# -----------------------

tss <- prepare_tss(args$gff, args$feature_type, args$gene_key, args$extra, args$verbose)

region_intx <- intersection(args$peaks, args$genome, tss, args$span, summit= args$summit, verbose= args$verbose)

xclosest <- closest(args$peaks, tss, args$summit, args$verbose)

intx <- merge_intx(region_intx, xclosest, args$max_rank)

if(identical(names(intx)[6], 'strand') & (all(intx$strand == '.') | all(is.na(intx$strand)))) {
    # If the input peak file contains a "strand" column in position 6 and this
    # column is only missing values, replace it with the GFF strand
    intx[, strand := NULL]
    setcolorder(intx, c(1:5, which(names(intx) == 'strand')))
}

# Sanity checks
peaks <- fread(args$peaks)
stopifnot(nrow(peaks) == nrow(intx[tss_distance_rank == 1]))
stopifnot(identical(sort(peaks[[2]]), sort(intx[tss_distance_rank == 1][[2]])))
stopifnot(names(peaks) %in% names(intx))

write.table(intx, file= stdout(), row.names= FALSE, sep= '\t', quote= FALSE)
quit()
