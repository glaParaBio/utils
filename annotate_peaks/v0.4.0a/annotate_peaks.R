#!/usr/bin/env Rscript

# New in v0.3.0: 
# * Accepts annotation in GTF format. We intentionally keep GFF as the only documented
#   one. 

suppressWarnings(library(data.table))
suppressWarnings(library(argparse))

parser <- ArgumentParser(description= 'Annotate peaks with genes in their vicinity')
parser$add_argument('--peaks', '-p', help= 'Peaks to annotate in bed-compatible format', required= TRUE, metavar= 'FILE')
parser$add_argument('--gff', '-gff', help= 'Annotation file in GFF format', required= TRUE, metavar= 'FILE')

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
parser$add_argument('--version', '-v', action= 'version', version= '0.4.0')

args <- parser$parse_args()

datestr <- function() {
    tm <- as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S")
    str <- strftime(tm , "%Y-%m-%dT%H_%M_%S")
    return(str)
}

peak_reader <- function(peak_file) {
    peaks <- fread(peak_file, sep= '\t')
    if(!grepl('^#', names(peaks)[1])) {
        setnames(peaks, names(peaks)[1], paste0('#', names(peaks)[1]))
    }
    return(peaks)
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
    xgff <- fread(cmd= sprintf('grep -v "^#" %s | awk -v FS="\t" \'$3 == "%s"\'', gff, feature_type), header= FALSE, sep= '\t')
    if(nrow(xgff) == 0) {
        stop(sprintf('There are no records of type "%s" in file "%s"', feature_type, gff))
    }
    setnames(xgff, c('V1', 'V4', 'V5', 'V7', 'V9'), c('chrom', 'start', 'end', 'tss_strand', 'gff_attr'))
    stopifnot(xgff$tss_strand %in% c('+', '-'))

    xgff[, tss_start := ifelse(tss_strand == '+', start - 1, end - 1)]
    
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
        xgff [, tag__ := sapply(tag__, function(x) ifelse(is.na(x), NA, URLdecode(x)))]
        setnames(xgff, 'tag__', xt)
    }
    outnames <- c('chrom', 'tss_start', 'tss_id', 'tss_strand', extra)
    xgff <- xgff[, outnames, with= FALSE]
    tss <- rbind(
            xgff[tss_strand == '+'][, list(tss_start= min(tss_start)), by= list(chrom, tss_id)],
            xgff[tss_strand == '-'][, list(tss_start= max(tss_start)), by= list(chrom, tss_id)]
            )
    xgff <- merge(xgff, tss, by= c('chrom', 'tss_id', 'tss_start'), sort= FALSE)
    xgff$tss_end <- xgff$tss_start + 1
    xgff$unused <- '.'
    xgff <- xgff[, c('chrom', 'tss_start', 'tss_end', 'tss_id', 'unused', 'tss_strand', extra), with= FALSE]

    # Prefix every column with a string that will make very unlikely to have
    # duplicate columns with the peak file (not great but good enough). 
    setnames(xgff, names(xgff), paste0('gff__', names(xgff))) 
    setnames(xgff, names(xgff)[1], paste0('#', names(xgff)[1])) # Comment out header line
    xgff <- unique(xgff[order(`#gff__chrom`, gff__tss_start, gff__tss_end)])
    
    if(verbose) {
        write(sprintf('Found %s features\n', nrow(xgff)), stderr())
    }
    stopifnot(length(xgff$gff__tss_id) == length(unique(xgff$gff__tss_id)))
    xgff <- xgff[order(`#gff__chrom`, gff__tss_start)]
    return(xgff)
}

closest <- function(peaks, tss, summit, tss_strand= '.', verbose= FALSE, tmpdir) {
    stopifnot(tss_strand %in% c('+', '-', '.'))
    if(tss_strand != '.') {
        tss <- tss[gff__tss_strand == tss_strand]
        ignore <- ' -id '
    } else {
        ignore <- ''
    }

    peaks <- copy(peaks)
    # Prepare peak file. Get the summit 
    #peaks <- peak_reader(peak_file)
    #peaks[, id__ := 1:nrow(peaks)]

    #if(summit != '') {
    #    offset <- peaks[[summit]]
    #} else {
    #    offset <- round((peaks[,3] - peaks[,2])/2)[[1]]
    #    offset <- ifelse(offset == 0, 1, offset)
    #}
    #peak_summit <- peaks[, list()]
    #setnames(peak_summit, names(peak_summit), c('#chrom', 'start'))
    #peak_summit[, ref_start := start + offset]
    #peak_summit[, ref_end := ref_start + 1]
    #peak_summit[, start := NULL]
    #peak_summit[, id__ := 1:nrow(peaks)]
    setcolorder(peaks, c('#chrom', 'summit_start', 'summit_end'))

    summit_file <- file.path(tmpdir, 'summit.bed') 
    write.table(peaks, summit_file, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- file.path(tmpdir, 'tss.bed')
    write.table(tss, tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <- file.path(tmpdir, 'intx_closest.bed') 

    cmd <- sprintf("set -e \
    set -o pipefail \
    closestBed %s -a %s -b %s -D b > %s", ignore, summit_file, tss_file, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    if(exit != 0) {
        stop(sprintf('Executing %s', cmd))
    }
    intx <- fread(out_intx, header= FALSE, col.names= c(names(peaks), names(tss), 'tss_distance'), sep= '\t')
    intx[, `#gff__chrom` := NULL]
    #intx[, summit_start := NULL]
    #intx[, summit_end := NULL]
    #intx[, gff__tss_start := NULL]
    #intx[, gff__tss_end := NULL]
    intx[, gff__unused := NULL]

    return(intx) 
}

intersection <- function(peaks, tss, span, verbose, tmpdir) {

    peaks <- copy(peaks)

    peaks[, extended_start := peak_start - span]
    peaks[, extended_start := ifelse(extended_start < 0, 0, extended_start)]
    peaks[, extended_end := peak_end + span]
    setcolorder(peaks, c('#chrom', 'extended_start', 'extended_end'))

    peak_file <- file.path(tmpdir, 'peaks.bed') 
    write.table(peaks, peak_file, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- file.path(tmpdir, 'tss.bed')
    write.table(tss, tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <-  file.path(tmpdir, 'extended_intersection.bed')

    cmd <- sprintf("set -e \
    intersectBed -a %s -wa -wb -b %s > %s", peak_file, tss_file, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    if(exit != 0) {
        stop(sprintf('Executing %s', cmd))
    }

    hdr <- c(names(peaks), names(tss))
    stopifnot(length(hdr) == length(unique(hdr)))

    if(file.size(out_intx) == 0) {
        return(NA)
    }
    intx <- fread(out_intx, header= FALSE, sep= '\t')
    setnames(intx, names(intx), hdr)

    intx[, extended_start := NULL]
    intx[, extended_end := NULL]
    peaks[, extended_start := NULL]
    peaks[, extended_end := NULL]

    stopifnot(intx$gff__tss_strand %in% c('+', '-'))

    intx[, tss_distance := ifelse(gff__tss_strand == "+", summit_start - gff__tss_start,
                                  gff__tss_start - summit_start)]
    
    intx[, `#gff__chrom` := NULL]
    intx[, gff__tss_start := NULL]
    intx[, gff__tss_end := NULL]
    intx[, gff__unused := NULL]
    
    intx <- intx[order(`#chrom`, peak_start, peak_end, abs(tss_distance))]
    
    if(verbose) {
        write(sprintf('Intersection returned %s features\n', nrow(intx)), stderr())
    }
    return(intx)
}

merge_intx <- function(region_intx, xclosest, max_rank) {
    hdr <- names(xclosest)[1:3]
    setnames(xclosest, 1:3, c('chrom', 'start', 'end'))
    xclosest[, id__ := paste(chrom, start, end, sep= '_')]

    if(is.data.table(region_intx)) {
        setnames(region_intx, 1:3, c('chrom', 'start', 'end'))
        region_intx[, id__ := paste(chrom, start, end, sep= '_')]
        stopifnot(identical(names(region_intx), names(xclosest)))
        stopifnot(region_intx$id__ %in% xclosest$id__)
        intx <- unique(rbind(region_intx, xclosest))
    } else {
        intx <- xclosest
    }

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

intervening_tss <- function(peaks, annotated_peaks, tss, verbose, tmpdir= tmpdir) {
    # We need to capture the TSSs within the most extreme TSSs already in output.
    annotated_peaks <- copy(annotated_peaks)
    annotated_peaks <- merge(annotated_peaks, tss[, list(gff__tss_start, gff__tss_end, gff__tss_id)], by= 'gff__tss_id')
    annotated_peaks <- annotated_peaks[, list(peak_start= min(gff__tss_start) - 1, peak_end= max(gff__tss_end)), by= list(`#chrom`, x_= peak_start, y_= peak_end)]

    xpeaks <- merge(peaks, annotated_peaks, by.x= c('#chrom', 'peak_start', 'peak_end'), by.y= c('#chrom', 'x_', 'y_'))
    xpeaks[,2] <- xpeaks$peak_start
    xpeaks[,3] <- xpeaks$peak_end
    xpeaks[, peak_start__ := NULL]
    xpeaks[, peak_end__ := NULL]
    xpeaks[['summit__']] <- xpeaks$summit_pos__ - as.integer(xpeaks[,2])

    xpeaks[, summit_pos__ := NULL]
    peaks[, summit_pos__ := NULL]

    intervene <- intersection(xpeaks, tss, span= 0, verbose= verbose, tmpdir)
    for(x in names(xpeaks)) {
        if(x != 'row_id__') {
            intervene[[x]] <- NULL
        }
    }
    intervene <- merge(peaks, intervene, by= 'row_id__')
    intervene[, row_id__ := NULL]
    return(intervene)
}

get_tmpdir <- function() {
    tmpdir <- tempfile(pattern= paste0('tmp_annotate_peaks.', datestr()), tmpdir= '.')
    dir.create(tmpdir)
    return(tmpdir)
}

# -----------------------

# Read peak file and add id__ column
if(args$peaks == '-') {
    peaks <- fread('file:///dev/stdin', sep= '\t')
} else {
    peaks <- fread(args$peaks, sep= '\t')
}
stopifnot(!grepl('id__', names(peaks)))
peaks[, id__ := 1:nrow(peaks)]
peaks[, 2] <- as.integer(unlist(peaks[, 2]))
peaks[, 3] <- as.integer(unlist(peaks[, 3]))

# Prepare the working bed file
wrk_peaks <- copy(peaks[, 1:3])
wrk_peaks[, id__ := peaks$id__]
setnames(wrk_peaks, 1:3, c('#chrom', 'peak_start', 'peak_end'))

# Add summit column
if(args$summit %in% names(peaks)) {
    peaks[[args$summit]] <- as.integer(peaks[[args$summit]])
    wrk_peaks[, summit := peaks[[args$summit]]]
    if(any(is.na(wrk_peaks$summit))) {
        stop(sprintf('Column %s contains non-integer values', args$summit))
    }
} else if(args$summit == '') {
    offset <- round((wrk_peaks$peak_end - wrk_peaks$peak_start)/2)
    offset <- ifelse(offset == 0, 1, offset)
    offset <- as.integer(offset)
    wrk_peaks[, summit := offset]
} else {
    msg <- sprintf('Summit column "%s" not found. Perhaps you need to adjust the argument to option --summit?', args$summit)
    stop(msg)
}
wrk_peaks[, summit_start := peak_start + summit]
wrk_peaks[, summit_end := summit_start + 1]
wrk_peaks <- wrk_peaks[order(`#chrom`, peak_start, peak_end)]

tmpdir <- get_tmpdir()

tss <- prepare_tss(args$gff, args$feature_type, args$gene_key, args$extra, args$verbose)

region_intx <- intersection(wrk_peaks, tss, args$span, verbose= args$verbose, tmpdir= tmpdir)

xclosest <- closest(wrk_peaks, tss, verbose= args$verbose, tmpdir= tmpdir)
xclosest_plus <- closest(wrk_peaks, tss, tss_strand= '+' ,verbose= args$verbose, tmpdir= tmpdir)
xclosest_minus <- closest(wrk_peaks, tss, tss_strand= '-' ,verbose= args$verbose, tmpdir= tmpdir)
xclosest <- unique(rbindlist(list(xclosest, xclosest_plus, xclosest_minus)))

print(xclosest)
quit()

intervene <- intervening_tss(wrk_peaks, xclosest, tss, verbose= args$verbose, tmpdir= tmpdir)

xclosest <- unique(rbind(xclosest, intervene))

intx <- merge_intx(region_intx, xclosest, args$max_rank)

if(identical(names(intx)[6], 'tss_strand') & (all(intx$tss_strand == '.') | all(is.na(intx$tss_strand)))) {
    # If the input peak file contains a "tss_strand" column in position 6 and this
    # column is only missing values, replace it with the GFF tss_strand
    intx[, tss_strand := NULL]
    setcolorder(intx, c(1:5, which(names(intx) == 'tss_strand')))
}

# Sanity checks
stopifnot(nrow(peaks) == nrow(intx[tss_distance_rank == 1]))
stopifnot(identical(sort(peaks[[2]]), sort(intx[tss_distance_rank == 1][[2]])))
stopifnot(names(peaks) %in% names(intx))

write.table(intx, file= stdout(), row.names= FALSE, sep= '\t', quote= FALSE)
quit()
