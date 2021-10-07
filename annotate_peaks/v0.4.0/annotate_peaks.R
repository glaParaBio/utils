#!/usr/bin/env Rscript

# New in v0.3.0: 
# * Accepts annotation in GTF format. We intentionally keep GFF as the only documented
#   one. 

suppressWarnings(library(data.table))
suppressWarnings(library(argparse))

options(scipen= 9)

parser <- ArgumentParser(description= 'Annotate peaks with genes in their vicinity')
parser$add_argument('--peaks', '-p', help= 'Peaks to annotate in bed-compatible format [%(default)s]', default= "-", metavar= 'FILE')
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
parser$add_argument('--summit', '-sm', help= "Column name or column index in the peak file giving the distance of the peak summit from the peak start. If set to '', use the peak mid-point", default= def,  metavar= sprintf('[%s]', def))

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

isGzip <- function(path){
    # Adapted from https://stackoverflow.com/questions/29493302/how-to-check-if-a-file-is-compressed-in-r
    f <- file(path)
    ext <- summary(f)$class
    close.connection(f)
    if(ext == 'gzfile') {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

prepare_tss <- function(gff, feature_type, gene_key, extra, verbose) {
    if(gff != '-' & isGzip(gff) == TRUE) {
        xcat <- 'gzip -cd'
    } else {
        xcat <- 'cat'
    }
    xgff <- fread(cmd= sprintf('%s %s | grep -v "^#" | awk -v FS="\t" \'$3 == "%s"\'', xcat, gff, feature_type), header= FALSE, sep= '\t')
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

    # setnames(xgff, names(xgff), paste0('', names(xgff))) 
    setnames(xgff, names(xgff)[1], paste0('#', names(xgff)[1])) # Comment out header line
    xgff <- unique(xgff[order(`#chrom`, tss_start, tss_end)])
    
    if(verbose) {
        write(sprintf('Found %s features\n', nrow(xgff)), stderr())
    }
    stopifnot(length(xgff$tss_id) == length(unique(xgff$tss_id)))
    xgff <- xgff[order(`#chrom`, tss_start)]
    return(xgff)
}

closest <- function(peaks, tss, summit, strand= '.', verbose= FALSE, tmpdir) {
    stopifnot(strand %in% c('+', '-', '.'))
    if(strand != '.') {
        tss <- tss[tss_strand == strand]
        ignore <- ' -id '
    } else {
        ignore <- ''
    }

    peaks <- copy(peaks)
    setcolorder(peaks, c('#chrom', 'summit_start', 'summit_end'))

    summit_file <- file.path(tmpdir, 'summit.bed') 
    write.table(peaks, summit_file, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- file.path(tmpdir, 'tss.bed')
    write.table(tss, tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <- file.path(tmpdir, 'intx_closest.bed') 

    cmd <- sprintf("set -e \
    closestBed %s -a %s -b %s -D b > %s", ignore, summit_file, tss_file, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    if(exit != 0) {
        stop(sprintf('Executing %s', cmd))
    }
    intx <- fread(out_intx, header= FALSE, col.names= c(names(peaks), names(tss), 'tss_distance'), sep= '\t')
    intx[, `#chrom` := NULL]
    intx[, unused := NULL]
    intx <- intx[tss_start >= 0] # Remove peaks with no feature assigned (e.g. on the edge of chroms)

    return(intx) 
}

intersection <- function(peaks, tss, span, verbose, tmpdir) {

    peaks <- copy(peaks)

    peaks[, extended_start := summit_start - span]
    peaks[, extended_start := ifelse(extended_start < 0, 0, extended_start)]
    peaks[, extended_end := summit_end + span]
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
    
    tss_names <- names(tss)
    tss_names[which(tss_names == '#chrom')] <- '#tss_chrom'
    
    hdr <- c(names(peaks), tss_names)

    if(file.size(out_intx) > 0) {
        intx <- fread(out_intx, header= FALSE, sep= '\t')
    } else {
        intx <- list() 
        for(x in hdr) {
            intx[[x]] <- NA
        }
        intx <- as.data.table(intx)[0, ]
    }
    setnames(intx, names(intx), hdr)
    intx[, `#tss_chrom` := NULL]
    stopifnot(length(names(intx)) == length(unique(names(intx))))

    intx[, extended_start := NULL]
    intx[, extended_end := NULL]
    peaks[, extended_start := NULL]
    peaks[, extended_end := NULL]

    stopifnot(intx$tss_strand %in% c('+', '-'))

    intx[, tss_distance := ifelse(tss_strand == "+", summit_start - tss_start,
                                  tss_start - summit_start)]
    
    intx[, unused := NULL]
    
    intx <- intx[order(`#chrom`, summit_start, summit_end, abs(tss_distance))]
    setcolorder(intx, c('#chrom', 'summit_start', 'summit_end')) 
    if(verbose) {
        write(sprintf('Intersection returned %s features\n', nrow(intx)), stderr())
    }
    return(intx)
}

intervening_tss <- function(annotated_peaks, tss, verbose, tmpdir= tmpdir) {
    # We need to capture the TSSs within the most extreme TSSs already in output.
    # Create intervals where each peak_id goes from the leftmost TSS to rightmost TSS
    extended_peaks <- annotated_peaks[, list(summit_start= min(tss_start), summit_end= max(tss_end)), by= list(`#chrom`, id__, xsummit_start= summit_start, xsummit_end= summit_end)]
    setcolorder(extended_peaks, c('#chrom', 'summit_start', 'summit_end'))
    
    intervene <- intersection(extended_peaks, tss, span= 0, verbose= verbose, tmpdir)

    # Replace extended coords with original summit_start/end
    intervene[, summit_start := xsummit_start]
    intervene[, summit_end := xsummit_end]
    intervene[, xsummit_start := NULL]
    intervene[, xsummit_end := NULL]

    # Recalculate distance
    intervene[, tss_distance := ifelse(tss_strand == "+", summit_start - tss_start,
                                  tss_start - summit_start)]

    return(intervene)
}

get_tmpdir <- function() {
    tmpdir <- tempfile(pattern= paste0('tmp_annotate_peaks.', datestr()), tmpdir= '.')
    dir.create(tmpdir)
    return(tmpdir)
}

prepare_peaks <- function(peaks, summit) {
    # Prepare the working bed file
    wrk_peaks <- copy(peaks[, 1:3])
    wrk_peaks[, id__ := peaks$id__]
    setnames(wrk_peaks, 1:3, c('#chrom', 'peak_start', 'peak_end'))

    # Add summit column
    if(grepl("^[[:digit:]]", summit) & ! summit %in% names(peaks)) {
        summit <- paste0('V', summit)
    }
    if(summit %in% names(peaks)) {
        peaks[[summit]] <- as.integer(peaks[[summit]])
        wrk_peaks[, summit := peaks[[summit]]]
        if(any(is.na(wrk_peaks$summit))) {
            stop(sprintf('Column %s contains non-integer values', summit))
        }
    } else if(summit == '') {
        offset <- round((wrk_peaks$peak_end - wrk_peaks$peak_start)/2)
        offset <- ifelse(offset == 0, 1, offset)
        offset <- as.integer(offset)
        wrk_peaks[, summit := offset]
    } else {
        msg <- sprintf('Summit column "%s" not found. Perhaps you need to adjust the argument to option --summit?', summit)
        stop(msg)
    }
    wrk_peaks[, summit_start := peak_start + summit - 1]
    wrk_peaks[, summit_end := summit_start + 1]
    # Keep only necessary columns
    wrk_peaks[, peak_start := NULL]
    wrk_peaks[, peak_end := NULL]
    wrk_peaks[, summit := NULL]
    setcolorder(wrk_peaks, c('#chrom', 'summit_start', 'summit_end'))
    wrk_peaks <- unique(wrk_peaks)
    wrk_peaks <- wrk_peaks[order(`#chrom`, summit_start, summit_end)]
    return(wrk_peaks)
}

# -----------------------

# Read peak file and add id__ column
if(args$peaks == '-') {
    peaks <- fread('file:///dev/stdin', sep= '\t')
} else {
    peaks <- fread(args$peaks, sep= '\t')
}
stopifnot(!grepl('id__', names(peaks)))
peaks[, 2] <- as.integer(unlist(peaks[, 2]))
peaks[, 3] <- as.integer(unlist(peaks[, 3]))
peaks[, id__ := sprintf('%s_%s_%s', peaks[[1]], peaks[[2]], peaks[[3]])]

wrk_peaks <- prepare_peaks(peaks, args$summit)


tss <- prepare_tss(args$gff, args$feature_type, args$gene_key, args$extra, args$verbose)

tryCatch({
        tmpdir <- get_tmpdir()
        
        region_intx <- intersection(wrk_peaks, tss, args$span, verbose= args$verbose, tmpdir= tmpdir)

        xclosest <- closest(wrk_peaks, tss, verbose= args$verbose, tmpdir= tmpdir)
        xclosest_plus <- closest(wrk_peaks, tss, strand= '+' ,verbose= args$verbose, tmpdir= tmpdir)
        xclosest_minus <- closest(wrk_peaks, tss, strand= '-' ,verbose= args$verbose, tmpdir= tmpdir)

        xclosest <- unique(rbindlist(list(xclosest, xclosest_plus, xclosest_minus, region_intx), use.names= TRUE))

        intervene <- intervening_tss(xclosest, tss, verbose= args$verbose, tmpdir= tmpdir)
    }, finally= {
        unlink(tmpdir, recursive= TRUE)
    })

xclosest <- unique(rbind(xclosest, intervene))
xclosest[, tss_start := NULL]
xclosest[, tss_end := NULL]

xclosest <- xclosest[order(id__, abs(tss_distance))]
xclosest[, tss_distance_rank := 1:nrow(.SD), by= id__]

xclosest[, `#chrom` := NULL]
xclosest[, summit_start := NULL]
xclosest[, summit_end := NULL]

for(x in names(xclosest)) {
    if(x == "id__") {
        next
    }
    peaks[[x]] <- NULL
}
peaks <- unique(peaks)
intx <- merge(peaks, xclosest, by= 'id__', sort= FALSE, all.x= TRUE)
intx <- intx[order(intx[[1]], intx[[2]], intx[[3]], tss_distance_rank)]
intx[, id__ := NULL]
peaks[, id__ := NULL]

# Sanity checks
stopifnot(nrow(peaks) == nrow(intx[tss_distance_rank == 1 | is.na(tss_distance_rank)]))
stopifnot(identical(sort(peaks[[2]]), sort(intx[tss_distance_rank == 1 | is.na(tss_distance_rank)][[2]])))
stopifnot(names(peaks) %in% names(intx))

if(!grepl('^#', names(intx)[1])) {
    setnames(intx, names(intx)[1], paste0('#', names(intx)[1]))
}
write.table(intx, file= stdout(), row.names= FALSE, sep= '\t', quote= FALSE)

quit()

