#!/usr/bin/env Rscript

suppressWarnings(library(data.table))
suppressWarnings(library(argparse))

options(scipen= 9)

parser <- ArgumentParser(description= 'Assign to each feature in GFF or BED file the nearest peak upstream and the nearest downstream')
parser$add_argument('--peaks', '-p', help= 'Peak file in bed-compatible format [%(default)s]', default= "-", metavar= 'FILE')
parser$add_argument('--gff', '-gff', help= 'Annotation file in GFF or BED format [%(default)s]', default= '-', metavar= 'FILE')

def <- 'mRNA'
parser$add_argument('--feature-type', '-f', help= 'Feature type from 3rd column in GFF identifying features to be annotated. Typically use "gene" or "mRNA". Ignored with annotation in BED format [%(default)s]', default= def, metavar= sprintf('[%s]', def))

def <- 'ID'
parser$add_argument('--feature-name', '-fn', help= 'With annotation in GFF: Use this attribute key to name features (typically: ID, gene_id). With annotation in BED format: Column index with the feature names (typically column 4) [%(default)s]', default= def, metavar= sprintf('[%s]', def))

def <- 'summit'
parser$add_argument('--summit', '-sm', help= "Column name or column index in the peak file giving the distance of the peak summit from the peak start. If set to '', use the peak mid-point [%(default)s]", default= def,  metavar= sprintf('[%s]', def))

def <- 1
parser$add_argument('--n-closest', '-n', help= 'Number of closest peaks upstream and downstream to assign to each gene [%(default)s]', default= def,  metavar= sprintf('[%s]', def), type= 'integer')

parser$add_argument('--add-column-prefix', '-a', help= "Add this prefix to the column names of the annotation to avoid collision with the column names of the peak file [%(default)s]", default= 'feature_')

parser$add_argument('--verbose', '-V', action= 'store_true', help= 'Verbose mode mostly for debugging')
parser$add_argument('--version', '-v', action= 'version', version= '0.2.0')

args <- parser$parse_args()

datestr <- function() {
    tm <- as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S")
    str <- strftime(tm , "%Y-%m-%dT%H_%M_%S")
    return(str)
}

peak_reader <- function(peak_file) {
    if(peak_file == '-') {
        peaks <- fread('file:///dev/stdin')
    } else if(grepl('^/dev/fd/', peak_file) == TRUE) {
        fin <- fifo(peak_file, open = "r")
        xpeaks <- readLines(fin)
        close(fin)
        fn <- tempfile() 
        writeLines(xpeaks, fn)
        peaks <- fread(fn)
        unlink(fn)
    } else {
        if(isGzip(peak_file) == TRUE) {
            xcat <- 'gzip -cd'
        } else {
            xcat <- 'cat'
        }
        peaks <- fread(cmd= sprintf('%s %s', xcat, peak_file), sep= '\t')
    }
    if(ncol(peaks) == 0) {
        peaks <- data.table(`#chrom`= character(), start= integer(), end= integer())
    }
    peaks[[1]] <- as.character(peaks[[1]])
    peaks <- unique(peaks)
    peaks <- peaks[order(peaks[[1]], peaks[[2]], peaks[[3]])]
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

get_data_type <- function(bedOrGff) {
    if(ncol(bedOrGff) == 6){
        if(is.numeric(bedOrGff[[2]]) & is.numeric(bedOrGff[[3]])) {
            if(all(bedOrGff[[6]] %in% c('+', '-'))) {
                return('bed')
            } else {
                stop('Invalid type of BED: strand column (6th) must contain only "+", "-"')
            }
        } else {
            stop('Invalid type of BED: columns 2 and/or 3 are not numeric')
        }
    } else if(ncol(bedOrGff) == 9) {
        return('gff')
    } else {
        stop('Data type is neither BED6 or GFF') 
    }
}

convert_bed_to_gff <- function(bed, feature_type, feature_name_idx) {
    if(feature_name_idx <= 0 | feature_name_idx > ncol(bed)) {
        stop(sprintf('Invalid column index to extract feature names: %s', feature_name_idx))
    }

    gff <- data.table(
        V1= bed[[1]], 
        V2= '.', 
        V3= feature_type, 
        V4= bed[[2]] + 1, 
        V5= bed[[3]], 
        V6= '.', 
        V7= bed[[6]], 
        V8= '.', 
        V9= sprintf('%s=%s', feature_name_idx, bed[[feature_name_idx]])
    )
    return(gff)
}


prepare_tss <- function(gff, feature_type, feature_name, verbose) {
    if(gff == '-') {
        xgff <- read.table(file("stdin"), sep= '\t', header= FALSE, comment.char= '#', stringsAsFactors = FALSE, quote= '')
        xgff <- as.data.table(xgff)
    } else if(grepl('^/dev/fd/', gff) == TRUE) {
        fin <- fifo(gff, open = "r")
        xgff <- as.data.table(read.table(fin, sep= '\t', header= FALSE, comment.char= '#', stringsAsFactors = FALSE, quote= ''))
        close(fin)
    } else {
        if(isGzip(gff) == TRUE) {
            xcat <- 'gzip -cd'
        } else {
            xcat <- 'cat'
        }
        xgff <- fread(cmd= sprintf('%s %s | grep -v "^#"', xcat, gff), header= FALSE, sep= '\t')
    }

    if(get_data_type(xgff) == 'bed') {
        feature_name <- ifelse(is.na(suppressWarnings(as.integer(feature_name))), 4, as.integer(feature_name))
        xgff <- convert_bed_to_gff(xgff, feature_type, feature_name)
    }

    xgff <- xgff[xgff[[3]] ==  feature_type]
    
    if(nrow(xgff) == 0) {
        stop(sprintf('There are no records of type [%s] in file "%s"', feature_type, gff))
    }
    setnames(xgff, c('V1', 'V4', 'V5', 'V7', 'V9'), c('chrom', 'start', 'end', 'strand', 'gff_attr'))
    stopifnot(xgff$strand %in% c('+', '-'))
    xgff[, start := start - 1]
    xgff[, tss_start := ifelse(strand == '+', start, end - 1)]
    
    xgff[, id := get_gff_attribute(gff_attr, feature_name)]
    xgff[, id := sprintf('%s_%s', 1:nrow(xgff), id)]

    outnames <- c('chrom', 'tss_start', 'id', 'strand', 'start', 'end')
    xgff <- xgff[, outnames, with= FALSE]
    tss <- rbind(
            xgff[strand == '+'][, list(tss_start= min(tss_start)), by= list(chrom, id)],
            xgff[strand == '-'][, list(tss_start= max(tss_start)), by= list(chrom, id)]
            )
    xgff <- merge(xgff, tss, by= c('chrom', 'id', 'tss_start'), sort= FALSE, allow.cartesian=TRUE)
    xgff$tss_end <- xgff$tss_start + 1
    xgff$unused <- '.'
    xgff <- xgff[, c('chrom', 'tss_start', 'tss_end', 'id', 'unused', 'strand', 'start', 'end'), with= FALSE]
    xgff <- unique(xgff[order(chrom, tss_start, tss_end)])
    xgff[, chrom := as.character(chrom)]
    setnames(xgff, names(xgff)[1], paste0('#', names(xgff)[1])) # Comment out header line

    if(verbose) {
        write(sprintf('Found %s features\n', nrow(xgff)), stderr())
    }
    stopifnot(length(xgff$id) == length(unique(xgff$id)))
    return(xgff)
}

closest <- function(peaks, tss, verbose= FALSE, tmpdir, ignore= '', nclosest= 1) {

    peaks <- copy(peaks)
    setcolorder(peaks, c('#chrom', 'summit_start', 'summit_end'))
    setnames(peaks, '#chrom', '#summit_chrom')

    summit_file <- file.path(tmpdir, 'summit.bed') 
    write.table(peaks[order(peaks[[1]], peaks[[2]])], summit_file, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- file.path(tmpdir, 'tss.bed')
    write.table(tss[order(tss[[1]], tss[[2]])], tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <- file.path(tmpdir, 'intx_closest.bed') 
    
    if(ignore == '') {
        
    } else if(ignore == 'downstream') {
        ignore <- '-id'
    } else if(ignore == 'upstream') {
        ignore <- '-iu'
    } else {
        stop(sprintf('Invalid value for ignore: %s', ignore))
    }

    cmd <- sprintf('set -e \
    closestBed -k %s -a %s -b %s -D a %s > %s', nclosest, tss_file, summit_file, ignore, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    if(exit != 0) {
        stop(sprintf('Executing %s', cmd))
    }

    header <- c(names(tss), names(peaks), 'distance')
    
    intx <- fread(out_intx, header= FALSE, sep= '\t')
    if(nrow(peaks) == 0) {
        # With no peaks closestBed returns only chrom, start, end, distance. We
        # need to add the id column
        stopifnot(ncol(intx) == length(header) - 1)
        stopifnot(intx[[ncol(intx) - 3]] == '.')
        stopifnot(intx[[ncol(intx) - 2]] == -1)
        stopifnot(intx[[ncol(intx) - 1]] == -1)
        stopifnot(intx[[ncol(intx) - 0]] == -1)
        intx[, id__ := '.']
        setcolorder(intx, c(1:(ncol(intx)-2), ncol(intx)))
    } 
    setnames(intx, names(intx), header)
    intx[, distance := ifelse(`#summit_chrom` == '.', NA, distance)]
    intx[, `#summit_chrom` := NULL]
    intx[, unused := NULL]
    intx[[1]] <- as.character(intx[[1]])
    return(intx) 
}

intervening_tss <- function(annotated_tss, tss, verbose, tmpdir= tmpdir) {
    # We want to know how many TSSs occur between each TSS and the nearest peak

    # Extend coordinates from TSS to PEAK
    xspan <- annotated_tss[!is.na(distance), list(`#chrom`, 
        start= ifelse(tss_start < summit_start, tss_start, summit_start),
        end= ifelse(tss_end > summit_end, tss_end, summit_end),
        id,
        id__)][order(`#chrom`, start, end)]

    span_file <- file.path(tmpdir, 'span.bed')
    write.table(xspan, span_file, row.names= FALSE, sep= '\t', quote= FALSE)

    tss_file <- file.path(tmpdir, 'tss.bed')
    write.table(tss, tss_file, row.names= FALSE, sep= '\t', quote= FALSE)

    out_intx <-  file.path(tmpdir, 'intervening_tss.bed')

    cmd <- sprintf('set -e \
    intersectBed -a %s -b %s -c > %s', span_file, tss_file, out_intx)

    if(verbose) {
        write(sprintf('%s\n', cmd), stderr())
    }

    exit <- system(cmd)
    if(exit != 0) {
        stop(sprintf('Executing %s', cmd))
    }
    intervene <- suppressWarnings(fread(out_intx, col.names= c(names(xspan), 'n_intra_tss')))
    if(ncol(intervene) == 0) {
        intervene <- data.table(id= character(), id__= character(), n_intra_tss= integer())
    }
    intervene <- intervene[, list(id, id__, n_intra_tss)]
    return(intervene)
}

get_tmpdir <- function() {
    tmpdir <- tempfile(pattern= paste0('tmp_assignPeaksToGenes.', datestr()), tmpdir= '.')
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
        offset <- ceiling((wrk_peaks$peak_end - wrk_peaks$peak_start)/2)
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

keep_me <- function(x) {
    # Return logical vector with TRUE for the elements of x to be kept
    if(length(x) == 1) {
        return(TRUE)
    }
    if(any(!is.na(x))) {
        return(!is.na(x))
    } 
    stop('Unexpected condition')
}

# -----------------------

if(args$n_closest <= 0) {
    stop(sprintf('Argument to --n-closest must be >= 1. Got %s', args$n_closest))
}

# Read peak file and add id__ column
peaks <- peak_reader(args$peaks)
stopifnot(!grepl('id__', names(peaks)))
peaks[, 2] <- as.integer(unlist(peaks[, 2]))
peaks[, 3] <- as.integer(unlist(peaks[, 3]))
peaks[, id__ := sprintf('%s_%s_%s', peaks[[1]], peaks[[2]], peaks[[3]])]
setnames(peaks, 1, '#chrom')

wrk_peaks <- prepare_peaks(peaks, args$summit)

tss <- prepare_tss(args$gff, args$feature_type, args$feature_name, args$verbose)

tryCatch({
        tmpdir <- get_tmpdir()
        closest_downstream <- closest(wrk_peaks, tss, verbose= args$verbose, tmpdir= tmpdir, ignore= 'upstream', nclosest= args$n_closest)
        closest_upstream <- closest(wrk_peaks, tss, verbose= args$verbose, tmpdir= tmpdir, ignore= 'downstream', nclosest= args$n_closest)
        xclosest <- unique(rbind(closest_downstream, closest_upstream))
        intervene <- intervening_tss(xclosest, tss, verbose= args$verbose, tmpdir= tmpdir)
    }, finally= {
        unlink(tmpdir, recursive= TRUE)
    })
xclosest[, `#chrom` := as.character(`#chrom`)]

xclosest <- xclosest[, .SD[keep_me(distance)], by= list(`#chrom`, tss_start, tss_end)][order(`#chrom`, tss_start)]
xclosest[, summit_start := NULL]
xclosest[, summit_end := NULL]
xclosest[, tss_start := NULL]
xclosest[, tss_end := NULL]

xclosest <- merge(xclosest, intervene, all.x= TRUE, by= c('id', 'id__'), sort= FALSE, allow.cartesian=TRUE)
xclosest[, id := sub('^\\d+_', '', id)]
xclosest[, n_intra_tss := ifelse(is.na(n_intra_tss), 0, n_intra_tss)]

name_order <- c('#chrom', 'start', 'end', 'id', 'distance', 'strand', 'n_intra_tss')
stopifnot(name_order %in% names(xclosest))
name_order <- sapply(name_order, function(x) sprintf('%s%s', args$add_column_prefix, x))
setnames(xclosest, names(name_order), name_order)
setcolorder(xclosest, name_order)

if(any(name_order %in% names(peaks))) {
    stop('Duplicate column names in output. Choose a different prefix with --add-column-prefix')
}
xclosest <- merge(xclosest, peaks, by.x= c('id__', sprintf('%s#chrom', args$add_column_prefix)), by.y= c('id__', '#chrom'), all.x= TRUE, sort= FALSE, allow.cartesian=TRUE)
xclosest[, id__ := NULL]
setnames(xclosest, sprintf('%s#chrom', args$add_column_prefix), sprintf('#%schrom', args$add_column_prefix))

stopifnot(grepl('^\\d+_', xclosest$id))

write.table(xclosest, file= stdout(), row.names= FALSE, sep= '\t', quote= FALSE)

quit()

