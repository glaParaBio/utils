#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GO.db))

VERSION = '0.2.0'

get_command_call <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    prog <- grepl('^--file=', cmdArgs)
    prog <- cmdArgs[prog]
    stopifnot(length(prog) == 1)
    prog <- sub('^--file=', '', prog)
    cmdCall <- c(prog, commandArgs(trailingOnly = TRUE))
    cmdCall <- paste(cmdCall, collapse=' ')
    return(cmdCall)
}

reader <- function(fileOrUrl, ...) {
    options(warn=2)
    grepcmd <- "awk '{if($0 ~ /^##FASTA/) {exit 0} else {print $0}}' | grep -v -P '^#|^!gaf-version'"
    if(grepl('^http', fileOrUrl) == TRUE) {
        cmd <- sprintf('curl -s -L %s | %s', fileOrUrl, grepcmd)
    } else if(grepl('*.gz$', fileOrUrl) == TRUE) {
        cmd <- sprintf('gzip -cd %s | %s', fileOrUrl, grepcmd) 
    } else {
        cmd <- sprintf("cat %s | %s", fileOrUrl,  grepcmd)
    }
    ff <- fread(cmd=cmd, ...)
    options(warn=0)
    return(ff)
}

get_package_name <- function(outdir, genus, species) {
    g <- toupper(substr(genus, 1, 1))
    pkg <- sprintf('org.%s%s.eg.db', g, species)
    path <- file.path(outdir, pkg)
    return(path)
} 

delete_existing_package <- function(path) {
    if(file.exists(path) == TRUE) {
        sqlite <- file.path(path, 'inst', 'extdata', '*.sqlite')
        options(warn=2)
        xcode <- system(sprintf('rm -f %s', sqlite))
        xcode <- system(sprintf('rm -r %s', path))
        options(warn=0)
    }
} 

autodetect_genes <- function(gff, include) {
    # Filter gff table to extract likely genes. We use column FEATURE_TYPE (3rd
    # column) to guess genes and we augment it with any gene in the `include` vector. 
    # Typically, `include` is the list of genes from the GAF file
    stopifnot(c('ID', 'FEATURE_TYPE') %in% names(gff))
    KEEP_FEATURES <- c('gene', 'pseudogene', 'ncRNA_gene', 'protein_coding_gene') 
    features <- which(gff$FEATURE_TYPE %in% KEEP_FEATURES)
    extra <- which(gff$ID %in% include)
    keep <- sort(unique(c(features, extra)))
    return(gff[keep])
}

gff_to_dbitable <- function(gff_file, feature_types='AUTO', keep_attributes='ALL', include_ids=NULL) {
    # Read gff file and return table suitable for annotationForge

    gff <- reader(gff_file, header=FALSE, sep='\t')
    stopifnot(names(gff) == sprintf('V%s', 1:9))
    setnames(gff, c('V1', 'V3', 'V4', 'V5', 'V7', 'V9'), c('CHROM', 'FEATURE_TYPE', 'GID_START', 'GID_END', 'GID_STRAND', 'ATTRIBUTES'))

    gff[, ID := get_gff_attribute(ATTRIBUTES, 'ID')]

    if(length(feature_types) == 1 && feature_types == 'AUTO') {
        gff <- autodetect_genes(gff, include_ids)
    } else {
        gff <- gff[FEATURE_TYPE %in% feature_types]
    }
    if(nrow(gff) == 0) {
        stop('No usable feature found in GFF file')
    }

    gff <- gff[, list(ID, FEATURE_TYPE, CHROM, GID_START, GID_END, GID_STRAND, ATTRIBUTES)]

    attributes <- attributes_to_table(gff$ATTRIBUTES, keep=keep_attributes)
    if(!is.null(attributes)) {
        stopifnot(attributes$ID == gff$ID)
        gff <- merge(gff, attributes, by='ID', sort=FALSE)
    }
    setnames(gff, 'ID', 'GID')
    gff[, ATTRIBUTES := NULL]
    gff[is.na(gff)] <- 'n/a'
    return(gff)
}

attributes_to_table <- function(gff_attr, keep='ALL') {
    # Convert the vector of GFF attributes (column 9 in GFF) to data.table
    if(is.na(keep) || is.null(keep) || length(keep) == 0 || keep == '') {
        keep <- 'ID'
    }
    
    keys <- unlist(strsplit(gff_attr, ';'))
    keys <- unique(sapply(keys, function(x) strsplit(x, '=')[[1]][1], USE.NAMES=FALSE))
    if(length(keep) == 1 && keep == 'ALL') {
        keep <- keys
    }

    if(!'ID' %in% keep) {
        keep <- c('ID', keep)
    }

    keys <- keys[keys %in% keep]
    attr_dt <- list()
    for(x in keys) {
        raw <- get_gff_attribute(gff_attr, x, na_value='n/a')
        decoded <- sapply(raw, URLdecode, USE.NAMES=FALSE)
        attr_dt[[x]] <- decoded
    }
    attr_dt <- as.data.table(attr_dt)
    return(attr_dt)
}

get_gff_attribute <- function(gff_attr, key, na_value=NA) {
    values <- rep(NA, length(gff_attr))
    for(i in 1:length(gff_attr)) {
        attr_list <- strsplit(gff_attr[i], ';')[[1]]
        attr_list <- sapply(attr_list, sub, pattern= '^ ', replacement= '', USE.NAMES= FALSE)
        keep <- grepl(sprintf('^%s( |=)', key), attr_list)
        if(sum(keep) > 1) {
            stop(sprintf('Attribute "%s" found more than once in\n%s', key, gff_attr[i]))
        }
        if(sum(keep) == 0) {
            value <- na_value
        } else {
            value <- attr_list[keep]
            value <- sub(sprintf('^%s( |=)', key), '', value)
            value <- gsub('^"|"$', '', value)
        }
        values[i] <- value
    }
    return(values)
}

edit_package_description <- function(description_file, cmdCall, version) {
    content <- readLines(description_file)
    descr_idx <- grep('^Description:', content)
    description <- content[descr_idx]
    description <- sub(', primarily based on mapping using Entrez Gene identifiers', '', description, fixed=TRUE)
    description <- paste(description, 'Created by: [', cmdCall, ']', ', version:', version, collapse=' ')
    content[descr_idx] <- description
    writeLines(content, description_file)
}

is_validate_name_for_package <- function(name) {
    # From https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Creating-R-packages:
    #
    # The mandatory ‘Package’ field gives the name of the package. This should
    # contain only (ASCII) letters, numbers and dot, have at least two
    # characters and start with a letter and not end in a dot. 
    # 
    # We need to check that species and genus contain only numbers, letters, and
    # dots. Check on length, start, and end of name are not necessary since the
    # package will always start with 'org.' and end in '.db'
    is_valid <- !grepl("[^A-Za-z0-9.]", name)
    return(is_valid)
}

parser <- ArgumentParser(description='Prepare a Bioconductor annotationDbi (i.e. a package like org.Hs.eg.db) given a GFF file and a GAF file of gene IDs and associated GO terms')
parser$add_argument('--gff', help='GFF local file or URL [required]', required=TRUE)
parser$add_argument('--gaf', help='GAF local file or URL linking GO terms to genes [required]', required=TRUE)
parser$add_argument('--genus', '-g', help='Genus name [required]', required=TRUE)
parser$add_argument('--species', '-s', help='Species name [required]', required=TRUE)
parser$add_argument('--pckg-version', '-p', help='Version you want to give to this package. Must be in decimal format e.g. "0.1.2" [%(default)s]', default='0.1')
parser$add_argument('--taxid', '-t', help='Taxonomy ID of this species. Find yours at https://www.ncbi.nlm.nih.gov/taxonomy [%(default)s]', default=0, type='integer')
parser$add_argument('--outdir', '-o', help='Output directory [%(default)s]', default='.')
parser$add_argument('--maintainer', '-m', help='Maintainer name and email in format "Some One <so@someplace.org>" [%(default)s]', default='NA <na@na.com>')
parser$add_argument('--author', '-a', help='Name and address of the author [%(default)s]', default='n/a')
parser$add_argument('--feature-types', '-f', help='Use GFF records where the feature type (i.e. column 3) is in this list. Use AUTO to autodetect [%(default)s]', nargs='+', default='AUTO')
parser$add_argument('--attributes', '-A', help='Include these GFF attributes. Use ALL to include all of them [%(default)s]', nargs='*', default='ALL')
def <- c(2, 5, 7)
parser$add_argument('--gaf-column-idx', '-x', help=sprintf('Index of columns in GAF input for, respectively: GID, GO, EVIDENCE [%s]', paste(def, collapse=' ')), nargs=3, default=def, type='integer')
parser$add_argument('--install', '-I', action='store_true', help='If set, also install the package in the default R library')
parser$add_argument('--version', '-v', action= 'version', version=VERSION)

if(sys.nframe() == 0){
    # Script is being executed from the command line

    xargs <- parser$parse_args()

    for(x in c(xargs$genus, xargs$species)) {
        if(is_validate_name_for_package(x) == FALSE) {
            write(sprintf('\nInvalid name: "%s"\nGenus and species name can contain only letters, numbers, and dots\n', x), stderr())
            quit(status=1)
        }
    }

    idx <- xargs$gaf_column_idx
    names(idx) <- c('GID', 'GO', 'EVIDENCE')
    gaf <- reader(xargs$gaf, header=FALSE, sep='\t', select=as.numeric(idx), col.names=names(idx))
    gaf <- unique(gaf)

    gff <- gff_to_dbitable(xargs$gff, feature_types=xargs$feature_types,
                           keep_attributes=xargs$attributes,
                           include_ids=gaf$GID)

    pkg_name <- get_package_name(xargs$outdir, xargs$genus, xargs$species)
    delete_existing_package(pkg_name)

    dir.create(xargs$outdir, showWarnings=FALSE, recursive=TRUE)
    
    sink(stderr(), type = "output")
    suppressMessages(
    name <- makeOrgPackage(gene_info=gff, go=gaf,
                   version=xargs$pckg_version,
                   maintainer=xargs$maintainer,
                   author=xargs$author,
                   outputDir=xargs$outdir,
                   tax_id=xargs$taxid,
                   genus=xargs$genus,
                   species=xargs$species,
                   verbose=FALSE,
                   goTable="go")
    )
    sink(NULL, type="output")
    stopifnot(path.expand(pkg_name) == path.expand(name))
    chmod <- sprintf('chmod --reference=%s/DESCRIPTION %s/inst/extdata/*.sqlite', pkg_name, pkg_name)
    system(chmod)

    cmd <- get_command_call()
    edit_package_description(file.path(pkg_name, 'DESCRIPTION'), cmd, VERSION)

    if(xargs$install == TRUE) {
        install.packages(pkg_name, repos=NULL, type='source', quiet=TRUE)
    } else {
        cmd <- sprintf("R CMD INSTALL %s", pkg_name)
        write(sprintf('To install execute:\n%s', cmd), stderr())
    }
}
