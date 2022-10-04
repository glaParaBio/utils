library(testthat)

context("Test functions")

source('../makeBioconductorAnnotationDbi.r')

test_that("Can convert GFF to dbi table", {
    STANDARD_COLUMNS <- c('GID', 'FEATURE_TYPE', 'CHROM', 'GID_START', 'GID_END', 'GID_STRAND')

    dbi <- gff_to_dbitable('data/PlasmoDB-59_PbergheiANKA.gff.gz', 'AUTO', '')
    expect_equal(names(dbi), STANDARD_COLUMNS)
    expect_equal(dbi$GID_START[1], 1175534)
    expect_equal(dbi$GID_START[nrow(dbi)], 276592)

    # From URL    
    dbi <- gff_to_dbitable('https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gff/data/PlasmoDB-55_PbergheiANKA.gff', 'AUTO', '')
    expect_equal(names(dbi), STANDARD_COLUMNS)
    expect_equal(dbi$GID_START[1], 1175534)
    expect_equal(dbi$GID_START[nrow(dbi)], 276592)

    dbi <- gff_to_dbitable('data/PlasmoDB-59_PbergheiANKA.gff.gz', 'AUTO', 'ALL')
    expect_true(all(STANDARD_COLUMNS %in% names(dbi)))
    expect_true(all(c('Name', 'description', 'ebi_biotype') %in% names(dbi)))
    expect_equal(dbi$GID_START[1], 1175534)
    expect_equal(dbi$GID_START[nrow(dbi)], 276592)

    dbi <- gff_to_dbitable('data/PlasmoDB-59_PbergheiANKA.gff.gz', 'AUTO', c('Name', 'description'))
    expect_equal(sort(names(dbi)), sort(c(STANDARD_COLUMNS, 'Name', 'description')))
    expect_equal(dbi$GID_START[1], 1175534)
    expect_equal(dbi$GID_START[nrow(dbi)], 276592)

    dbi <- gff_to_dbitable('data/PlasmoDB-59_PbergheiANKA.gff.gz', 'AUTO', 'Name')
    expect_equal(sort(names(dbi)), sort(c(STANDARD_COLUMNS, 'Name')))
    expect_equal(dbi$GID_START[1], 1175534)
    expect_equal(dbi$GID_START[nrow(dbi)], 276592)

    expect_true(sum(is.na(dbi)) == 0)
    expect_true('n/a' %in% dbi$Name)

    dbi <- gff_to_dbitable('data/PlasmoDB-59_PbergheiANKA.gff.gz', 'protein_coding_gene', '')
    expect_equal(nrow(dbi), 4945)
    dbi <- gff_to_dbitable('data/PlasmoDB-59_PbergheiANKA.gff.gz', c('protein_coding_gene', 'pseudogene'), '')
    expect_equal(nrow(dbi), 5077)

})

test_that("Can autodetect genes", {
    dbi <- gff_to_dbitable('data/short.gff', 'AUTO', '', include_ids=NULL)
    expect_equal(dbi$FEATURE_TYPE, rep('protein_coding_gene', 6))
    expect_true('PBANKA_1307600' %in% dbi$GID)
    
    dbi <- gff_to_dbitable('data/short.gff', 'AUTO', '', include_ids='exon_PBANKA_1112300.1-E1')
    expect_true('PBANKA_1307600' %in% dbi$GID)
    expect_true('exon_PBANKA_1112300.1-E1' %in% dbi$GID)

    dbi <- gff_to_dbitable('data/short.gff', 'AUTO', '', include_ids=c('exon_PBANKA_1112300.1-E1', 'exon_PBANKA_0832100.1-E1', 'FOOBAR'))
    expect_true('PBANKA_1307600' %in% dbi$GID)
    expect_true('exon_PBANKA_1112300.1-E1' %in% dbi$GID)
    expect_true('exon_PBANKA_0832100.1-E1' %in% dbi$GID)
})
