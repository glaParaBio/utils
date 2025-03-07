<!-- vim-markdown-toc GFM -->

* [Quick start for the inpatient](#quick-start-for-the-inpatient)
* [Description](#description)
* [Installation](#installation)
* [Usage](#usage)
    * [Input](#input)
* [Notes](#notes)

<!-- vim-markdown-toc -->

You need an [annotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) package to use some Bioconductor package (e.g.
[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html))
but your species doesn't have one. Here's a script to make it.

Quick start for the inpatient 
===========

The minimal input to create an
[annotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
package is a gff file and a tab separated file assigning GO terms to gene IDs:

```
./makeBioconductorAnnotationDbi.r \
    --gff https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gff/data/PlasmoDB-55_PbergheiANKA.gff \
    --gaf https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gaf/PlasmoDB-55_PbergheiANKA_GO.gaf \
    -g Plasmodium \
    -s bergheiANKA \
    --install
```

the above creates package `org.PbergheiANKA.eg.db` in the local directory and it installs it
in the default R library.

Description
===========

Create an
[AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
package (like e.g.
`org.Hs.eg.db`)
suitable for the many Bioconductor packages that require an AnnotationDbi to work.
`makeBioconductorAnnotationDbi.r` is designed around GFF and
[GAF](http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/) files
from [VEuPathDB](https://veupathdb.org/veupathdb/app/) although any GFF
and a tab-delimited file should also work.

If you work with species in VEuPathDB (and if things have not changed), you can
find the data files in:

```
https://{database}.org/common/downloads/release-{N}/{species}/{gff|gaf}/data/
```

E.g.:

```
https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gff/data/PlasmoDB-55_PbergheiANKA.gff 
https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gaf/PlasmoDB-55_PbergheiANKA_GO.gaf
```

Installation
============

It depends on R and the packages in [requirements.txt](requirements.txt). You
can use mamba/conda to install them. *E.g.*

```
mamba install --file requirements.txt
```

any reasonably recent version of R and these packages should do. However, it
appears that package AnnotationForge version < 1.36 throws some (harmless)
warnings about SQL statements.

Usage
=====

```
usage: ./makeBioconductorAnnotationDbi.r [-h] --gff GFF --gaf GAF --genus
                                         GENUS --species SPECIES
                                         [--pckg-version PCKG_VERSION]
                                         [--taxid TAXID] [--outdir OUTDIR]
                                         [--maintainer MAINTAINER]
                                         [--author AUTHOR]
                                         [--feature-types FEATURE_TYPES [FEATURE_TYPES ...]]
                                         [--attributes [ATTRIBUTES ...]]
                                         [--gaf-column-idx GAF_COLUMN_IDX GAF_COLUMN_IDX GAF_COLUMN_IDX]
                                         [--install] [--version]

Prepare a Bioconductor annotationDbi (i.e. a package like org.Hs.eg.db) given
a GFF file and a GAF file of gene IDs and associated GO terms

options:
  -h, --help            show this help message and exit
  --gff GFF             GFF local file or URL. Gzip input ok provided the
                        filename ends with .gz [required]
  --gaf GAF             GAF local file or URL linking GO terms to genes. Gzip
                        input ok provided the filename ends with .gz
                        [required]
  --genus GENUS, -g GENUS
                        Genus name [required]
  --species SPECIES, -s SPECIES
                        Species name [required]
  --pckg-version PCKG_VERSION, -p PCKG_VERSION
                        Version you want to give to this package. Must be in
                        decimal format e.g. "0.1.2" [0.1]
  --taxid TAXID, -t TAXID
                        Taxonomy ID of this species. Find yours at
                        https://www.ncbi.nlm.nih.gov/taxonomy [0]
  --outdir OUTDIR, -o OUTDIR
                        Output directory [.]
  --maintainer MAINTAINER, -m MAINTAINER
                        Maintainer name and email in format "Some One
                        <so@someplace.org>" [NA <na@na.com>]
  --author AUTHOR, -a AUTHOR
                        Name and address of the author [n/a]
  --feature-types FEATURE_TYPES [FEATURE_TYPES ...], -f FEATURE_TYPES [FEATURE_TYPES ...]
                        Use GFF records where the feature type (i.e. column 3)
                        is in this list. Use AUTO to autodetect [AUTO]
  --attributes [ATTRIBUTES ...], -A [ATTRIBUTES ...]
                        Include these GFF attributes. Use ALL to include all
                        of them [ALL]
  --gaf-column-idx GAF_COLUMN_IDX GAF_COLUMN_IDX GAF_COLUMN_IDX, -x GAF_COLUMN_IDX GAF_COLUMN_IDX GAF_COLUMN_IDX
                        Index of columns in GAF input for, respectively: GID,
                        GO, EVIDENCE [2 5 7]
  --install, -I         If set, also install the package in the default R
                        library
  --version, -v         show program's version number and exit
```

(Printout above maybe outdated)


See [test.py](test.py) and the `test` directory for examples that should
actually complete without errors.

Input
-----

* **GFF** files reasonably well behaved should work (GTF input is not supported, if
  GTF is what you have google for GTF to GFF converters). As per GFF format,
  genes should have attribute `ID` which identifies the feature ID. By default,
  gene feature are autodetected see function `autodetect_genes` for details.
  Alternatively, use the `--feature-types` argument to specify which features
  should be used. 

  If all you have is a list of genes, you can fake a GFF file by making a
  9-column, tab-separated file with "gene" in column 3 and `ID=<gene_id>` in
  column 9. Add any other gene information, like "description" as an attribute. E.g.

```
NA | NA | gene | NA | NA | NA | NA | NA | ID=PBANKA_0932200;description=A Fake gene
NA | NA | gene | NA | NA | NA | NA | NA | ID=PBANKA_1112300
NA | NA | gene | NA | NA | NA | NA | NA | ID=PBANKA_0832100
NA | NA | gene | NA | NA | NA | NA | NA | ID=PBANKA_1307600
NA | NA | gene | NA | NA | NA | NA | NA | ID=PBANKA_0818100
NA | NA | gene | NA | NA | NA | NA | NA | ID=PBANKA_1329900
```

* **Gene annotation file**: A tab-separated file with no header with three
  columns (additional columns ignored) for gene id, go, and evidence.

Input files can be gzip compressed, but they must have extension .gz for the
script to recognise them as such. Lines starting with `#` or `!` are skipped.

Notes
=====

* To run tests:

```
./test/test.py
./test/testthat.r
```

* Functions to parse attributes in GFF file is pretty slow - consider refactoring it
