<!-- vim-markdown-toc GFM -->

* [Quick start for the inpatient](#quick-start-for-the-inpatient)
* [Description](#description)
* [Installation & Usage](#installation--usage)

<!-- vim-markdown-toc -->

Quick start for the inpatient 
===========

The minimal input to create an annotationDbi package is a gff file and tab
separated file assigning GO terms to gene IDs:

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
or GAF should also work.

If using VEuPathDB (and if things have not changed), you can find the data files in:

```
https://{database}.org/common/downloads/release-{N}/{species}/{gff|gaf}/data/
```

E.g.:

```
https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gff/data/PlasmoDB-55_PbergheiANKA.gff 
https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gaf/PlasmoDB-55_PbergheiANKA_GO.gaf
```

Installation & Usage
=====

It depends on R and the packages in [requirements.txt](requirements.txt). You
can use mamba/conda to install them. *E.g.*

```
mamba install --file requirements.txt
```

any reasonably recent version of R and these packages should do.

```
./makeBioconductorAnnotationDbi.r -h
usage: ./makeBioconductorAnnotationDbi.r [-h] --gff GFF --gaf GAF
                                  [--pckg-version PCKG_VERSION] --genus GENUS
                                  --species SPECIES [--taxid TAXID]
                                  [--outdir OUTDIR] [--maintainer MAINTAINER]
                                  [--author AUTHOR] [--install] [--version]

Prepare a Bioconductor annotationDbi (i.e. a package like org.Hs.eg.db) for
organisms in VEuPathDB or for any data files formatted accordingly

optional arguments:
  -h, --help            show this help message and exit
  --gff GFF             GFF local file or URL
  --gaf GAF             GAF local file or URL to extract GO annotation
  --pckg-version PCKG_VERSION, -p PCKG_VERSION
                        Version string. Must be in decimal format e.g.
                        "0.1.2". Consider using the VEuPathDB version of the
                        source data [0.1]
  --genus GENUS, -g GENUS
                        Genus name [None]
  --species SPECIES, -s SPECIES
                        Species name [None]
  --taxid TAXID, -t TAXID
                        Taxonomy ID of this species. See
                        https://www.ncbi.nlm.nih.gov/taxonomy [0]
  --outdir OUTDIR, -o OUTDIR
                        Output directory [.]
  --maintainer MAINTAINER, -m MAINTAINER
                        Maintainer in format "Some One <so@someplace.org>" [NA
                        <na@na.com>]
  --author AUTHOR, -a AUTHOR
                        Name and address of the author[n/a]
  --install, -I         If set, also install the package
  --version, -v         show program's version number and exit
```

(Printout above maybe outdated)

See [test.py](test.py) and the `test` directory for examples that should
actually complete without errors.
