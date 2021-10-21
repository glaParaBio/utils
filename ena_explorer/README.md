<!-- vim-markdown-toc GFM -->

* [Requirements & Installation](#requirements--installation)
* [Usage](#usage)
    * [The `query` command](#the-query-command)

<!-- vim-markdown-toc -->

Simple command line program to query the [ENA
archive](https://www.ebi.ac.uk/ena/browser/home) and download fastq files.

# Requirements & Installation

`ena` is a single python script depending only on
[pandas](https://pandas.pydata.org/).
[Install](https://pandas.pydata.org/docs/getting_started/install.html) pandas
with your favourite strategy, *e.g.* `mamba install pandas` or `conda install
pandas`. 

Add `ena` to your `PATH` (e.g. `~/bin/`) or execute by giving the full
path to the script (*e.g.* `/path/ena -h`).

# Usage

`ena` has two subcommands, `query` and `download`. As usual, use
`-h/--help` for help.

## The `query` command

Get information for one or more accession IDs in a friendly tabular format.
With no accession, return a table of all the availble columns and their
description.

```
ena query -f markdown

| columnId                   | description                                                                                       |
|:---------------------------|:--------------------------------------------------------------------------------------------------|
| study_accession            | study accession number                                                                            |
| secondary_study_accession  | secondary study accession number                                                                  |
| sample_accession           | sample accession number                                                                           |
| secondary_sample_accession | secondary sample accession number                                                                 |
| experiment_accession       | experiment accession number                                                                       |
| run_accession              | run accession number                                                                              |
| submission_accession       | submission accession number                                                                       |
| tax_id                     | taxonomic ID                                                                                      |
| scientific_name            | scientific name                                                                                   |
| instrument_platform        | instrument platform used in sequencing experiment                                                 |
| instrument_model           | instrument model used in sequencing experiment                                                    |
| library_name               | sequencing library name                                                                           |
| library_layout             | sequencing library layout                                                                         |
| nominal_length             | average fragmentation size of paired reads                                                        |
| library_strategy           | sequencing technique intended for the library                                                     |
| library_source             | source material being sequenced                                                                   |
| library_selection          | method used to select or enrich the material being sequenced                                      |
| read_count                 | number of reads                                                                                   |
| base_count                 | number of base pairs                                                                              |
| center_name                | Submitting center                                                                                 |
...Etc
```

Get the **full** metadata table for project `PRJNA433164`

```
ena query -f markdown PRJNA433164 | less -S
```

Select only some columns:

```
ena query -f markdown -i 'run_accession|sample_title|fastq_ftp|read_count' PRJNA433164 

| run_accession   |   read_count | fastq_ftp                                                                                                                                         | sample_title   |
|:----------------|-------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------|:---------------|
| SRR6676668      |     20485394 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/008/SRR6676668/SRR6676668_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/008/SRR6676668/SRR6676668_2.fastq.gz | 0h_R+_1        |
| SRR6676669      |     19382256 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/009/SRR6676669/SRR6676669_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/009/SRR6676669/SRR6676669_2.fastq.gz | 12h_R-_1       |
| SRR6676670      |     22405889 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/000/SRR6676670/SRR6676670_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/000/SRR6676670/SRR6676670_2.fastq.gz | 24h_R-_1       |
| SRR6676671      |     22694892 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/001/SRR6676671/SRR6676671_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/001/SRR6676671/SRR6676671_2.fastq.gz | 30h_R-_1       |
| SRR6676672      |     26032595 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/002/SRR6676672/SRR6676672_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/002/SRR6676672/SRR6676672_2.fastq.gz | 12h_R+_1       |
| SRR6676673      |     20111872 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/003/SRR6676673/SRR6676673_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/003/SRR6676673/SRR6676673_2.fastq.gz | 24h_R+_1       |
| SRR6676674      |     24608940 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/004/SRR6676674/SRR6676674_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/004/SRR6676674/SRR6676674_2.fastq.gz | 30h_R+_1       |
| SRR6676675      |     23258076 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/005/SRR6676675/SRR6676675_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/005/SRR6676675/SRR6676675_2.fastq.gz | 0h_R-_1        |
| SRR6676676      |     20558602 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/006/SRR6676676/SRR6676676_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/006/SRR6676676/SRR6676676_2.fastq.gz | 24h_R-_2       |
| SRR6676677      |     23057299 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/007/SRR6676677/SRR6676677_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/007/SRR6676677/SRR6676677_2.fastq.gz | 0h_R+_2        |
| SRR6676678      |      4493715 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/008/SRR6676678/SRR6676678_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/008/SRR6676678/SRR6676678_2.fastq.gz | 12h_R-_2       |
...Etc
```
