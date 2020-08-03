metaGSEA
========

An R package for I/O of GenePattern datatypes, and running meta-GSEA analyses

Description
===========
A collection of R and bash code to simplify
the down-stream analysis of GSEA results. Efficient
methods for importing and storing GSEA outputs,
comparing genesets within a single GSEA run, and
comparing genesets between different GSEA runs. Useful
visualisation is the key strength of this package,
utilising GSEA bar plots, correspondence at the top
(CAT) plots, and hierarchical clustering plots of
geneset similarities. Also the ability to clean up GSEA
output on unix-alike systems, by improving the GSEA
pre-ranked output, and the gene tables in each GSEA
report by including the gene symbol, description and
hyperlinks to Entrez Gene. This has evolved into a suite
of code for importing/exporting/subsetting many of the
GenePattern file types (GCT, CLS, CLM, CHIP, RES, ODF)

See also
========
This code empowers the metaGSEA GenePattern module, available at 
http://pwbc.garvan.unsw.edu.au/gp

Future Plans
============
- split the package into genepatternIO and metaGSEA
- add S3 or S4 datatypes for at least GCT, CHIP, CLS objects.

Installation
============

    install.packages(c("XML", "gplots", "stringr", "devtools"))
    source("http://bioconductor.org/biocLite.R")
    biocLite(c("Biobase", "GSEABase"))
    library(devtools)
    install_github("excelIO", "drmjc")
    install_github("mjcbase", "drmjc")
    install_github("mjcgraphics", "drmjc")
    install_github("microarrays", "drmjc")
    install_github("mjcstats", "drmjc")
    install_github("mjcstats", "genomics")

Usage
=====
Extensive package documentations is available via:

	library(metaGSEA)
	?metaGSEA
