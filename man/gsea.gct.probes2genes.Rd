\name{gsea.gct.probes2genes}
\alias{gsea.gct.probes2genes}
\title{Take a table of expression data (at the probe-level), and convert it to a
table of expression data at the unique gene-symbol-level. This chooses the
best probe for each gene by looking at the stats vector which should be
equivalent to the contents of a rnk file (todo: allow a rnk file to be
specified) There's 2 ways of chooseing the best probe - Mark's way
(Best_probe), or the default GSEA way which is called "Max_probe".}
\usage{
  gsea.gct.probes2genes(data, stats, chip, mode =
  c("Best_probe", "Max_probe", "Median_of_probes"))
}
\arguments{
  \item{data}{a data.frame of expression data}

  \item{stats}{a vector of statistics}

  \item{chip}{either the filename of a chip file, or the
  result of import.gsea.chip}

  \item{mode}{Best_probe or Max_probe}
}
\value{
  a data.frame of expression data with 1 row per unique
  gene symbol, and gene symbols as the rownames see also
  gsea.genesets2gct.
}
\description{
  Take a table of expression data (at the probe-level), and
  convert it to a table of expression data at the unique
  gene-symbol-level. This chooses the best probe for each
  gene by looking at the stats vector which should be
  equivalent to the contents of a rnk file (todo: allow a
  rnk file to be specified) There's 2 ways of chooseing the
  best probe - Mark's way (Best_probe), or the default GSEA
  way which is called "Max_probe".
}
\author{
  Mark Cowley, 2009-07-27
}

