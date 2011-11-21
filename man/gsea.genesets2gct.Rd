\name{gsea.genesets2gct}
\alias{gsea.genesets2gct}
\title{Create gct files for all gseaPreranked runs comparing 1 rnk to lots of
gmt's.}
\usage{
  gsea.genesets2gct(gsea.dir, gct.file,
    leading.edge = TRUE, which.sets = c("best", "all"))
}
\arguments{
  \item{gsea.dir}{either a single gsea directory, or a
  parent-dir that contains lots of gsea dirs, (each vs a
  different collection of genesets)}

  \item{gct.file}{a unique gene-level gct file (very
  important to see gsea.gct.probes2genes)}

  \item{leading.edge}{if TRUE then only the leading edge
  genes are exported, or else, the genes in the same order
  as given by edb/gene_sets.gmt are given.}

  \item{which.sets}{do you want to make a gct for every
  single leading edge ("all"), or just those that have been
  created in the dir ("best") which are usually the top 20
  or 50.}
}
\description{
  By default, the leading edge genes are identified by the
  'leading.edge' slot from running import.gsea(); this then
  creates an edb/leading_edge.gmt file. From a gct.file at
  the unique genesymbol level (see gsea.gct.probes2genes),
  the relevant rows are retrieved, keeping the same order
  as the genes in the leading edge itself.
}
\author{
  Mark Cowley, 2009-07-27
}

