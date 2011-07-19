\name{gsealist2adjmat}
\alias{gsealist2adjmat}
\title{Convert a GSEA list into an adjacency matrix at the gene-level.}
\usage{gsealist2adjmat(x)
}
\description{Convert a GSEA list into an adjacency matrix at the gene-level.}
\details{From a list of GSEA objects (like what you get from gsea.split), extract
the genes from the leading edges of each GSEA object, and then convert into
an adjacency matrix.}
\value{a data.frame of 0 or 1 depending on whether each gene was found in
each element of the GSEA list.  See also: list2adjmat}
\author{Mark Cowley, 2009-10-29}
\arguments{\item{x}{a list of GSEA objects.}
}

