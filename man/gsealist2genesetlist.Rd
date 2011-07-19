\name{gsealist2genesetlist}
\alias{gsealist2genesetlist}
\title{Extract the geneset names from a list of GSEA objects.}
\usage{gsealist2genesetlist(x)
}
\description{Extract the geneset names from a list of GSEA objects.}
\details{From a list of GSEA objects (like what you get from gsea.split), extract
the geneset names, and create a gmx-like object. i.e. a table where each
column contains a vector of geneset names mentioned in each of the GSEA
objects in the list.}
\value{a data.frame with colnames == names(x), no rownames, where each
column contains the genesets which were in the relevant object inside x.}
\author{Mark Cowley, 2009-10-29}
\arguments{\item{x}{a list of GSEA objects}
}

