\name{gsea.combine}
\alias{gsea.combine}
\title{Combine N gsea runs together...}
\usage{gsea.combine(..., gsea.runs, verbose=TRUE)
}
\description{Combine N gsea runs together}
\details{This merges N independent GSEA results into a single much larger GSEA result
object. So that you don't end up with N copies of each geneset, you need to
rename the genesets before running this code to uniquify them, eg adding
a short prefix tag (eg \dQuote{MCF7} or \dQuote{T47D}) is enough to uniquify the
geneset names.
You should probably also run some sort of geneset filtering, to keep the result
size manageble, and since filtering on this combined 
list may fail.}
\value{a GSEA object with larger tt, and larger leading.edge objects.
Currently, the rnk and edb objects are set to NA -- WARNING: this may
make routines that work on GSEA objects fail. We recommend treating
results from gsea.combine as immutable objects.}
\author{Mark Cowley, 2010-08-16}
\seealso{\code{\link{gsea.rename.genesets}} \code{\link{gsea.filter}}}
\arguments{\item{\dots}{at least 2 individual GSEA objects, or}
\item{gsea.runs}{a list of GSEA objects. This overrides \dots}
\item{verbose}{logical.}
}

