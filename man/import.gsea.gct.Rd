\name{import.gsea.gct}
\alias{import.gsea.gct}
\title{Import a gct file...}
\usage{import.gsea.gct(f, clm.file, check.names=FALSE)
}
\description{Import a gct file}
\details{Import a GenePattern gct file into R as a data.frame object. 
If a clm.file is also provided, then the columns of the gct are optionally subsetted, reordered, and renamed as per the clm.}
\seealso{\url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#gct}, \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#clm}}
\value{a \code{data.frame} of usually gene expression data. columns 1&2 are annotation, 3+ are data.}
\author{Mark Cowley, 2009-07-27}
\arguments{\item{f}{[character] the path to a gct file}
\item{clm.file}{[optional] path to a clm.file}
\item{check.names}{logical: R likes to make sure that the column names are 'syntactically valid', If you don't care, then choose FALSE. default is FALSE.}
}

