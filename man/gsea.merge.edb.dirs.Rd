\name{gsea.merge.edb.dirs}
\alias{gsea.merge.edb.dirs}
\title{Function to merge the edb files from 2 GSEA runs.}
\usage{gsea.merge.edb.dirs(gsea.dir1, gsea.dir2, prefix, suffix)
}
\description{Function to merge the edb files from 2 GSEA runs.}
\details{Why?
so you can examine them both together in the Leading Edge analysis tool of
course!
Warnings:
The resultant edb must have the same RANKED_LIST attribute for all nodes.
Thus the RANKED_LIST
attribute for the 2nd GSEA run is changed to the name of the first rnk
which is obviously a bit
of a hack.}
\value{none. merges 2 directories}
\author{Mark Cowley, 2009-10-06}
\arguments{\item{gsea.dir1}{paths to two GSEA directories (the ones that contain
index.html)}
\item{gsea.dir2}{paths to two GSEA directories (the ones that contain
index.html)}
\item{prefix}{both should be character(2) indicating how you want the
genesets to be renamed so that they reflect which experimental run they
came from.}
\item{suffix}{both should be character(2) indicating how you want the
genesets to be renamed so that they reflect which experimental run they
came from.}
}

