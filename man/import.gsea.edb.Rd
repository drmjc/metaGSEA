\name{import.gsea.edb}
\alias{import.gsea.edb}
\title{Import a GSEA edb file.}
\usage{import.gsea.edb(x)
}
\description{Import a GSEA edb file.
The edb file is a machine readable format of a GSEA run, which the GSEA
desktop tool can open. 
This function is only really indented to be run by \code{\link{import.gsea}} with 
\code{edb=TRUE}, and it's a good idea to import this file if you want to
start manipulating, subsetting, or mergeing GSEA runs and still wish to use
the GSEA tool (such as the GSEA Leading Edge tool) in the future.}
\details{These edb files are quite large, eg 22 MB for c2_all with 1000 permutations. Size
is proportional to the number of permutations.

The edb file is quite simple XML, utilising one node per geneset, and lots
of attributes but not sub-terms for each geneset. We use the \code{XML} package,
originally from Duncan Temple Lang.}
\value{The edb file in the form of an XMLDocument. See \code{\link[XML]{xmlTreeParse}},
This object appears to be immutable.}
\author{Mark Cowley, 2009-10-12}
\seealso{\code{\link[XML]{xmlTreeParse}}, \code{\link{import.gsea}},
\code{\link{gsea.rename.genesets.edb}}}
\references{Omega hat XML: \url{http://www.omegahat.org/RSXML/} 
Quick Guide to the XML package: \url{http://www.omegahat.org/RSXML/Tour.pdf}}
\keyword{file}
\keyword{IO}
\arguments{\item{x}{A GSEA directory.}
}
\examples{# not run
# edbT <- import.gsea.edb("T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710")
# edbM <- import.gsea.edb("MCF7_TreatedVsUntreated/c2_all.GseaPreranked.1236217933720")
# edb.merged <- merge.gsea.edb("T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710", "MCF7_TreatedVsUntreated/c2_all.GseaPreranked.1236217933720", c("T47D_", "MCF7_"), c("", ""))
# export.gsea.edb(edb.merged, "./merged")}

