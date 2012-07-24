\name{import.gsea.edb}
\alias{import.gsea.edb}
\title{Import a GSEA edb file.
The edb file is a machine readable format of a GSEA run, which the GSEA
desktop tool can open.
This function is only really indented to be run by \code{\link{import.gsea}} with
\code{edb=TRUE}, and it's a good idea to import this file if you want to
start manipulating, subsetting, or mergeing GSEA runs and still wish to use
the GSEA tool (such as the GSEA Leading Edge tool) in the future.}
\usage{
  import.gsea.edb(x)
}
\arguments{
  \item{x}{A GSEA directory.}
}
\value{
  The edb file in the form of an XMLDocument. See
  \code{\link[XML]{xmlTreeParse}}, This object appears to
  be immutable.
}
\description{
  These edb files are quite large, eg 22 MB for c2_all with
  1000 permutations. Size is proportional to the number of
  permutations.
}
\details{
  The edb file is quite simple XML, utilising one node per
  geneset, and lots of attributes but not sub-terms for
  each geneset. We use the \code{XML} package, originally
  from Duncan Temple Lang.
}
\note{
  Here's some code I found useful in deciphering the XML
  document: \code{ edb.file <-
  "T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710/edb/results.edb"
  doc <- xmlTreeParse(edb.file) r <- xmlRoot(doc);
  xmlAttrs(r); xmlSize(r); xmlSize(r[[1]])
  str(xmlAttrs(r[[1]])) xmlAttrs(r[[1]])[3] # modify the
  attributes... a <- xmlAttrs(r[[1]]) a[3] <- paste(a[3],
  "_MJC", sep="") n <- xmlNode("DTG", attrs=a) r[[1]] <- n
  }
}
\examples{
\dontrun{
edbT <- import.gsea.edb("T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710")
edbM <- import.gsea.edb("MCF7_TreatedVsUntreated/c2_all.GseaPreranked.1236217933720")
edb.merged <- merge.gsea.edb("T47D_TreatedVsUntreated/c2_all.GseaPreranked.1236217933710", "MCF7_TreatedVsUntreated/c2_all.GseaPreranked.1236217933720", c("T47D_", "MCF7_"), c("", ""))
export.gsea.edb(edb.merged, "./merged")
}
}
\author{
  Mark Cowley, 2009-10-12
}
\references{
  Omega hat XML: \url{http://www.omegahat.org/RSXML/} Quick
  Guide to the XML package:
  \url{http://www.omegahat.org/RSXML/Tour.pdf}
}
\seealso{
  \code{\link[XML]{xmlTreeParse}},
  \code{\link{import.gsea}},
  \code{\link{gsea.rename.genesets.edb}}
}
\keyword{file}
\keyword{IO}

