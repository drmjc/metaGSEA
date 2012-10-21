\name{gsea.topTable.sort}
\alias{gsea.topTable.sort}
\title{sort a GSEA topTable}
\usage{
  gsea.topTable.sort(tt,
    sort.by = c("NES", "absNES", "FDR", "P", "SIZE"))
}
\arguments{
  \item{tt}{a GSEA topTable}

  \item{sort.by}{one of c("NES", "absNES", "FDR", "P",
  "SIZE"), controlling how to sort the top table.}
}
\value{
  a sorted GSEA topTable
}
\description{
  sort a GSEA topTable
}
\author{
  Mark Cowley, 2012-10-12
}

