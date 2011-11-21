\name{gsea.gct.select.features.columns}
\alias{gsea.gct.select.features.columns}
\title{Restrict the number of columns in a gct & cls file.}
\usage{
  gsea.gct.select.features.columns(gct.file, cls.file,
    gct.out, cls.out, columns = NULL, classes = NULL,
    cls.start.value = 0)
}
\arguments{
  \item{gct.file}{the path to a GCT file}

  \item{cls.file}{the path to a CLS file}

  \item{gct.out}{the path to a GCT output file}

  \item{cls.out}{the path to a CLS output file}

  \item{columns}{the column indices to include, or
  \code{NULL}}

  \item{classes}{the classes to include, or \code{NULL}}

  \item{cls.start.value}{most cls files start from 0, but
  that can be changed if necessary}
}
\value{
  none. creates a gct and cls file containing the specified
  subset of classes or columns.
}
\description{
  You can do this either by indexing the columns ([1:n]),
  or by specifying the names of the classes. The order of
  columns, or classes is respected.
}
\author{
  Mark Cowley, 2010-07-06
}

