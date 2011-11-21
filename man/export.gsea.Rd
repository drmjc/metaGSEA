\name{export.gsea}
\alias{export.gsea}
\title{Export a gsea object to the filesystem.}
\usage{
  export.gsea(x, location, as.zip = FALSE,
    overwrite = FALSE)
}
\arguments{
  \item{x}{a GSEA object.}

  \item{location}{the dir that will be created. eg:
  \dQuote{./GSEA/c2_all}, \dQuote{c2_all.zip}, and set
  \code{as.zip=TRUE}}

  \item{as.zip}{logical: create a zip file? If \code{TRUE},
  then make sure the \code{location} parameter is set to a
  zip file name, rather than a directory path.}

  \item{overwrite}{logical: overwrite an existing file or
  directory?}
}
\value{
  none. Creates either a directory, or a zip file of a
  directory. It will only contain the information that was
  loaded by import.gsea, which means that individual
  geneset html resports and pictures are not exported.
}
\description{
  This is mostly used to export the edb portion of a gsea
  run, since copying over the individual geneset reports
  etc sounds far too tedious, plus these elements are not
  stored inside the GSEA object in R.
}
\author{
  Mark Cowley, 2009-10-15
}

