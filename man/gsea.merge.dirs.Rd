\name{gsea.merge.dirs}
\alias{gsea.merge.dirs}
\title{Merge 2 gsea runs, by referencing only the GSEA directories.}
\usage{
  gsea.merge.dirs(gsea.dir1, gsea.dir2, prefix, suffix,
    outdir)
}
\arguments{
  \item{gsea.dir1}{paths to two GSEA directories (the ones
  that contain index.html)}

  \item{gsea.dir2}{paths to two GSEA directories (the ones
  that contain index.html)}

  \item{prefix}{both should be character(2) indicating how
  you want the genesets to be renamed so that they reflect
  which experimental run they came from.}

  \item{suffix}{both should be character(2) indicating how
  you want the genesets to be renamed so that they reflect
  which experimental run they came from.}

  \item{outdir}{The path to an edb directory that will be
  created.}
}
\description{
  This only merges the contents of the edb directory.
  Merging the individual html geneset reports is very
  complex... Why? so you can examine them both together in
  the Leading Edge analysis tool of course! See also
  gsea.merge
}
\author{
  Mark Cowley, 2009-10-06
}

