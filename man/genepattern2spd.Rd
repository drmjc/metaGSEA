\name{genepattern2spd}
\alias{genepattern2spd}
\title{GenePattern to SPD}
\usage{
  genepattern2spd(gct.file, cls.files, cls.names = "Class",
    out.file = "spd.mat")
}
\arguments{
  \item{gct.file}{The path to a GenePattern formatted GCT
  file}

  \item{cls.files}{The path to a GenePattern formatted CLS
  file. Optional.}

  \item{cls.names}{What do the classes in cls.file
  represent? Tumour Type? Cell Line? Specify a short label,
  eg \dQuote{Tumour Type}}

  \item{out.file}{The path to the output file. It should
  end in .mat.}
}
\value{
  Creates a version 5 matlab file which you can then open
  in Matlab
}
\description{
  Convert from GenePattern GCT and CLS files to the input
  file required for Sample Progression Discovery (SPD) The
  input file is a matlab binary file, containing probe
  names, data matrix and samples names, and optionally
  class labels, from a cls file. This makes use of the
  R.matlab library by Henrik Bengtsson & creates v5 matlab
  files (which work fine with matlab v7.)
}
\section{TODO}{
  support >1 cls files
}
\examples{
\dontrun{
setwd("/Volumes/marcow-1/projects/SPD/analysisV1")
genepattern2spd("norm.rsn.medianNorm.gct", out.file="norm.rsn.medianNorm.SPD.mat")
genepattern2spd("norm.rsn.medianNorm.gct", "norm.rsn.medianNorm.cls", "Tumour Grade", "norm.rsn.medianNorm.SPD.mat")
}
}
\author{
  Mark Cowley, 2011-04-20
}
\seealso{
  \url{http://icbp.stanford.edu/software/SPD},
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#gct},
  \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}
}

