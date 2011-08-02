\name{gsea.hardwire.rpt.paths}
\alias{gsea.hardwire.rpt.paths}
\title{Rewire the file locations referenced within a GSEA rpt file.}
\usage{gsea.hardwire.rpt.paths(dir, chip.dir="/pwbc/data/GSEA/chips",
    gmx.dir="/pwbc/data/GSEA/gmt/v2.5")
}
\description{Rewire the file locations referenced within a GSEA rpt file.}
\details{Why? If you run a GSEA on a GenePattern server, you can end up with files
in a rpt file that you can
no longer access on the machine where you're analysing the GSEA results...
Useful if you have a directory full of chip files, and gmx files.}
\value{none. the rpt file will be modified.}
\author{Mark Cowley, 2009-12-11}
\arguments{\item{dir}{a GSA directory. The rpt file within this dir will be modified.}
\item{chip.dir}{the directory containing your GSEA .chip files}
\item{gmx.dir}{the directory containing your GSEA .gmx/.gmt files}
}

