\name{import.gsea.rpt}
\alias{import.gsea.rpt}
\title{Import the rpt file from a GSEA, or GseaPreRanked analysis run...}
\usage{import.gsea.rpt(x)
}
\description{Import the rpt file from a GSEA, or GseaPreRanked analysis run}
\details{If the GSEA was run on a GenePattern server, then some of the values (eg
chip, gmx, cls/rnk) will point to
locations that are unavailable on the end users machine. chip and gmx files
are replaced with URL's to Broad Institute. cls, rnk, out, file references
are changed to locations relative to their current location on the end
users hard drive.}
\author{Mark Cowley, 2009-04-28}
\value{a list representation of a GSEA rpt file, where names = the gsea flag,
values = the setting.}
\arguments{\item{x}{the path to either a rpt file, or the top-level dir that contains
index.html.}
}

