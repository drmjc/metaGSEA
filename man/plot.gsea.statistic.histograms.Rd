\name{plot.gsea.statistic.histograms}
\alias{plot.gsea.statistic.histograms}
\title{plot a panel of gsea histograms...}
\usage{plot.gsea.statistic.histograms(x, main="", which=c("p", "q", "FWER"))
}
\description{plot a panel of gsea histograms
From a GSEA object, plot a histogram of the nominal p values,
FDR q values, or FWER p values.}
\value{none. creates a panel of histograms.}
\author{Mark Cowley, 2011-10-20}
\arguments{\item{x}{a gsea object}
\item{main}{the plot titles}
\item{which}{one or more of \dQuote{p}, \dQuote{q}, or \dQuote{FWER}, 
corresponding to these columns in a GSEA report: \dQuote{NOM.p.val}, 
\dQuote{FDR.q.val}, \dQuote{FWER.p.val}, respectively.}
}

