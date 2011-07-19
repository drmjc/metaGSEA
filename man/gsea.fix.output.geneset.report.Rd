\name{gsea.fix.output.geneset.report}
\alias{gsea.fix.output.geneset.report}
\title{fix an individual GSEA GeneSet html report.}
\usage{gsea.fix.output.geneset.report(htmlfile, symbol2description, addHeatMap=TRUE, debug=FALSE)
}
\description{fix an individual GSEA GeneSet html report.}
\value{none. fixes html files.}
\author{Mark Cowley}
\seealso{\code{\link{gsea.fix.output}}}
\arguments{\item{htmlfile}{the path to the gsea report file}
\item{symbol2description}{a 2 column data.frame mapping from symbols to description}
\item{addHeatMap}{logical: generate & include a heatmap into the html page}
\item{debug}{logical: debugging mode?}
}

