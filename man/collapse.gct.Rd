\name{collapse.gct}
\alias{collapse.gct}
\title{Collapse a GCT object...}
\usage{collapse.gct(gct, chip, rnk, method=c("var", "mean", "median"), reverse=FALSE, filter=FALSE)
}
\description{Collapse a GCT object}
\details{It's common for microarrays to have multiple probes per gene. They tend to
represent different isoforms.
Most geneset testing is done at the gene symbol level & ignores isoforms,
so you need to choose 1 probe
for each gene.
How?
2 common approaches are to take the most abundant probe, or the most
variable probe, considered
across the cohort.
I quite like doing t-stats on each gene & selecting the best performing
probe - ie the one with the largest t-stat in
either direction. Why? On the Affy 133+2 array, there can be lots of poor
probes for each gene. If 5 probes for a gene
have these t-stats: 1.2, 0.9, 0.1, -0.1, -10; then IMO, the one that scored
-10 is the best probe, since it had a really
strong t-stat score. thus method="maxabs" combined with a rnk.file (NB NOT
implemented right now.)}
\value{A gct object with 1 row per gene symbol & now the \sQuote{probe ids} in
column 1 are actually gene symbols.}
\author{Mark Cowley, 2011-02-27}
\arguments{\item{gct}{a GCT object}
\item{chip}{a CHIP object}
\item{rnk}{[optional] path to a rnk file (eg a t-statistic for each probe,
where you want to select best probe from this score) NB currently UNUSED}
\item{method}{\dQuote{mean}, \dQuote{median} select the probe with highest 
average/median level, or \dQuote{var}: select the probe with highest variance 
across samples; \dQuote{maxabs} select the probe with the large absolute score in the rnk
file (see details).}
\item{reverse}{[default=FALSE] reverse the ordering selected by method arg.
so instead of most variable, it would be least variable.}
\item{filter}{Filter out (ie exclude) those probes that don't have a gene
symbol (as determined by probes that have a gene symbol of \code{NA}, \dQuote{},
\dQuote{---}, or \dQuote{NA}.)}
}

