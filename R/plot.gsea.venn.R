#' Compare the significant genesets from 2 or 3 GSEA runs using Venn diagrams
#' 
#' @param \dots either 2 or 3 GSEA objects by name, or as elements of gsea.list
#' @param gsea.list either 2 or 3 GSEA objects by name, or as elements of
#'   gsea.list
#' @param names a character vector of names, 1 per GSEA run
#' @param main plot title prefix, eg "Venn"
#' @param N see gsea.filter. It doesn't make sense to do a Venn diagram on all
#'   genesets, so you need to filter to find the significant ones.
#' @param NES see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param FDR see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param P see gsea.filter. It doesn't make sense to do a Venn diagram on all
#'   genesets, so you need to filter to find the significant ones.
#' @param FWER see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param direction see gsea.filter. It doesn't make sense to do a Venn
#'   diagram on all genesets, so you need to filter to find the significant
#'   ones.
#' @return a 2D or 3D venn diagram
#' @author Mark Cowley, 2010-10-15
#' @export
#' @importFrom mjcgraphics plot.venn
plot.gsea.venn <- function(..., gsea.list=NULL, names=NULL, main="",
	N=NULL, NES=NULL, FDR=NULL, P=NULL, FWER=NULL, direction=c("either", "up", "down")
	) {
	if( is.null(gsea.list) )
		gsea.list <- list(...)
	if( !length(gsea.list) %in% c(2,3) ) {
		stop("Need to supply 2 or 3 GSEA objects, either by name, or as elements of gsea.list.\n")
	}
	gmx <- sapply(gsea.list, function(x) basename(x$rpt$gmx))
	if( !alleq(gmx) )
		stop("You did not compare all GSEA runs to the same gmx/gmt file.\n")
	
	if( is.null(names) )
		names <- LETTERS[1:length(gsea.list)]

	pop <- unique(unlist(sapply(gsea.list, function(x) x$tt$NAME)))
	gsea.list <- lapply(gsea.list, gsea.filter, N=N, NES=NES, FDR=FDR, P=P, FWER=FWER, direction=direction)
	genesets <- lapply(gsea.list, function(x) x$tt$NAME)
	if( length(genesets) == 2 )
		plot.venn(genesets[[1]], genesets[[2]], names=names, population=pop, main=main)
	else
		plot.venn(genesets[[1]], genesets[[2]], genesets[[3]], names=names, population=pop, main=main)
}


#' Print a Venn diagram of overlaps between genesets from >= 2 GSEA runs
#' 
#' @param \dots >=2 GSEA results; or
#' @param gsea.list a list of GSEA objects
#' @param names the GSEA run names for labelling the Venn diagram
#' @param N see gsea.filter. It doesn't make sense to do a Venn diagram on all
#'   genesets, so you need to filter to find the significant ones.
#' @param NES see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param FDR see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param P see gsea.filter. It doesn't make sense to do a Venn diagram on all
#'   genesets, so you need to filter to find the significant ones.
#' @param FWER see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param direction see gsea.filter. It doesn't make sense to do a Venn
#'   diagram on all genesets, so you need to filter to find the significant
#'   ones.
#' @return none. prints a venn diagram
#' @author Mark Cowley
#' @export
print.gsea.venn <- function(..., gsea.list=NULL, names=NULL,
	N=NULL, NES=NULL, FDR=NULL, P=NULL, FWER=NULL, direction=c("either", "up", "down")
	) {
	if( is.null(gsea.list) )
		gsea.list <- list(...)
	if( !length(gsea.list) %in% c(2,3) ) {
		stop("Need to supply 2 or 3 GSEA objects, either by name, or as elements of gsea.list.\n")
	}
	gmx <- sapply(gsea.list, function(x) basename(x$rpt$gmx))
	if( !alleq(gmx) )
		stop("You did not compare all GSEA runs to the same gmx/gmt file.\n")
	
	if( is.null(names) )
		names <- LETTERS[1:length(gsea.list)]

	pop <- unique(unlist(sapply(gsea.list, function(x) x$tt$NAME)))
	gsea.list <- lapply(gsea.list, gsea.filter, N=N, NES=NES, FDR=FDR, P=P, FWER=FWER, direction=direction)
	genesets <- lapply(gsea.list, function(x) x$tt$NAME)
	if( length(genesets) == 2 )
		genesets[[3]] <- NULL

	# c <- gtools::combinations(length(genesets),2)
	c <- choose.pairs(length(genesets))
	for(i in 1:nrow(c)) {
		cat(names[ c[i,1] ], "vs", names[ c[i,2] ], "\n")
		print(print.venn(genesets[[ c[i,1] ]], genesets[[ c[i,2] ]]))
	}
}


#' Plot a Venn diagram of overlaps between genesets from >= 2 GSEA runs across
#' a range of thresholds.
#' By default, this will plot 4 venn diagrams: top 50 up, top 50 down, FDR<0.25 up,
#' FDR<0.25 down
#' 
#' @param \dots >=2 GSEA results; or
#' @param gsea.list a list of GSEA objects
#' @param names the GSEA run names for labelling the Venn diagram
#' @param file the path to the output pdf file, or \code{NULL}.
#' @param main the plot title
#' @param N see gsea.filter. It doesn't make sense to do a Venn diagram on all
#'   genesets, so you need to filter to find the significant ones.
#' @param FDR see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @param P see gsea.filter. It doesn't make sense to do a Venn diagram on all
#'   genesets, so you need to filter to find the significant ones.
#' @param FWER see gsea.filter. It doesn't make sense to do a Venn diagram on
#'   all genesets, so you need to filter to find the significant ones.
#' @return none. prints a venn diagram
#' @author Mark Cowley
#' @export
plot.gsea.venn.auto <- function(..., gsea.list=NULL, names=NULL, file=NULL, main="",
	N=c(50), FDR=c(0.25), P=NULL, FWER=NULL ) {

	if( !is.null(file) ) {
		pdf.A4(file)
		par(mfcol=c(2, length(N)+length(FDR)+length(P)+length(FWER)), mar=c(1,1,3,1), pty="s")
		on.exit(dev.off())
	}
	if( is.null(gsea.list) )
		gsea.list <- list(...)

	# N
	for(n in N) { # skips the loop if N is NULL
		for(direction in c("up", "down")) {
			tmp.main <- paste(main, " - ", direction, " - N=", n, sep="")
			plot.gsea.venn(gsea.list=gsea.list, names=names, N=n, direction=direction, main=tmp.main)
		}
	}
	# FDR
	for(fdr in FDR) {
		for(direction in c("up", "down")) {
			tmp.main <- paste(main, " - ", direction, " - FDR<", fdr, sep="")
			plot.gsea.venn(gsea.list=gsea.list, names=names, FDR=fdr, direction=direction, main=tmp.main)
		}
	}
	# P
	for(p in P) {
		for(direction in c("up", "down")) {
			tmp.main <- paste(main, " - ", direction, " - P<", p, sep="")
			plot.gsea.venn(gsea.list=gsea.list, names=names, P=p, direction=direction, main=tmp.main)
		}
	}
	# FWER
	for(fwer in FWER) {
		for(direction in c("up", "down")) {
			tmp.main <- paste(main, " - ", direction, " - FWER<", fwer, sep="")
			plot.gsea.venn(gsea.list=gsea.list, names=names, FWER=fwer, direction=direction, main=tmp.main)
		}
	}
}
