#' GSEA Leading Edge analysis: Barplot of number of genesets each gene occured
#' in.
#' 
#' This tool emulates the barplot found in the GSEA Leading Edge Tool.
#' Typically, the user selects a subset of
#' genesets of interest (use \code{\link{gsea.filter}}, or 
#' \code{\link{gsea.split}} to do this) and then
#' wants to find genes that are present
#' in lots of those genesets. This tool gives you a barplot indicating the
#' number of genesets that each gene was
#' found in.
#' 
#' @param x a GSEA object. Usually this has been subset to a managable number
#'   of genesets
#' @param min.count filter out genes found in only a few genesets. min.count=2
#'   is useful in practice.
#' @param xlab see \code{\link{plot}}
#' @param ylab see \code{\link{plot}}
#' @param main see \code{\link{plot}}
#' @param horiz logical: if TRUE, then the bars will be horizontal; FALSE they are
#'   vertical (default).
#' @param col the colour of the bars
#' @param cex.names the character expansion factor of the gene names. see \code{\link{plot}}
#' @param \dots additional arguments passed to barplt
#' @return invisibly returns the named vector of counts.
#' @author Mark Cowley, 2009-10-28
#' @seealso \code{\link{barplot}}
#' @export
plot.gsea.leadingedge.barplot <- function(x, min.count=0, 
	xlab="", ylab="Number Of Gene Sets", main="Leading Edge Barplot", 
	horiz=FALSE, col="blue", cex.names=0.6, ...) {

	# genes <- lapply(x, "[", "leading.edge")
	counts <- table(unlist(x$leading.edge))
	counts <- sort(counts, decreasing=TRUE)
	if( min.count > 0 ) {
		counts <- counts[counts >= min.count]
	}
	
	if( length(counts) == 0 ) {
		plot.blank(main=main, box=TRUE, message="**** too few genes pass filters ****")
	}
	else {
		opar <- par(no.readonly=TRUE)
		on.exit(par(opar))
		if( horiz ) {
			tmp <- xlab
			xlab <- ylab
			ylab <- tmp
			par(mar=c(5,8,4,2), las=1)
			counts <- rev(counts)
		}
		else {
			par(mar=c(8,4,4,2), las=2, mgp=c(3,1,0))
		}

		barplot(counts, main=main, xlab=xlab, ylab=ylab, border=NA, horiz=horiz, width=1, space=10, col=col, cex.names=cex.names, ...)

		if( horiz ) {
			vgrid(col="lightgrey", lty="dashed")
		}
		else {
			hgrid(col="lightgrey", lty="dashed")
		}
		box()
		if( min.count > 0 ) {
			mtext(side=3, outer=FALSE, adj=0.99, paste("min count >=", min.count), line=0, font=3, las=1)
		}
	}

	invisible(counts)
}


#' Function to automatically plot barplot's of filtered GSEA data.
#' 
#' @note Usually, I filter by up or down, by top 50, and by FDR<0.25. This generates 4 plots
#' @param x a GSEA list
#' @param file the path to an output PDF file. if \code{NULL}, then plots are printed to current device.
#' @param N if not \code{NULL}, then choose the top N most significant genesets
#' @param FDR if not \code{NULL}, then choose those genesets with FDR < the
#'   threshold
#' @param P if not \code{NULL}, then choose those genesets with nominal P <
#'   the threshold
#' @param FWER if not \code{NULL}, then choose those genesets with FWER < the
#'   threshold
#' @param min.count filter out genes found in only a few genesets. min.count=2
#'   is useful in practice.
#' @param horiz logical: horizontal or vertical plot?
#' @param col the colour of the bars
#' @param cex.names the character expansion factor of the gene names. see \code{\link{plot}}
#' @param \dots arguments passed to \code{\link{plot.gsea.leadingedge.barplot}}
#' @param write.delim logical: if \code{TRUE}, then write the data plotted as a tab
#' delimited txt file, named \code{file + ".txt"}. Default is \code{FALSE}. If 
#' \code{file=NULL}, then write.delim will be set to \code{FALSE}.
#' @author Mark Cowley, 2010-10-14
#' @export
#' @return none. Generates a number of plots, either to a pdf if \code{file!=NULL}, or 
#'   \code{\link{dev.cur}}
#' 
plot.gsea.leadingedge.barplot.auto <- function(x, file=NULL, 
	N=c(50), FDR=c(0.25), P=NULL, FWER=NULL, 
	min.count=2,
	horiz=FALSE, col="blue", cex.names=0.6, write.delim=FALSE, ...) {

	if( "leading.edge" %in% names(x) )
		x <- list(gsea.collection=x)

	if( !is.null(file) ) {
		pdf.A4(file)
		on.exit(dev.off())
	}
	else {
		write.delim <- FALSE
	}

	collections <- names(x)
	for( collection in collections ) {
		cat(".")
		
		main <- "Leading Edge Barplot"
		if( collection != "gsea.collection")
			main <- paste(main, collection, sep=" - ")
			
		# N
		for(n in N) { # skips the loop if N is NULL
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - topN=", n, sep="")
				tmp <- gsea.filter(x[[collection]], N=n, direction=direction)
				if( length(tmp$leading.edge) > 0 ) {
					tmp.counts <- plot.gsea.leadingedge.barplot(tmp, main=tmp.main, min.count=min.count, horiz=horiz, col=col, cex.names=cex.names, ...)
					if( write.delim ) {
						f <- paste(file, " - ", tmp.main, ".txt", sep="")
						tmp.counts <- data.frame(GeneSymbol=names(tmp.counts), Frequency=tmp.counts)
						write.delim(tmp.counts, f)
					}
				}
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
		# FDR
		for(fdr in FDR) {
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - FDR<", fdr, sep="")
				tmp <- gsea.filter(x[[collection]], FDR=fdr, direction=direction)
				if( length(tmp$leading.edge) > 0 ) {
					tmp.counts <- plot.gsea.leadingedge.barplot(tmp, main=tmp.main, min.count=min.count, horiz=horiz, col=col, cex.names=cex.names, ...)
					if( write.delim ) {
						f <- paste(file, " - ", tmp.main, ".txt", sep="")
						tmp.counts <- data.frame(GeneSymbol=names(tmp.counts), Frequency=tmp.counts)
						write.delim(tmp.counts, f)
					}
				}
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
		# P
		for(p in P) {
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - P<", p, sep="")
				tmp <- gsea.filter(x[[collection]], P=p, direction=direction)
				if( length(tmp$leading.edge) > 0 ) {
					tmp.counts <- plot.gsea.leadingedge.barplot(tmp, main=tmp.main, min.count=min.count, horiz=horiz, col=col, cex.names=cex.names, ...)
					if( write.delim ) {
						f <- paste(file, " - ", tmp.main, ".txt", sep="")
						tmp.counts <- data.frame(GeneSymbol=names(tmp.counts), Frequency=tmp.counts)
						write.delim(tmp.counts, f)
					}
				}
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
		# FWER
		for(fwer in FWER) {
			for(direction in c("up", "down")) {
				tmp.main <- paste(main, " - ", direction, " - FWER<", fwer, sep="")
				tmp <- gsea.filter(x[[collection]], FWER=fwer, direction=direction)
				if( length(tmp$leading.edge) > 0 ) {
					tmp.counts <- plot.gsea.leadingedge.barplot(tmp, main=tmp.main, min.count=min.count, horiz=horiz, col=col, cex.names=cex.names, ...)
					if( write.delim ) {
						f <- paste(file, " - ", tmp.main, ".txt", sep="")
						tmp.counts <- data.frame(GeneSymbol=names(tmp.counts), Frequency=tmp.counts)
						write.delim(tmp.counts, f)
					}
				}
				else
					plot.blank(main=tmp.main, box=TRUE, message="**** too few genesets pass thresholds ****")
			}
		}
	}
	cat("\n")
}
# CHANGELOG
# 2010-10-14: v1
# 2012-02-19: added write.delim option.