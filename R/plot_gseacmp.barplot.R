#' Taking a gseacmp object (see gsea.compare.runs.1gmt), create barplots of
#' the NES scores for each geneset
#' 
#' @param x a gseacmp object (see gsea.compare.runs.1gmt)
#' @param fdr.thresh Only keep genesets where at least 1 result is < this
#'   threshold.
#' @param do.par setup the plotting region? defaults to TRUE
#' @param sub see ?title
#' @param max.label.length max number of characters in a name. Middle
#'   characters will be replaced with ... to make up a name of at most
#'   max.label.length. Hint: choose an odd number.
#' @param col the barplot's background colour. ?see barplot
#' @param las see ?par
#' @param nrow how many plots per page. Passed to par(mfrow) if do.par=TRUE
#' @param ncol how many plots per page. Passed to par(mfrow) if do.par=TRUE
#' @param legend.pos Where to put the legend ONLY on the first plot. Try "top"
#'   or "bottom"
#' @param label.bars stars: statistically significant bars are labelled using
#'   stars. See label.fun
#' @param FDR the FDR value is written on the bar, only if FDR < 0.25
#' @param none no FDR labelling is done.
#' @param label.fun Function used to work out how many stars to plot for each
#'   bar. Defaults to pvalue.stars.
#' @author Mark Cowley, 2009-03-23
#' @export
plot_gseacmp.barplot <- function(x, fdr.thresh=0.05, do.par=TRUE, sub="", max.label.length=13, col="lightblue", las=2, nrow=3, ncol=4, legend.pos="none", label.bars=c("stars", "FDR", "none"), label.fun=pvalue.stars) {
	if( do.par ) {
		opar=par(no.readonly=TRUE)
		on.exit(par(opar))
		par(mfrow=c(nrow, ncol), mar=c(7,4,3,1)+0.1, las=las, cex.main=0.8, cex.lab=0.8)
	}
	
	label.bars <- label.bars[1]
	
	col <- recycle(col, ncol(x))
	
	YLIM=c(-3.1,3.1)
	WIDTH=2
	SPACE=0.5

	# filter out the poor genesets in all comparisons.
	if( fdr.thresh < 1.0 ) {
		minFDR <- apply(x[,grep("FDR", colnames(x))], 1, min, na.rm=TRUE)
		idx <- which(minFDR < fdr.thresh)
		x <- x[idx, ]
	}
	
	NES <- x[,grep("NES", colnames(x))]
	FDR <- x[,grep("FDR", colnames(x))]
	colnames(NES) <- sub("NES.", "", colnames(NES))
	colnames(FDR) <- sub("FDR.q.val.", "", colnames(FDR))

	NES[is.na(NES)] <- 0
	FDR[is.na(FDR)] <- 1

	# rename the GSEA comparison names if they are too long.
	if( !is.null(max.label.length) && !is.na(max.label.length) ) {
		colname.len <- nchar(colnames(NES))
		if( any(colname.len > max.label.length) ) {
			idx <- which(colname.len > max.label.length)
			cn <- colnames(NES)
			start.len <- floor((max.label.length-3)/2)
			cn[idx] <- paste(substr(cn[idx], 1, start.len),  "...", 
							 substr(cn[idx], colname.len[idx]-start.len+1, colname.len[idx]), sep="")
			colnames(NES) <- colnames(FDR) <- cn
		}
	}
	
	# 1 barplot per geneset.
	for(i in 1:nrow(x)) {
		nes.i <- as.numeric(NES[i,])
		fdr.i <- as.numeric(FDR[i,])
		
		barplot(nes.i, ylab="NES", main=x$NAME[i], ylim=YLIM, names.arg=colnames(NES), col=col, width=WIDTH, space=SPACE)
		if( sub != "" )
			title(sub=sub)
		abline(h=0)
		box()
		if( label.bars != "none" ) {
			bar.x <- barplot.getx(nes.i, width=WIDTH, space=SPACE)
			bar.y <- barplot.gety(nes.i, space=0.02)
			if( label.bars == "stars" ) {
				stars <- label.fun(fdr.i, c(0.05, 0.01, 0.001))
				# text(bar.x, bar.y, stars)
				text(bar.x, nes.i/2, stars, cex=2, col="white")
			}
			else if( label.bars == "FDR" && fdr.thresh < 0.25 ) {
				text(bar.x, nes.i/2, prettyNum(fdr.i), cex=0.8, col="black")
			}
		}
		
		if( i == 1 && legend.pos != "none" && label.bars == "stars" ) {
			legend(legend.pos, paste(c("*", "**", "***"), c(0.05, 0.01, 0.001), sep=" = P<"), horiz=TRUE, inset=0.02, cex=0.7, x.intersp=0.5)
		}
	}
}
# CHANGELOG
# 2010-10-26: name changed from plot_gseacmp.barplot
