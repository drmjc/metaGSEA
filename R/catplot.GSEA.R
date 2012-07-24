#' Function to generate CAT plots for >= 2 GSEA comparisons
#' 
#' This function expects all GSEA runs to be generated vs the \emph{same} gmt file,
#' and with no filtering done. Aside from the filtering upon geneset size that GSEA
#' performs, there should be a large overlap in terms of genesets
#' 
#' @param \dots at least 2 GSEA objects, or:
#' @param gsea.list an optional list of >= 2 GSEA objects. This overrides \dots
#' @param labels the name of of each GSEA run. must be same length as \dots, or
#'   gsea.list.
#' @param legend.pos where to add the names legend. default = \dQuote{bottomright}, 
#' set to \dQuote{none} to disable
#' @param proportion if \code{TRUE}, x-axis is proportion of all genesets, if \code{FALSE}
#'   (default), the number of genesets is plotted.
#' @param main the plot title.
#' @param FDRthresh If you'd like the \code{pch} to change depending on which of A or
#'   B is significant, then choose a threshold of 0.25 for example. Pch will
#'   change depending on whether A and B are both significant, only one is
#'   significant, or neither are.
#' @param method \dQuote{1-vs-all} means compare the first GSEA run to all subseqent
#'   runs; \dQuote{all-pairs} plots all pairwise comparisons of genesets.
#' @author Mark Cowley, 2008-12-08 (edited date)
#' @export
catplot.GSEA <- function(..., gsea.list=NULL, labels=LETTERS, legend.pos="bottomright", proportion=FALSE, FDRthresh=1.01, main="", method=c("1-vs-all", "all-pairs")) {
	# parse the ..., which could be a N objects, or a list of gsea.lists.
	if( is.null(gsea.list) )
		gsea.list <- list(...)
	if( length(gsea.list)==1 ) {
		gsea.list <- gsea.list[[1]]
		if( is.list(gsea.list) && !is.gsea.list(gsea.list) && is.gsea.list(gsea.list[[1]]) ) {
			#nothing to do, each element is a gsea list
			random <- 0
			# names(gsea.list) <- labels
		}
		else if( is.gsea.list(gsea.list) ) {
			# N gsea.list should be N GSEA's vs 1 gmt...
			gmx <- sapply(gsea.list, function(x) basename(x$rpt$gmx))
			if( !alleq(gmx) )
				stop("You did not compare all GSEA runs to the same gmx/gmt file.\n")
			# ... convert into a gmt list of length 1.
			res <- as.list(rep(NA,length(gsea.list)))
			for(i in 1:length(gsea.list)) {
				res[[i]] <- gsea.list[[i]]
				# names(res[[i]]) <- gmx[i]
			}
			names(res) <- names(gsea.list)
			gsea.list <- res
		}
		else {
			stop("Arguments should either be (1) a single gsea list or multiple gsea lists of GSEA runs vs the same gmx/gmt file, or (2) a list of gsea lists of GSEA runs vs muliple gmx/gmt files.\n")
		}
	}

	# work out which sets to compare to
	method <- method[1]
	if( method == "all-pairs")
		combinations <- combinations(length(gsea.list), 2)
	else if( method == "1-vs-all" )
		combinations <- cbind(1, 2:length(gsea.list))
	colours <- 1:nrow(combinations)
	
	for(i in 1:nrow(combinations)) {
		lhs <- combinations[i,1]
		rhs <- combinations[i,2]
		catplot.GSEA.2way(gsea.list[[lhs]], gsea.list[[rhs]], 
			add=!(lhs==1&rhs==2), col=colours[i], proportion=proportion, 
			FDRthresh=FDRthresh, main=main)
	}

	if( !is.null(labels) ) {
		# labels <- recycle(labels, length(gsea.list))
		# legend.text <- paste(rep(labels[1], length(gsea.list)-1), labels[2:length(gsea.list)], sep=" vs ")
		legend.text <- paste(labels[combinations[,1]], labels[combinations[,2]], sep=" -vs- ")
		
		legend(legend.pos, legend.text, pch=1, col=1:length(legend.text), inset=0.01)
	}
}

#' is \code{x} a gsea list, that is, a list containing gsea results?
#' @param x an object
#' @author Mark Cowley, 2008-05-13
#' @return logical: \code{TRUE} if x is a list, and contains only GSEA objects, \code{FALSE} otherwise.
#' @export
is.gsea.list <- function(x) {
	is.list(x) && all(sapply(x, is.gsea))
}




# catplot.GSEA.2way <- function(x, y, col=1, proportion=FALSE, add=FALSE, GSEAsets=NULL, FDRthresh=1.01, main="") {
# 
# 	if( is.data.frame(x) && is.data.frame(y) ) {
# 		x <- list(x)
# 		y <- list(y)
# 		if( is.null(GSEAsets) ) {
# 			names(x) <- names(y) <- GSEAsets <- "geneset"
# 		}
# 		else {
# 			names(x) <- names(y) <- GSEAsets <- GSEAsets[1]
# 		}
# 	}
# 	
# 	if( is.null(GSEAsets) )
# 		GSEAsets <- intersect(names(x), names(y))
# 	else if( is.numeric(GSEAsets) )
# 		GSEAsets <- names(x)[GSEAsets]
# 	
# 	if( !add ) par(mfcol=c(2, length(GSEAsets)), las=1, oma=c(0,3,3,0))
# 
# 	.plot <- function(x, y, direction, main, plot.column) {
# 		if( direction == "down" ) {
# 			x <- x[nrow(x):1,]
# 			y <- y[nrow(y):1,]
# 		}
# 		x.is.signif <- x$DIRECTION == direction & x$FDR < FDRthresh
# 		y.is.signif <- y$DIRECTION == direction & y$FDR < FDRthresh
# 
# 		if(add) par(mfg=c(ifelse(direction=="up", 1, 2), plot.column))
# 		catplot(x$NAME, y$NAME, 
# 			main=main,
# 			proportion=proportion, col=col, add=add, 
# 			xTrue=x.is.signif, yTrue=y.is.signif, 
# 			sketch=1.0)
# 		if(!add) grid()
# 		
# 	}
# 
# 	for(i in 1:length(GSEAsets)) {
# 		GSEAset <- GSEAsets[i]
# 		.plot(x[[ GSEAset ]], y[[ GSEAset ]], "up", GSEAset, i)
# 		.plot(x[[ GSEAset ]], y[[ GSEAset ]], "down", GSEAset, i)
# 	}
# 
# 	mtext(side=3, outer=TRUE, main, cex=2)
# 	mtext(side=2, outer=TRUE, at=c(0.75, 0.25), c("+NES -> -NES", "-NES -> +NES"), las=0, line=1, cex=1.5)
# }


#' Plot a CAT plot upon 2 GSEA objects.
#' 
#' This is a helper function, you should probably consider using \code{\link{catplot.GSEA}}.
#'
#' @param x a GSEA object
#' @param y a GSEA object
#' @param col the plot character foreground colour. default = 1 = \dQuote{black}
#' @param proportion if \code{TRUE}, x-axis is proportion of all genesets, if \code{FALSE}
#'   (default), the number of genesets is plotted.
#' @param add logical: add this to existing plot? if \code{FALSE}, then create a new plotting device.
#' @param FDRthresh If you'd like the \code{pch} to change depending on which of A or
#'   B is significant, then choose a threshold of 0.25 for example. Pch will
#'   change depending on whether A and B are both significant, only one is
#'   significant, or neither are.
#' @param main the plot title. set to \dQuote{} to ignore
#' @return none. makes a plot
#' @export
#' @author Mark Cowley, 2011-07-19
#'
catplot.GSEA.2way <- function(x, y, col=1, proportion=FALSE, add=FALSE, FDRthresh=1.01, main="") {
	## Function is optimised to work with a list of GSEA objects. IF x and y are just GSEA objects, then make them lists, and give them an appropriate name.
	if( is.gsea(x) && is.gsea(y) ) {
		x.gmx <- basename(x$rpt$gmx)
		y.gmx <- basename(y$rpt$gmx)
		if( x.gmx != y.gmx ) {
			stop("x and y were compared to different gmt/gmx files.\n")
		}
		x <- list(x)
		y <- list(y)
		names(x) <- names(y) <- x.gmx
	}
	else if( is.gsea.list(x) && is.gsea.list(y)) {
		x.gmx <- sapply(x, function(gsea) basename(gsea$rpt$gmx))
		y.gmx <- sapply(y, function(gsea) basename(gsea$rpt$gmx))
		if( length(x.gmx) != length(y.gmx) || !all(x.gmx == y.gmx))
			stop("The GSEA runs within x and y differ, either in their order, or length.\n")
		if( is.null(names(x.gmx)) ) names(x.gmx) <- x.gmx
		if( is.null(names(y.gmx)) ) names(y.gmx) <- y.gmx
	}
	else {
		stop("Unsupported format for x. Should be a gsea object, or a list of gsea objects.\n")
	}
	#
	# now x and y should be a list of at least 1 gsea objects.
	#####
	
	if( !add ) par(mfcol=c(2, length(x)), las=1, oma=c(0,3,3,0))

	.plot <- function(x, y, direction, main, plot.column) {
		x <- x$tt
		y <- y$tt
		# make sure tables are sorted by NES
		x <- x[order(x$NES, decreasing=TRUE), ]
		y <- y[order(y$NES, decreasing=TRUE), ]
		if( direction == "down" ) {
			x <- x[nrow(x):1,]
			y <- y[nrow(y):1,]
		}
		x.is.signif <- x$DIRECTION == direction & x$FDR < FDRthresh
		y.is.signif <- y$DIRECTION == direction & y$FDR < FDRthresh

		if(add) par(mfg=c(ifelse(direction=="up", 1, 2), plot.column))
		catplot(x$NAME, y$NAME, 
			main=main,
			proportion=proportion, col=col, add=add, 
			xTrue=x.is.signif, yTrue=y.is.signif, 
			sketch=1.0)
		if(!add) grid()
		
	}

	for(i in 1:length(x)) {
		GSEAset <- names(x)[i]
		.plot(x[[ i ]], y[[ i ]], "up", GSEAset, i)
		.plot(x[[ i ]], y[[ i ]], "down", GSEAset, i)
	}

	mtext(side=3, outer=TRUE, main, cex=2)
	mtext(side=2, outer=TRUE, at=c(0.75, 0.25), c("+NES -> -NES", "-NES -> +NES"), las=0, line=1, cex=1.5)
}



#' Compare multiple GSEA runs via CAT plots & estimate the expected overlap due to random chance.
#' 
#' This function generates a normal X vs Y CAT plot, but then randomises Y, B times to get a permuted distribution
#' of overlaps.
#' 
#' Controlling the granularity of the CAT plot with sizes and sizes.random\cr
#' \dQuote{sizes} controls the granularity of the CAT plot. 
#' If sizes=\code{NULL}, then sizes = 1,2,3,4,\dots,N, where N is the
#'   total number of genesets. This gives a very detailed CAT plot, but if this is too slow, 
#'   try setting \code{sizes=seq(0,N,5)}.
#' Same goes for sizes.random, but since there are B randomisations, this setting has more of an impact
#' upon the speed of this code. typically with just ~1400 genesets, it's fast enough to calculate on all
#' values from 1 to N, however: \code{sizes=seq(0,N,5)}, or  \code{sizes=seq(0,N,10)} will be much quicker.
#' 
#' @param \dots at least 2 GSEA objects
#' @param labels the name of of each GSEA run. must be same length as \dots, or
#'   gsea.list.
#' @param B the number of permutations. default=100
#' @param sizes Controls the granularity of the CAT plot. see details
#' @param sizes.random Controls the granulatity of the randomised CAT plot. Since this gets done B times,
#'    this setting can really impact the speed of this function. See details.
#' @param random.col the fill colour for the random background
#' @param ylim see \code{\link{par}}
#' @param proportion if \code{TRUE}, x-axis is proportion of all genesets, if \code{FALSE}
#'   (default), the number of genesets is plotted.
#' @param pch the print character. see \code{\link{par}}
#' @param legend.pos the legend position. see \code{\link{legend}}
#' @seealso \code{\link{catplot.vs.random}}, \code{\link{catplot.GSEA}}
#' @author Mark Cowley
#' @export
catplot.GSEA.random <- function(..., labels=LETTERS, B=100, sizes=NULL, sizes.random=sizes, random.col="grey", ylim=c(0,1), proportion=TRUE, pch=1, legend.pos="bottomright") {
	warning("untested since updating catplot.GSEA.2way.vs.random. MJC 2010-10-27\n")
	gsea.list <- list(...)

	if( length(gsea.list)==1 && is.list(gsea.list[[1]]) && !is.gsea.list(gsea.list[[1]]) && is.gsea.list(gsea.list[[1]][[1]]) && length(gsea.list[[1]]) >= 2 )
		# then only 1 arg was supplied, which is itself a list.
		gsea.list <- gsea.list[[1]]

	for(i in 2:length(gsea.list)) {
		catplot.GSEA.2way.vs.random(gsea.list[[1]], gsea.list[[i]], add=i>2,col=i-1, proportion=proportion, B=B, sizes=sizes, sizes.random=sizes.random, random.col=random.col, pch=pch)
	}

	if( !is.null(labels) ) {
		labels <- recycle(labels, length(gsea.list))
		legend.text <- paste(rep(labels[1], length(gsea.list)-1), labels[2:length(gsea.list)], sep=" vs ")
		legend(legend.pos, legend.text, pch=1, col=1:length(legend.text), inset=0.01)
	}
}


#' Compare 2 GSEA runs via CAT plots & estimate the expected overlap due to random chance.
#' 
#' You should probably use \code{\link{catplot.GSEA.random}} instead
#' 
#' @param x a GSEA object
#' @param y a GSEA object
#' @param B the number of permutations. default=100
#' @param sizes Controls the granularity of the CAT plot. see details
#' @param sizes.random Controls the granulatity of the randomised CAT plot. Since this gets done \code{B} times,
#'    this setting can really impact the speed of this function. See details.
#' @param random.col the fill colour for the random background
#' @param col the pch colour
#' @param proportion if \code{TRUE}, x-axis is proportion of all genesets, if \code{FALSE}
#'   (default), the number of genesets is plotted.
#' @param add logical: add to the current plot (\code{TRUE}), or create a new plot (\code{FALSE}).
#' @param plot.column if \code{add=TRUE}, which column of plots does the current plot fit 
#'     into? see mfg setting in \code{\link{par}}
#' @param pch the print character. see \code{\link{par}}
#' @param main the plot title
#' @seealso \code{\link{catplot.vs.random}}, \code{\link{catplot.GSEA}}
#' @author Mark Cowley
#' @export
catplot.GSEA.2way.vs.random <- function(x, y, B=100, sizes=NULL, sizes.random=sizes, random.col="lightgrey", col=1, proportion=FALSE, add=FALSE, plot.column=1, pch=1, main="") {

	stopifnot( is.gsea(x), is.gsea(y) )
	
	if( !add ) par(mfcol=c(2, 1), las=1, oma=c(0,0,3,0), mar=c(4,4,3,0.5))

	if( is.null(sizes) )
		sizes <- 1:length(x$tt$NAME)
	if( is.null(sizes.random) )
		sizes.random <- sizes

	if(add) par(mfg=c(1, plot.column))
	catplot.vs.random(x$tt$NAME, y$tt$NAME, 
		B=B, sizes=sizes, sizes.random=sizes.random,
		main="CAT plot overlap in up-regulated genesets",
		random.col=random.col, proportion=proportion, col=col, add=add, pch=pch
	)
	if(!add) grid()

	if(add) par(mfg=c(2, plot.column))
	catplot.vs.random(rev(x$tt$NAME), rev(y$tt$NAME), 
		B=B, sizes=sizes, sizes.random=sizes.random,
		main="CAT plot overlap in down-regulated genesets",
		random.col=random.col, proportion=proportion, col=col, add=add, pch=pch
	)
	if(!add) grid()

	mtext(side=3, outer=TRUE, main, cex=2)
	mtext(side=2, outer=TRUE, at=c(0.75, 0.25), c("up-regulated", "down-regulated"), las=0, line=1, cex=1.5)
}
# CHANGELOG
# ???: v1
# 2010-10-26; important update. Only operates on GSEA objects, NOT lists of
# GSEA objects.
#