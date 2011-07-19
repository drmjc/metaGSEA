#' modifiction of plclust for plotting hclust objects *in colour*!
#'
#' @param hclust an hclust object
#' @param lab a character vector of labels of the leaves of the tree
#' @param lab.col colour for the labels; \code{NA}=default device foreground colour
#' @param hang see \code{\link{plclust}}
#' @param \dots args passed to \code{\link{text}}
#' @return none. Generates A display of hierarchical cluster with coloured leaf labels. 
#' @author Eva Chan, 2009
#' @export
plclust_in_colour <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, ... ) {

	y <- rep(hclust$height,2)
	x <- as.numeric(hclust$merge)

	y <- y[which(x<0)]
	x <- x[which(x<0)]

	x <- abs(x)

	y <- y[order(x)]
	x <- x[order(x)]

	plot( hclust, labels=FALSE, hang=hang, ... )
	text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), 
		  labels=lab[hclust$order], col=lab.col[hclust$order], 
		  srt=90, adj=c(1,0.5), xpd=NA, ... )
}



# 
# plclust_in_colour2 <- function( hclust, lab=hclust$labels, col=rep(1,length(hclust$labels)), hang=0.1, ... ) {
# 	## modifiction of plclust for plotting hclust objects *in colour*!
# 	## Copyright Eva KF Chan 2009
# 	## Parameters:
# 	##	hclust: hclust object
# 	##	lab:		a character vector of labels of the leaves of the tree
# 	##	lab.col:	colour for the labels; NA=default device foreground colour
# 	##	hang:	as in hclust & plclust
# 	## Side effect:
# 	##	A display of hierarchical cluster with coloured leaf labels. 
# 
# 	y <- rep(hclust$height,2)
# 	x <- as.numeric(hclust$merge)
# 
# 	y <- y[which(x<0)]
# 	x <- x[which(x<0)]
# 
# 	x <- abs(x)
# 
# 	y <- y[order(x)]
# 	x <- x[order(x)]
# 
# 	plclust( hclust, labels=lab, hang=hang, ... )
# 	text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), 
# 		  labels=rep("\255", length(hclust$height)), col=col[hclust$order], 
# 		  srt=90, adj=c(1,0.5), xpd=NA )
# }
# 
# 
# plclust_in_colour3 <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, ... ) {
# 	## modifiction of plclust for plotting hclust objects *in colour*!
# 	## Copyright Eva KF Chan 2009
# 	## Parameters:
# 	##	hclust: hclust object
# 	##	lab:		a character vector of labels of the leaves of the tree
# 	##	lab.col:	colour for the labels; NA=default device foreground colour
# 	##	hang:	as in hclust & plclust
# 	## Side effect:
# 	##	A display of hierarchical cluster with coloured leaf labels. 
# 
# 	plclust2 <- function (tree, hang = 0.1, unit = FALSE, level = FALSE, hmin = 0, 
# 	    square = TRUE, labels = NULL, plot. = TRUE, axes = TRUE, 
# 	    frame.plot = FALSE, ann = TRUE, main = "", sub = NULL, xlab = NULL, 
# 	    ylab = "Height") {
# 	    if (!missing(level) && level) 
# 	        .NotYetUsed("level", error = FALSE)
# 	    if (!missing(hmin) && hmin != 0) 
# 	        .NotYetUsed("hmin", error = FALSE)
# 	    if (!missing(square) && !square) 
# 	        .NotYetUsed("square", error = FALSE)
# 	    if (!missing(plot.) && !plot.) 
# 	        .NotYetUsed("plot.", error = TRUE)
# 	    if (!missing(hmin)) 
# 	        tree$height <- pmax(tree$height, hmin)
# 	    if (unit) 
# 	        tree$height <- rank(tree$height)
# 	    plot.hclust2(x = tree, labels = labels, hang = hang, axes = axes, 
# 	        frame.plot = frame.plot, ann = ann, main = main, sub = sub, 
# 	        xlab = xlab, ylab = ylab)
# 	}
# 	plot.hclust2 <- function (x, labels = NULL, hang = 0.1, axes = TRUE, frame.plot = FALSE, 
# 	    ann = TRUE, main = "Cluster Dendrogram", sub = NULL, xlab = NULL, 
# 	    ylab = "Height", lab.col=rep("black", length(x$labels), ...) {
# 	    merge <- x$merge
# 	    if (!is.matrix(merge) || ncol(merge) != 2) 
# 	        stop("invalid dendrogram")
# 	    if (any(as.integer(merge) != merge)) 
# 	        stop("'merge' component in dendrogram must be integer")
# 	    storage.mode(merge) <- "integer"
# 	    n <- nrow(merge)
# 	    height <- as.double(x$height)
# 	    labels <- if (missing(labels) || is.null(labels)) {
# 	        if (is.null(x$labels)) 
# 	            paste(1L:(n + 1))
# 	        else as.character(x$labels)
# 	    }
# 	    else {
# 	        if (is.logical(labels) && !labels) 
# 	            character(n + 1)
# 	        else as.character(labels)
# 	    }
# 	    plot.new()
# 	    .Internal(dend.window(n, merge, height, hang, labels, ...))
# 	    .Internal(dend(n, merge, height, order(x$order), hang, labels, 
# 	        ...))
# 	    if (axes) 
# 	        axis(2, at = pretty(range(height)))
# 	    if (frame.plot) 
# 	        box(...)
# 	    if (ann) {
# 	        if (!is.null(cl <- x$call) && is.null(sub)) 
# 	            sub <- paste(deparse(cl[[1L]]), " (*, \"", x$method, 
# 	                "\")", sep = "")
# 	        if (is.null(xlab) && !is.null(cl)) 
# 	            xlab <- deparse(cl[[2L]])
# 	        title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
# 	            ...)
# 	    }
# 	    invisible()
# 	}
# 	
# 	plclust2(hclust, hang=hang, ...)
# 	y <- rep(hclust$height,2)
# 	x <- as.numeric(hclust$merge)
# 
# 	y <- y[which(x<0)]
# 	x <- x[which(x<0)]
# 
# 	x <- abs(x)
# 
# 	y <- y[order(x)]
# 	x <- x[order(x)]
# 
# 	plot( hclust, labels=FALSE, hang=hang, ... )
# 	text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), 
# 		  labels=lab[hclust$order], col=lab.col[hclust$order], 
# 		  srt=90, adj=c(1,0.5), xpd=NA, ... )
# }
