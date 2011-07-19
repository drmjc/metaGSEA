#' Export a cls file for GenePattern.
#' Export a cls file for GenePattern. cls files represent the class labels - a
#' way of telling GenePattern which experimental group each array belongs to.
#' They can be discrete, or continuous.
#' 
#' In the discrete case, each cls can contain only only parameter, thus the
#' argument must be a vector of labels. If you'd like control over which is
#' gets assigned to 0, 1, 2... then use pass in a factor of your making.
#' 
#' In the continuous case, you can have >= 1 parameter.  either pass in a
#' numeric vector; or a matrix/data.frame where each row represents a
#' different parameter; or a list where each element is a vector, and each
#' element represents a parameter.
#' 
#' @param labels the labels. see description.
#' @param file the output file name.
#' @param continuous are the labels continuous or discrete?
#'   \code{TRUE}/\code{FALSE}, respectively.
#' @param start.value most cls files start from 0, however you can change this if you like.
#' @return none. creates a cls file.
#' @note update 2009-07-08: cls labels can't have spaces in them.
#' @author Mark Cowley, 2009-06-26
#' @seealso \code{\link{import.gsea.cls}}
#' @references
#'   \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}
#' @keywords file IO
#' @examples
#' cls <- c(rep("wt",5), rep("mut", 5))
#' f <- tempfile()
#' export.broad.cls(cls, f, FALSE, 0)
#' cls2 <- rnorm(10)
#' f2 <- tempfile()
#' export.broad.cls(cls2, f2, TRUE, 0)
#' @export
export.broad.cls <- function(labels, file, continuous=FALSE, start.value=0) {

	if( !continuous ) {
		if( start.value < 0 )
			stop("start.value must be >= 0.\n")
		# if( is.factor(labels) )
		# 	labels <- as.character(labels)
		# line1 <- paste(length(labels), length(unique(labels)), "1")
		# line2 <- paste(c("#",unique(labels)), collapse=" ")
		# line3 <- paste(labels, collapse=" ")
		# writeLines(c(line1,line2,line3), file)
		
		if( is.factor(labels) && !all(unique(labels) == levels(labels)) ) {
			warning("You supplied a factor but the values must correspond with the order of the levels. reordering.\n")
			#
			# if a factor is supplied, the order of the values/labels MUST correspond with the order of the levels
			# line 3 of cls must go 0, 1, 2, 3, ... You can't have 3, 1, 2, 0, 5, ... 
			# Thus the unique labels on line 2 MUST be in the same order as the supplied labels
			#
			labels <- as.character(labels)
		}
		
		#
		# convert characters [back] to a factor, keeping the ordering correct.
		#
		if( !is.factor(labels) )
			labels <- factor(labels, levels=unique(labels))
		levels(labels) <- gsub(" ", ".", levels(labels))
		
		line1 <- paste(length(labels), length(levels(labels)), "1")
		line2 <- paste(c("#",levels(labels)), collapse=" ")
		line3 <- as.numeric(labels) - 1 + start.value
		line3 <- paste(line3, collapse=" ")
		writeLines(c(line1,line2,line3), file)
	}
	else {
		if( is.vector(labels) && !is.list(labels) ) {
			labels <- matrix(labels,nrow=1)
			rownames(labels) <- "Phenotype"
		}
		else if( is.list(labels) && !is.data.frame(labels) ) {
			# elements of list become rows of table.
			labels <- t(list2df(labels))
		}
		else if( is.matrix.like(labels) ) {
			if( is.null(rownames(labels)) )
				rownames(labels) <- paste("Phenotype", 1:nrow(labels), sep="")
		}
		else {
			stop("Unknown format of labels.\n")
		}
		OUT <- file(file, "w")
		writeLines("#numeric", OUT)
		for(i in 1:nrow(labels)) {
			line1 <- paste("#", rownames(labels)[i], sep="")
			line2 <- paste(as.numeric(labels[i,]), collapse=" ")
			writeLines(c(line1, line2), OUT)
		}
		close(OUT)
	}
}
# CHANGELOG
# 2009-06-26
# update 2009-07-08: cls labels can't have spaces in them.
# update 2009-11-20: cls should have numeric values, 0, 1, 2, ...; also made sure that if a factor is supplied that its levels are in the same order as the labels.
