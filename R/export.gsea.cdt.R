#' Cluster leading edge genesets, and export in Stanford cdt/gtr format.
#' 
#' cdt files = Clustered Data Table; gtr files = gene tree file. These files
#' describe the ordering of usually a matrix of gene expression data, however
#' we have used this notation to represent the relationships between GSEA
#' leading edge genesets. You can view these files using the
#' \code{HierarchicalClusteringViewer} tool in \code{GenePattern}.
#' 
#' The order of the tree should be identical to \code{\link{plot_gsea.leadingedge.HCL}}, so
#' you could compare a static pdf image (from \code{\link{plot_gsea.leadingedge.HCL}}), or
#' interrogate a dynamically updating image using the \code{GenePattern} tool.
#' 
#' @param x a gsea object (not the leading.edge)
#' @param file the path to a cdt file. This should have the .cdt extension
#' @param gtr if TRUE, then a gtr file will also be made. It gets its filename
#'   from the cdt. if \code{FALSE}, then it will not be made (can't see why that
#'   would be useful to you though?)
#' @param \dots Unprocessed arguments for export.gsea.gtr.
#' @param rename logical: rename the genesets to at most 40 characters? see \code{\link{gsea.rename.genesets}}
#' @author Mark Cowley, 2009-10-12
#' @seealso \code{\link{plot_gsea.leadingedge.HCL}}
#' @references Stanford File formats:
#'   \url{http://smd.stanford.edu/help/formats.shtml} 
#' GenePattern file formats:
#'   \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html}
#' @keywords manip IO file
#' @export
#' @return nothing. Creates a cdt file, and optionally a gtr file.
#' 
export.gsea.cdt <- function(x, file, gtr=TRUE, rename=FALSE) {
	cdt.file <- file
	if(gtr)
		gtr.file <- sub("cdt", "gtr", cdt.file)
	
	d <- gsea.leadingedge.distance(x$leading.edge)
	hc <- hclust(d)
	
	#
	# using the examples shown at Stanford:
	# http://smd.stanford.edu/help/formats.shtml#cdt
	# generate a pcl object, then a cdt and gtr.

	#
	# column 1 = geneset name, but column 2 = description.
	# Use the gsea.rename.genesets to append rank, fdr and direction to the geneset name.
	# (I exclude the edb object to make the code run quite a bit faster.)
	names <- x$tt$NAME
	if( rename ) {
		MAX.LABEL.LENGTH <- 40
		tmp <- gsea.rename.genesets(x[setdiff(names(x), "edb")], direction=TRUE, fdr=TRUE, rank=TRUE, maxlen=MAX.LABEL.LENGTH)
		names <- tmp$tt$NAME[match(hc$labels, x$tt$NAME)]
	}

	pcl <- cbind(
		ID=hc$labels, 
		NAME=names, 
		GWEIGHT=1, 
		NES=x$tt$NES[match(hc$labels, x$tt[,1])]
	)

	ids <- paste("GENESET", 1:nrow(pcl), "X", sep="")
	cdt <- cbind(GID=ids, pcl)
	cdt <- cdt[hc$order,]
	write.delim(cdt, cdt.file, row.names=FALSE, col.names=TRUE)
	
	if( gtr ) {
		tmp <- colapply(hc$merge, function(x) {
			idx <- x < 0
			res <- rep("", length(x))
			res[idx] <- sprintf("GENESET%dX",-x[idx])
			res[!idx] <- sprintf("NODE%dX",x[!idx])
			res
		})
		gtr <- data.frame(
			paste("NODE", 1:nrow(hc$merge), "X", sep=""),
			tmp,
			round(1-hc$height, 6),
			stringsAsFactors=FALSE
		)
	
		write.delim(gtr, gtr.file, row.names=FALSE, col.names=FALSE)
	}
}

#' You should use export.gsea.cdt instead.
#' @param \dots unused
#' @author Mark Cowley, 2009-10-12
#' @export
export.gsea.gtr <- function(...) {
	stop("You should use export.gsea.cdt.\n")
}

# Here's some example output from an hclust object....
#
# > hc$merge
#       [,1] [,2]
#  [1,]  -22  -23
#  [2,]   -4   -9
#  [3,]  -21    1
#  [4,]  -12  -26
#  [5,]   -3  -11
#  [6,]   -2  -19
#  [7,]   -7  -15
#  [8,]  -10  -14
#  [9,]   -1   -5
# [10,]  -18    7
# [11,]  -13  -25
# [12,]  -27    8
# [13,]  -16    5
# [14,]    6    9
# [15,]    3    4
# [16,]    2   14
# [17,]   -6  -29
# [18,]  -28   12
# [19,]   13   16
# [20,]  -17   15
# [21,]   18   19
# [22,]   10   11
# [23,]   17   22
# [24,]   20   21
# [25,]   23   24
# [26,]   -8   25
# [27,]  -20   26
# [28,]  -24   27
# [29,]  -30   28
# > hc$height
#  [1] 0.231 0.387 0.431 0.630 0.714 0.727 0.757 0.779 0.808 0.818 0.820 0.842
# [13] 0.842 0.855 0.860 0.905 0.917 0.932 0.940 0.943 0.961 0.961 0.979 0.981
# [25] 1.000 1.000 1.000 1.000 1.000
# > hc$order
#  [1] 30 24 20  8  6 29 18  7 15 13 25 17 21 22 23 12 26 28 27 10 14 16  3 11  4
# [26]  9  2 19  1  5
#
