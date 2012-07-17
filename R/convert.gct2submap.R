# Prepare data for use in SubMap.
#
# Converts a gct file and cls file so that:
# - one row per gene symbol & gene symbol becomes the Name of each row, rather than the probe id
# -- different methods of collapsing rows
# - cls file starts from 1, rather than 0
# - choose the most variable genes, either by a coefficient of variation threshold (eg 50)
#   or the top N most variable genes, as defined by CV.
#
# Parameters:
#	gct.file, cls.file: standard input files to be modified.
#	chip.file: used to map probe ID's to gene symbols
#	gct.out, cls.out: output files.
#	collapse.mode: how to choose the value for each gene.
#		max.avg: calculate the average for each probeset & choose the values from the probes with largest average value
#		mean, median: for N probesets mapping to the same gene, calculate the mean or median of those N values for each sample.
#		max.cv: calculate the CV for each probe & choose the values from the probe with highest cv
#	cv.thresh: [optional] CV threshold to use. CV's are percentages, so 50, or 100 is a good place to start. Thresholding is done after collapsing to unique genes. This overrides topN.thresh.
#	topN.thresh: [optional] take the top N genes when ranked by their CV's. Must leave cv.thresh=NULL. Thresholding is done after collapsing to unique genes.s
#
# Warning:
#	If the cv.thresh is too strict, an error is thrown
#	The cls file starts from '1' instead of '0' which means it may not work with other GenePattern modules!
#
# Value:
#	creates a new gct & cls file. The gct file has 1 row per gene which passed the filters & has gene symbols in column 1 (the Name column). The cls file starts from '1' instead of '0' which means it may not work with other GenePattern modules!
#
# Mark Cowley, 2010-07-06


##' Prepare data for use in SubMap.
##' 
##' Converts a gct file and cls file so that:
##' - one row per gene symbol & gene symbol becomes the Name of each row,
##' rather than the probe id
##' -- different methods of collapsing rows
##' - cls file starts from 1, rather than 0
##' - choose the most variable genes, either by a coefficient of variation
##' threshold (eg 50)
##' or the top N most variable genes, as defined by CV.
##' 
##' @param gct.file standard input files to be modified.
##' @param cls.file standard input files to be modified.
##' @param chip.file used to map probe ID's to gene symbols
##' @param gct.out output files.
##' @param cls.out output files.
##' @param collapse.mode how to choose the value for each gene.
##' @param max.avg calculate the average for each probeset & choose the values
##'   from the probes with largest average value
##' @param mean for N probesets mapping to the same gene, calculate the mean or
##'   median of those N values for each sample.
##' @param median for N probesets mapping to the same gene, calculate the mean
##'   or median of those N values for each sample.
##' @param max.cv calculate the CV for each probe & choose the values from the
##'   probe with highest cv
##' @param cv.thresh [optional] CV threshold to use. CV's are percentages, so
##'   50, or 100 is a good place to start. Thresholding is done after
##'   collapsing to unique genes. This overrides topN.thresh.
##' @param topN.thresh [optional] take the top N genes when ranked by their
##'   CV's. Must leave cv.thresh=NULL. Thresholding is done after collapsing to
##'   unique genes.s
##' @return creates a new gct & cls file. The gct file has 1 row per gene which
##'   passed the filters & has gene symbols in column 1 (the Name column). The
##'   cls file starts from '1' instead of '0' which means it may not work with
##'   other GenePattern modules!
##' @author Mark Cowley, 2010-07-06
##' @export
convert.gct2submap <- function(
	gct.file,
	cls.file,
	chip.file,
	gct.out,
	cls.out,
	collapse.mode=c("max.avg", "mean", "median", "max.cv"),
	cv.thresh=NULL,
	topN.thresh=NULL
	) {

	collapse.mode <- collapse.mode[1]
	
	# stopifnot(all(file.exists(gct.file, cls.file)))
	stopifnot(
		!missing(gct.file), !missing(cls.file), !missing(chip.file),
		file.exists(gct.file), file.exists(cls.file), file.exists(chip.file),
		!missing(gct.out), !missing(cls.out)
	)
	
	gct <- import.gsea.gct(gct.file, check.names=FALSE)
	chip <- import.gsea.chip(chip.file)
	chip <- chip[match(gct$Name, chip$Probe.Set.ID), ]
	chip$Gene.Symbol <- sub(" /// .*", "", chip$Gene.Symbol)
	# trim out the rows with no Gene.Symbol
	idx <- !is.na(chip$Gene.Symbol ) | nchar(chip$Gene.Symbol)>0
	gct <- gct[idx, ]
	chip <- chip[idx, ]

	gct$Name <- chip$Gene.Symbol
	gct$Description <- chip$Gene.Title
	
	#
	# Now we need to make gct$Name unique (as per file fmt definition)
	# -- ie need to collapse multiple rows for the same gene symbol -> 1 row
	#
	
	# wrapper which can run numeric functions, like mean, median, var, cv... on numeric
	# and character data. characters are collapsed into a single comma sep string
	.safe.FUN <- function(x, FUN) {
		if(is.character(x)) {
			if( alleq(x) )
				x[1]
			else
				paste(unique(x), collapse=", ")
		}
		else
			FUN(x)
	}

	# calculate the coefficient of variance.
	cv <- function(x, na.rm=TRUE) { 
		mean(x, na.rm=na.rm)/sd(x, na.rm=na.rm)
	}

	# collapse the rows of the gct file
	switch(collapse.mode,
		mean={
			gct <- aggregate(gct, list(Gene=gct$Name), .safe.FUN, mean)
			gct$Name <- gct$Gene
			gct$Gene <- NULL
		},
		median={
			gct <- aggregate(gct, list(Gene=gct$Name), .safe.FUN, median)
			gct$Name <- gct$Gene
			gct$Gene <- NULL
		},
		max.avg={
			# reorder gct from highest avg expression to lowest
			o <- order(apply(gct[,3:ncol(gct)], 1, mean, na.rm=TRUE), decreasing=TRUE)
			gct <- gct[o,]
			# then select the first instance of each gene symbol using match,
			# which returns only the first hit.
			idx <- match(sort(unique(gct$Name)), gct$Name)
			gct <- gct[idx, ]
		},
		max.cv={
			# reorder gct from highest variability to lowest
			o <- order(apply(gct[,3:ncol(gct)], 1, cv, na.rm=TRUE), decreasing=TRUE)
			gct <- gct[o,]
			# then select the first instance of each gene symbol using match,
			# which returns only the first hit.
			idx <- match(sort(unique(gct$Name)), gct$Name)
			gct <- gct[idx, ]
		}
	)
	
	#
	# now, gct has one row per unique gct$Name
	# - resort it by gene symbol
	#
	gct <- gct[order(gct$Name), ]
	rownames(gct) <- gct$Name
	
	#
	# threshold by cv or top N genes
	#
	
	row.cv <- apply(gct[,3:ncol(gct)], 1, cv)
	if( !is.null(cv.thresh) ) {
		idx <- which(row.cv > cv.thresh)
		if( length(idx)>0 ) {
			cat( sprintf("Filtering out %d genes, with cv < %.1f, leaving %d genes.\n",
				length(row.cv) - length(idx),
				cv.thresh,
				length(idx)
			) )
			gct <- gct[idx, ]
		}
		else {
			stop(sprintf("CV threshold is too strict, all genes filtered out. The 90th percentile of CV is: %.2f.\n", quantile(row.cv, 0.9)) )
		}
	}
	else if( !is.null(topN.thresh) && topN.thresh < nrow(gct) ) {
		cat(sprintf("Filtering out genes outside of the top %d.\n", topN.thresh))
		o <- order(row.cv, decreasing=TRUE)
		gct <- gct[o,]
		gct <- gct[1:topN.thresh, ]
	}
	
	export.gsea.gct(gct, file=gct.out)

	# ensure the cls file starts at '1'
	if( file.exists(cls.file) ) {
		cls <- import.gsea.cls(cls.file)
		# cls
		export.broad.cls(cls, cls.out, continuous=FALSE, start.value=1)
	}
}
