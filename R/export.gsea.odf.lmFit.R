#' Export a limma lmFit to a GSEA ODF file.
#' From two linear model files (see details), the coefficient of interest, export an
#' ODF file, which can be viewed in the GenePattern ComparativeMarkerSelectionViewer 
#' tool. The ODF file that is created is heavily modelled on the ODF file that the
#' ComparativeMarkerSelection GenePattern module creates.
#' 
#' This is a key function for the LimmaGP GenePattern module. 
#' In LimmaGP, we always perform 
#' two linear model fits, the first is a treatment means 
#' parameterisation (\code{fit1}, hereafter), 
#' the second is a contrast fit after constructing the appropriate contrasts (\code{fit2}, hereafter).
#' From \code{fit1}, we can get the group means, and standard errors, and from \code{fit2}, 
#' we can get the logFC,
#' moderated t-statistics, pvalues, and FDR/FWER values.
#'
#' @param fit1 a linear model fit, \emph{before} fitting contrasts, and running eBayes
#' @param fit2 a linear model fit, \emph{after} fitting contrasts, and running eBayes
#' @param coef the coefficient to be exported. see \code{\link[limma]{topTable}}
#' @param file the path to the ODF file
#' @param gct.file the name of the gct file that the linear models were fitted to. 
#'   This is written to the ODF header
#' @param cls.file the name of the cls file that the linear models were fitted to. 
#'   This is written to the ODF header
#' @param description a character vector, which will end up in the description column.
#'    if \code{NULL}, then the probe ID's will be used.
#' @param collapse logical: whether to collapse probes to genes? if \code{TRUE}, then 
#'     \code{probe2gene} must be set.
#' @param probe2gene if \code{collapse=TRUE}, then a 2 column data.frame mapping from
#'    probe to gene symbol should be supplied. default: \code{NULL}.
#' @return nothing. writes an ODF file.
#' @export
#' @author Mark Cowley, 2011-07-19
#'
export.gsea.odf.lmFit <- function(fit1, fit2, coef=1, file,
	gct.file, cls.file, description, collapse=FALSE, probe2gene=NULL) {

	if( ! "contrasts" %in% names(fit2) ) {
		stop("The fit2 must have used contrasts (so I can work out class 0 and class 1 means and names).\n")
	}
	if ( collapse && is.null(probe2gene) ) {
		collapse <- FALSE
		cat("you must supple a probe2gene table if settings collapse=TRUE.\n")
	}
	
	require(qvalue)

	tt <- topTable(fit2, coef=coef, number=nrow(fit2), adjust.method="none")
	tt <- tt[order(abs(tt$t), decreasing=TRUE), ]
	if( all(abs(tt$t)<0.000001) )
		tt$t <- tt$logFC
		
	tt$Rank <- rank( (tt$t-max(tt$t, na.rm=TRUE))*-1, ties.method="first" )
	tt$Feature <- tt$ID
	# tt$Description <- ""
	tt$Score <- as.double(tt$t)
	tt$"Fold Change" <- as.double(2^tt$logFC)
	tt$"log2 Fold Change" <- as.double(tt$logFC)
	tt$"Feature P" <- as.double(tt$P.Value)
	
	#
	# adjust for multiple testing
	#
	idx <- !is.na(tt$"Feature P")
	tt$"FDR(BH)"[idx] <-  as.double(p.adjust(tt$"Feature P"[idx], "BH", n=nrow(tt)))
	tt$Bonferroni[idx] <- as.double(p.adjust(tt$"Feature P"[idx], "bonferroni", n=nrow(tt)))

	q <- qvalue2(tt$"Feature P"[idx], lambda=seq(0,0.95,0.05)) # defined in pwbc, as a safer alternative to qvalue
	qdata <- .qplot.data(q, rng = c(0, 0.1), smooth.df = 3, smooth.log.pi0 = FALSE)
		
	tt$"Q Value"[idx] <-  as.double(q$qvalues)
	if(any(is.na(tt[,c("Score", "Feature P", "FDR(BH)", "Q Value", "Bonferroni")])))
		warning("There are NA's in the P-value column(s)!")

	#
	# determine a description, either from the supplied description argument, or the chip file.
	#
	if( is.null(description) ) {
		description <- rep("", nrow(tt))
	}
	else {
		description <- description[tt$Feature]
	}
	tt$Description <- description

	#
	# collapse probes to gene symbols AFTER doing multiple testing.
	# why? we had to do the t-stats in the first place in order to determine the most DE genes,
	# so you have to mult correct for those tests.
	#
	if( collapse ) {
		if ( !is.null(probe2gene) && is.matrix.like(probe2gene) && (ncol(probe2gene) >= 2) ) {
			colnames(probe2gene) <- c("ID", "Gene.Symbol")
			best.probe <- calc.best.probe.topTable(tt, probe2gene, toupper=FALSE, only.genes=FALSE)
			tt <- merge(tt, best.probe, by.x="ID", by.y="ID", all.x=FALSE, all.y=TRUE, sort=FALSE)
			tt <- move.column(tt, "Gene.Symbol", 0)
			colnames(tt)[colnames(tt) == "ID"] <- "ProbeID"
			tt$Feature <- tt$Gene.Symbol
		}
		else {
			collapse <- FALSE
		}
	}
	tt$Gene.Symbol <- NULL
	# if( !collapse ) {
	# 	tt$Gene.Symbol <- rep("", nrow(tt))
	# }

	#
	# get the class 0 and class 1 mean and std
	# (might want to consider subset.MArrayLM)
	# NB: this comes from fit1, not fit2 (since contrasts would have already been fit)
	#
	if( all(fit2$contrasts[, coef] %in% c(0,1)) ) {
		# then there was no comparison of classes -- probably because the model matrix was a unity vector, and thus the contrast matrix is the same value as the coefficient (ie no comparison of 2 classes was made.)
		class0 <- rownames(fit2$contrasts)[which(fit2$contrasts[, coef] == 1)]
		class1 <- "Baseline"
		
		if( ! all(c(class0) %in% colnames(fit1$coefficients)) )
			stop("The names of the coefficients do not match the names of the contrasts.\n")

		o <- match(tt$ID, fit2$genes[,1])
		tt$"class0 Mean" <- fit1$coefficients[o, class0]
		tt$"class1 Mean" <- rep(0, nrow(tt))

		## these are the unmoderated SE's
		# errors <- fit1$stdev.unscaled * fit1$sigma
		# 2011-06-20: how did this ever work? s2.post is added by eBayes & thus only in fit2
		# ## these are the moderated SE's update: 2010-11-12
		# errors <- fit1$stdev.unscaled * sqrt(fit1$s2.post)
		errors <- fit1$stdev.unscaled * sqrt(fit2$s2.post)
		
		if( any(is.na(errors)) ) { # this can happen if there was a 1vs1 analysis for instance.
			errors[is.na(errors)] <- 0
		}
		tt$"class0 Std"  <- errors[o, class0]
		tt$"class1 Std" <- rep(0, nrow(tt))
	}
	else if( all(fit2$contrasts[, coef] %in% c(-1,0,1)) ) {
		class0 <- rownames(fit2$contrasts)[which(fit2$contrasts[, coef] == 1)]
		class1 <- rownames(fit2$contrasts)[which(fit2$contrasts[, coef] == -1)]

		if( ! all(c(class0, class1) %in% colnames(fit1$coefficients)) )
			stop("The names of the coefficients do not match the names of the contrasts.\n")

		o <- match(tt$ID, fit2$genes[,1])
		tt$"class0 Mean" <- fit1$coefficients[o, class0]
		tt$"class1 Mean" <- fit1$coefficients[o, class1]

		## these are the unmoderated SE's
		# errors <- fit1$stdev.unscaled * fit1$sigma
		# 2011-06-20: how did this ever work? s2.post is added by eBayes & thus only in fit2
		# ## these are the moderated SE's update: 2010-11-12
		# errors <- fit1$stdev.unscaled * sqrt(fit1$s2.post)
		errors <- fit1$stdev.unscaled * sqrt(fit2$s2.post)

		if( any(is.na(errors)) ) { # this can happen if there was a 1vs1 analysis for instance.
			errors[is.na(errors)] <- 0
		}
		tt$"class0 Std"  <- errors[o, class0]
		tt$"class1 Std"  <- errors[o, class1]
	}
	else {
		# custom contrasts have been used
		#
		# class 0 is some sort of average of the positive cofficients, 
		# class 1 is some sort of average of the negative coefficients... 
		# ... do a seperate fit.contrasts for just the positive, or just 
		# the negative values to get the class 0 and class 1 means and errors.
		#
		class0 <- "class0"
		class1 <- "class1"
		# can we improve on these names?
		# the GenePattern each vs rest module has a single 1.0 coef, and the rest are all negative and equal
		x <- fit2$contrasts[, coef]
		a <- which(x==1)
		if( length(a)==1 && all(x[-a]<0) && alleq(x[-a]) ) {
			class0 <- rownames(fit2$contrasts)[a]
			class1 <- paste("non-", class0, sep="")
		}

		pos.contrasts <- neg.contrasts <- fit2$contrasts
		pos.contrasts[pos.contrasts<0] <- 0
		neg.contrasts[neg.contrasts>0] <- 0
		neg.contrasts <- neg.contrasts * -1 # make everything positive
		fit.pos <- contrasts.fit(fit1, pos.contrasts)
		fit.neg <- contrasts.fit(fit1, neg.contrasts)
		errors.pos <- fit.pos$stdev.unscaled * fit.pos$sigma
		errors.neg <- fit.neg$stdev.unscaled * fit.neg$sigma
		o <- match(tt$ID, fit1$genes[,1])
		tt$"class0 Mean" <- fit.pos$coefficients[o, coef]
		tt$"class1 Mean" <- fit.neg$coefficients[o, coef]
		tt$"class0 Std"  <- errors.pos[o, coef]
		tt$"class1 Std"  <- errors.neg[o, coef]
	}
	tt$"class0 Mean" <- as.double(tt$"class0 Mean")
	tt$"class1 Mean" <- as.double(tt$"class1 Mean")
	tt$"class0 Std"  <- as.double(tt$"class0 Std" )
	tt$"class1 Std"  <- as.double(tt$"class1 Std" )

	#
	# select and order the appropriate columns.
	#
	res <- tt[, c("Rank", "Feature", "Description", "Score", "Feature P", "FDR(BH)", "Q Value", "Bonferroni", "Fold Change", "log2 Fold Change", "class0 Mean", "class0 Std", "class1 Mean", "class1 Std")]
	colnames(res) <- sub("class0", class0, colnames(res))
	colnames(res) <- sub("class1", class1, colnames(res))
	
	headerLines <-  c(
		"Model=Comparative Marker Selection",
		paste("Dataset File=", basename(gct.file), sep=""),
		paste("Class File=", basename(cls.file), sep=""),
		"Test Direction=2 Sided",
		paste("Class 0=", class0, sep=""),
		paste("Class 1=", class1, sep=""),
		"Test Statistic=TRUE-Test",
		paste("pi0", round(qdata$pi0.estimate, 7), sep="="),
		paste("lambda",paste(qdata$lambda, collapse=" "), sep="="),
		paste("pi0(lambda)",paste(round(qdata$pi0, 7), collapse=" "), sep="="),
		paste("cubic spline(lambda)",paste(round(qdata$spline, 7), collapse=" "), sep="=")
	)
	
	# are any values Infinite?
	idx <- is.infinite(res)
	if( sum(idx) > 0 ) {
		cat("Warning, infinite values were created. Did you use unlogged data with Limma?\n")
		res[idx] <- 2^31
	}

	classes <- c("int", "String", "String", rep("double", ncol(res)-3))
	export.gsea.odf(res, file, headerLines, colclasses=classes)
}
# CHANGELOG
# 2010-11-12: realised that the errors I was producing were the unmoderated SE's. Fixed to
# moderated SE's. So if moderated.T = logFC/moderated.SE, accessing the stdev.unscaled & sqrt(s2.post) gets you what you need.
# semi-confirmed here: http://permalink.gmane.org/gmane.science.biology.informatics.conductor/31688
# 2011-02-22: added the ability to collapse probes 2 genes via collapse & probe2gene arguments.
# 2011-07-06: dropped the empty Gene.Symbol column.
# 2011-10-25: now uses qvalue2, which prevents the errors due to pi_0 estimation fail
# 2011-11-03: try to improve on class0, class1 labels when using each-vs-rest in GenePattern.
# - for 3 groups, the contrasts look like this: [[c(1,-0.5,-0.5), c(-0.5,1,-0.5), c(-0.5,-0.5,1)]]
# ... ie always a '1' and the rest are all negative and identical.
# 

# export.gsea.odf.topTable <- function(tt, fit, coef,
# 	file,
# 	gct.file, cls.file, chip.file, 
# 	class0, class1) {
# 	
# 	if( ! "contrasts" %in% names(fit) ) {
# 		stop("The fit must have used contrasts (so I can work out class 0 and class 1 means and names).\n")
# 	}
# 	
# 	tt$Rank <- rank( abs(tt$t) )
# 	tt$Feature <- tt$ID
# 	# tt$Description <- ""
# 	tt$Score <- tt$t
# 	tt$"Feature P" <- tt$P.Value
# 	# tt$"FDR(BH)" <- 1.0
# 	# tt$"Q Value" <- 1.0
# 	# tt$Bonferroni <- 1.0
# 	# tt$"class0 Mean" <- 0
# 	# tt$"class0 Std" <- 0
# 	# tt$"class1 Mean" <- 0
# 	# tt$"class1 Std" <- 0
# 	
# 	#
# 	# adjust for multiple testing
# 	#
# 	idx <- !is.na(tt$"Feature P")
# 	tt$"FDR(BH)"[idx] <-  p.adjust(tt$"Feature P"[idx], "BH", n=nrow(tt))
# 	tt$"Q Value"[idx] <-  p.adjust(tt$"Feature P"[idx], "Q", n=nrow(tt))
# 	tt$Bonferroni[idx] <- p.adjust(tt$"Feature P"[idx], "bonferroni", n=nrow(tt))
# 	
# 	#
# 	# import the chip, and add the annotations.
# 	#
# 	chip <- import.gsea.chip(chip.file)
# 	tt$Description <- chip$Description[match(tt$Feature, chip[,1])]
# 	
# 	#
# 	# get the class 0 and class 1 mean and std
# 	#
# 	# tt$"class0 Mean" <- 0
# 	# tt$"class0 Std" <- 0
# 	# tt$"class1 Mean" <- 0
# 	# tt$"class1 Std" <- 0
# 
# 
# 	#
# 	# select and order the appropriate columns.
# 	#
# 	res <- tt[, c("Rank", "Feature", "Description", "Score", "Feature P", "FDR(BH)", "Q Value", "Bonferroni", "Fold Change", "class0 Mean", "class0 Std", "class1 Mean", "class1 Std")]
# 	colnames(res) <- sub("class0", class0, colnames(res))
# 	colnames(res) <- sub("class1", class1, colnames(res))
# 	
# 	headerLines <-  c(
# 		"Model=Comparative Marker Selection",
# 		paste("Dataset File=", basename(gct.file), sep=""),
# 		paste("Class File=", basename(cls.file), sep=""),
# 		"Test Direction=2 Sided",
# 		paste("Class 0=", class0, sep=""),
# 		paste("Class 1=", class1, sep=""),
# 		"Test Statistic=TRUE-Test"
# 		)
# 	
# 	export.gsea.odf(res, file, headerLines)
# }