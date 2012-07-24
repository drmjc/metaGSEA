#' GenePattern to SPD
#'
#' Convert from GenePattern GCT and CLS files to the input file required for Sample Progression Discovery (SPD)
#' The input file is a matlab binary file, containing probe names, data matrix and samples names, and optionally
#' class labels, from a cls file.
#' This makes use of the R.matlab library by Henrik Bengtsson & creates v5 matlab files (which work fine with matlab v7.)
#'
#' @param gct.file The path to a GenePattern formatted GCT file
#' @param cls.files The path to a GenePattern formatted CLS file. Optional.
#' @param cls.names What do the classes in cls.file represent? Tumour Type? Cell Line? Specify a short label, eg \dQuote{Tumour Type}
#' @param out.file The path to the output file. It should end in .mat.
#' @return Creates a version 5 matlab file which you can then open in Matlab
#' @author Mark Cowley, 2011-04-20
#' @section TODO: support >1 cls files
#' @seealso 
#' \url{http://icbp.stanford.edu/software/SPD},
#' \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#gct},
#' \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#cls}
#' @examples
#' \dontrun{
#' setwd("/Volumes/marcow-1/projects/SPD/analysisV1")
#' genepattern2spd("norm.rsn.medianNorm.gct", out.file="norm.rsn.medianNorm.SPD.mat")
#' genepattern2spd("norm.rsn.medianNorm.gct", "norm.rsn.medianNorm.cls", "Tumour Grade", "norm.rsn.medianNorm.SPD.mat")
#' }
#' @export
genepattern2spd <- function(gct.file, cls.files, cls.names="Class", out.file="spd.mat") {
	require(R.matlab) || stop("required package 'R.matlab' is not installed")
	!missing(gct.file) || stop("gct.file MUST be specified")
	!missing(out.file) || stop("out.file MUST be specified")
	
	gct <- import.gsea.gct(gct.file)

	probe_names <- matrix(as.character(gct[,1]))
	exp_names <- colnames(gct[3:ncol(gct)])
	dat <- as.matrix(gct[,3:ncol(gct)])
	
	if( !missing(cls.files) ) {
		!missing(cls.names) || stop ("cls.names required if your supply a cls.file. It should be just a descriptive name.")
		# @TODO support >1 cls files
		length(cls.files) == 1 || stop("Can't handle >1 cls.file at the moment.")
		cls.files <- cls.files[1]
		cls.names <- cls.names[1]
		
		cls <- import.gsea.cls(cls.files)
		length(cls) == ncol(dat) || stop("mismatch between num samples in GCT and num samples in CLS.")
		
		color_code_vectors <- t(matrix(as.numeric(cls)-1))
		R.matlab::writeMat(out.file, probe_names=probe_names, exp_names=exp_names, data=dat, color_code_names=cls.names, color_code_vectors=color_code_vectors)
	}
	else {
		R.matlab::writeMat(out.file, probe_names=probe_names, exp_names=exp_names, data=dat)
	}
}
