#' Import a gct file
#'
#' Import a GenePattern gct file into R as a data.frame object. 
#' If a clm.file is also provided, then the columns of the gct are optionally subsetted, reordered, and renamed as per the clm.
#'
#' @param f [character] the path to a gct file
#' @param clm.file [optional] path to a clm.file
#' @param check.names logical: R likes to make sure that the column names are 'syntactically valid', If you don't care, then choose FALSE. default is FALSE.
#' @seealso \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#gct}, \url{http://www.broadinstitute.org/cancer/software/genepattern/tutorial/gp_fileformats.html#clm}
#' @return a \code{data.frame} of usually gene expression data. columns 1&2 are annotation, 3+ are data.
#' @author Mark Cowley, 2009-07-27
#' @export
#' @examples
#' \dontrun{
#' f <- "ftp://ftp.broadinstitute.org/pub/genepattern/datasets/all_aml/all_aml_test.gct"
#' import.gsea.gct(f)
#' }
import.gsea.gct <- function(f, clm.file=NULL, check.names=FALSE) {
	header <- readLines(f,2)[2]
	header <- strsplit(header,"\t")[[1]]
	NROW <- as.numeric(header[1]); NCOL <- as.numeric(header[2])
	res <- read.delim(f, skip=2, check.names=check.names, stringsAsFactors=FALSE, quote="")
	colnames(res)[1:2] <- c("Name", "Description")
	if( max(table(res$Name)) > 1 ) {
		stop("Names in column 1 of your gct file must be unique.\n")
	}
	
	# blank tabs at the bottom don't get skipped by read.table(..., blank.lines.skip=TRUE)
	allNA <- apply(is.na(res), 1, sum) == ncol(res)
	if( any(allNA) ) res <- res[!allNA, ]
	
	# GCT file spec says column 1 must be unique.
	if( any(duplicated(res[,1])) ) stop("Gene or Probe names in column 1 MUST be unique")
	rownames(res) <- res[,1]

	if( !is.null(clm.file) ) {
		clm <- import.gsea.clm(clm.file)
		# clm$CEL <- make.names(clm$CEL) # MJC edited out 2010-07-20
		if( !all(clm$CEL %in% colnames(res)) ) {
			warning("Not all files mentioned in column 1 of the CLM could be found in the gct. skipping.\n")
		}
		else {
			res <- res[,c(colnames(res)[1:2], clm$CEL)]
			colnames(res) <- c(colnames(res)[1:2], clm$Sample)
		}
	}
	if( nrow(res) != NROW ) warning(sprintf("GCT header specified %d rows; found %d.\n", NROW, nrow(res)))
	if( ncol(res)-2 != NCOL ) warning(sprintf("GCT header specified %d samples; found %d.\n", NCOL, ncol(res)-2))
	res
}
# 2010-07-20: changed default for check.names to FALSE
# 2011-09-09: turned quote="" in the read.delim to allow comments with " in them (as per the ftp://ftp.broadinstitute.org/pub/genepattern/datasets/all_aml/all_aml_test.gct file)