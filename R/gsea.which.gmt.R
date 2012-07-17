#' What gmt file was used?
#' 
#' What gmt file was used? This parses the gmx file name from within the rpt
#' file
#' 
#' @param dir the dir that contains index.html
#' @return eg "c2_cgp"
#' @author Mark Cowley, 2009-04-07
#' @export
gsea.which.gmt <- function(dir) {
	rpt <- import.gsea.rpt(dir)
	tmp <- rpt$gmx
	# f <- dir(dir, pattern=".*\\.rpt", full.names=TRUE)
	# tmp <- readLines(f)
	# tmp <- tmp[grep("param.gmx", tmp)]
	# tmp <- strsplit(tmp, "\t")[[1]][3]
	res <- basename(tmp)
	res <- strsplit(res, "\\.")[[1]]
	res <- paste(res[1:2], collapse="_")
	res
}
