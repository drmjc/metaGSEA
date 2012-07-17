#' Convert a res object to a gct object.
#' 
#' @param res a res object. see import.gsea.res
#' @return a gct object.
#' @author Mark Cowley, 2009-12-18
#' @export
convert.gsea.res2gct <- function(res) {
	idx <- grep("\\.Signal", colnames(res))
	gct <- res[,c(2,1,idx)]
	colnames(gct) <- sub("\\.Signal", "", colnames(gct))
	colnames(gct)[1:2] <- c("Name", "Description")
	
	return(gct)
}
