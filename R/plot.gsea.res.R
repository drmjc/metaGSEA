#' Boxplot the intensities for an array split by the detection call (P/M/A).
#' 
#' @param res a gsea res object. see import.gsea.res
#' @param array which array to plot?
#' @return a boxplot
#' @author Mark Cowley, 2010-01-07
#' @export
plot.gsea.res <- function(res, array=1) {
	signal.columns <- grep("Signal$", colnames(res), value=TRUE)
	signal <- res[,signal.columns[array]]
	log <- ""
	if(max(signal, na.rm=TRUE) > 16)
		log <- "y"
	calls <- factor(res[,sub("Signal$", "Call", signal.columns[array])], levels=c("A","M","P"))
	boxplot(signal ~ calls, main=sub(".Signal$", "", signal.columns[array]), ylab="Intensity", log=log)
}
