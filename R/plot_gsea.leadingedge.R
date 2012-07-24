#' Combination GSEA Leading Edge plot.
#' 
#' @note it's only worth using on filtered GSEA data - 
#' see \code{gsea.filter(*, N=50, direction="up")}
#' 
#' @param x a GSEA list
#' @param main the main plotting title
#' 
#' @return none. creates multiple plots.
#' 
#' @author Mark Cowley, 2009-10-29
#' 
#' @seealso \code{\link{gsea.filter}} \code{\link{plot_gsea.leadingedge.HCL}} \code{\link{plot_gsea.leadingedge.heatmap}}
#' \code{\link{plot_gsea.leadingedge.barplot}}
#' @export
plot_gsea.leadingedge <- function(x, main) {
	plot_gsea.leadingedge.HCL(x, main=main)
	plot_gsea.leadingedge.heatmap(x, main=main)
	plot_gsea.leadingedge.barplot(x, main=main, min.count=2)
	plot_gsea.leadingedge.adjmat(x, main=main)
}
