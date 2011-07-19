# Combination GSEA Leading Edge plot.
#
# NB, it's only worth using on filtered GSEA data - see gsea.filter(<>, N=50, direction="up")
# See also plot.gsea.leadingedge.HCL, plot.gsea.leadingedge.heatmap, plot.gsea.leadingedge.barplot
#
# Parameters:
#	x: a GSEA list
#	main: the main plotting title
#
# Value:
#	creates multiple plots.
#
# Mark Cowley, 2009-10-29


##' Combination GSEA Leading Edge plot.
##' 
##' NB, it's only worth using on filtered GSEA data - see gsea.filter(<>, N=50,
##' direction="up")
##' See also plot.gsea.leadingedge.HCL, plot.gsea.leadingedge.heatmap,
##' plot.gsea.leadingedge.barplot
##' 
##' @param x a GSEA list
##' @param main the main plotting title
##' @return creates multiple plots.
##' @author Mark Cowley, 2009-10-29
##' @export
plot.gsea.leadingedge <- function(x, main) {
	plot.gsea.leadingedge.HCL(x, main=main)
	plot.gsea.leadingedge.heatmap(x, main=main)
	plot.gsea.leadingedge.barplot(x, main=main, min.count=2)
	plot.gsea.leadingedge.adjmat(x, main=main)
}
