#
# The Jacquard coefficient as a measure of overlap in terms of genes in 2 gene sets.
# Motivated by the leading edge analysis of Broad's GSEA.
# 
# Parameters:
# 	a,b: a vector of entities, be they characters or integers
#
# Value:
# 	a decimal value in [0,1] where 0 represents no overlap, and 1.0 indicating 100% overlap - ie a=b
#
# Examples:
#	jacquard(letters[1:5], letters[3:5])
# 	# [1] 0.6
#
# Mark Cowley, 2009-04-06
# 


##' The Jacquard coefficient as a measure of overlap in terms of genes in 2
##' gene sets.
##' 
##' Motivated by the leading edge analysis of Broad's GSEA.
##' 
##' @param a a vector of entities, be they characters or integers
##' @param b a vector of entities, be they characters or integers
##' @return a decimal value in [0,1] where 0 represents no overlap, and 1.0
##'   indicating 100% overlap - ie a=b
##' @author Mark Cowley, 2009-04-06
##' @examples
##' jacquard(letters[1:5], letters[3:5])
##' # [1] 0.6
##' @export
jacquard <- function(a,b) {
	length(intersect(a,b)) / length(union(a,b))
}
