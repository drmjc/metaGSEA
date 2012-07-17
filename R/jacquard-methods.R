#' The Jacquard coefficient as a measure of overlap in terms of genes in 2
#' gene sets.
#' 
#' Motivated by the leading edge analysis of Broad's GSEA.
#' 
#' @param a a vector of entities, be they \code{character}s, \code{integer}s or \code{GeneSet}
#' @param b a vector of entities, be they \code{character}s, \code{integer}s or \code{GeneSet}
#' @return a decimal value in [0,1] where 0 represents no overlap, and 1.0
#'   indicating 100% overlap - ie a=b
#' @author Mark Cowley, 2009-04-06
#' 
#' @exportMethod jacquard
#' @importClassesFrom GSEABase GeneSet
#' @rdname jacquard-methods
#' @docType methods
#' @examples
#' jacquard(letters[1:5], letters[3:5])
#' jacquard(1:5, 3:5)
setGeneric(
	"jacquard",
	function(a, b) {
		standardGeneric("jacquard")
	}
)


#' @rdname jacquard-methods
#' @aliases jacquard,character,character-method
setMethod(
	"jacquard",
	signature=signature("character", "character"),
	function(a, b) {
		length(intersect(a, b)) / length(union(a, b))
	}
)

#' @rdname jacquard-methods
#' @aliases jacquard,numeric,numeric-method
setMethod(
	"jacquard",
	signature=signature("numeric", "numeric"),
	function(a, b) {
		jacquard(as.character(a), as.character(b))
	}
)

#' @rdname jacquard-methods
#' @aliases jacquard,GeneSet,GeneSet-method
setMethod(
	"jacquard",
	signature=signature("GeneSet", "GeneSet"),
	function(a, b) {
		jacquard(geneIds(a), geneIds(b))
	}
)

# Defunct: non-S4 method
# jacquard <- function(a,b) {
# 	length(intersect(a,b)) / length(union(a,b))
# }
