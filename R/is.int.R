# Checks which elements of a vector are integers.
# This differs from is.integer which only checks the class of the entire vector.
#
# Parameters:
#	x: a vector. usually this will be a numeric vector, but can also be a character
#
# Value:
#	a logical vector, same length as x, with TRUE if each element is an integer
#
# Mark Cowley, 2009-11-30
#


##' Checks which elements of a vector are integers.
##' 
##' This differs from is.integer which only checks the class of the entire
##' vector.
##' 
##' @param x a vector. usually this will be a numeric vector, but can also be a
##'   character
##' @return a logical vector, same length as x, with TRUE if each element is an
##'   integer
##' @author Mark Cowley, 2009-11-30
##' @export
is.int <- function(x) {
	suppressWarnings( as.integer(x) == x )
}
