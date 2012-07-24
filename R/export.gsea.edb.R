#' Export a GSEA edb object to an \code{XML} file.
#' 
#' @param edb an edb object.
#' @param file the name of an xml file. Usually this is inside a dir called
#'   edb, and ends in \dQuote{.edb}
#' @return creates an XML file
#' @author Mark Cowley, 2009-10-12
#' @seealso \code{\link{import.gsea.edb}}, \code{\link[XML]{saveXML}}
#' @keywords file IO
#' @export
#' @importFrom XML saveXML
export.gsea.edb <- function(edb, file) {
	saveXML(edb, file, encoding="UTF-8")
	perl.oneliner("perl -pi -e 's|<DTG|\n  <DTG|g'", file)
}
