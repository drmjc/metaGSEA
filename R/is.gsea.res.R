is.gsea.gct.file <- function(f) {
	grepl("gct$", f, ignore.case=TRUE)
}
is.gsea.res.file <- function(f) {
	grepl("res$", f, ignore.case=TRUE)
}
