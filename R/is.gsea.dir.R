is.gsea.dir <- function(x) {
	!is.null(x) & !is.na(x) & file.exists(x) & is.dir(x) & file.exists(file.path(x, "index.html"))
}
# is.gsea.dir.dir <- function(x) {
# 	if( file.exists(x) && is.dir(x) && !is.gsea.dir(x) ) {
# 		dirs <- dir(x, full.path=TRUE)
# 		dirs <- dirs[is.dir(dirs)]
# 		any(sapply(dirs, is.gsea.dir))
# 	}
# 	FALSE
# }
# 