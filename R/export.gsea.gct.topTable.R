#' Export a GCT file, from a GSEA topTable plus data table. 
#' This is a companion function to \code{\link{export.gsea.rnk.topTable}} which makes a gct
#' file which can then be used for annotating GSEA output with heatmaps.
#' 
#' @param tt a GSEA top-table at the probe-level
#' @param rma a data.frame of expression data with rownames being the
#'   probesetID
#' @param chip the chip file that will be used in the gseaPreRanked
#'   analysis. Either the path to the file name, or the imported version as a
#'   data.frame. see \code{\link{import.gsea.chip}}
#' @param file the fully specified path for the gct file
#' @param missing how would you like
#' @param verbose logical: print messages?
#' @param \dots further arguments passed to \code{\link{export.gsea.gct}}
#' @seealso \code{\link{export.gsea.gct}}
#' @author Mark Cowley, 2009-08-13
#' @export
export.gsea.gct.topTable <- function(tt, rma, chip, file, missing="", verbose=TRUE, ...) {
	if( !is.data.frame(chip) && is.character(chip) && file.exists(chip) )
		chip <- import.gsea.chip(chip)

	# probesetids <- intersect(tt$ID, chip[,1])
	# tt <- tt[tt$ID %in% probesetids]
	# chip <- chip[chip[,1] %in% probesetids]
	
	if( all(c("Nprobes", "Possible.Probes") %in% colnames(tt)) ) {
		ttg <- tt
	}
	else {
		ttg <- collapse.topTable(tt, chip[,1:2], toupper=TRUE, only.genes=TRUE, verbose=verbose)
	}
	chip <- chip[match(ttg$ID, chip[,1]), ]
	gct <- rma[ttg$ID, ]
	rn <- ttg[,1]
	if( any(is.na(rn) | nchar(rn) == 0 | is.null(rn) | rn == "---") ) {
		idx <- is.na(rn) | nchar(rn) == 0 | is.null(rn) | rn == "---"
		rn[idx] <- ttg$ID[idx]
	}
	rownames(gct) <- make.names(rn)

	export.gsea.gct(gct, chip[,3], file, missing=missing, ...)
}
