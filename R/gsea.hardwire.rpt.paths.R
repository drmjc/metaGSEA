# Rewire the file locations referenced within a GSEA rpt file.
# Why? If you run a GSEA on a GenePattern server, you can end up with files in a rpt file that you can
# no longer access on the machine where you're analysing the GSEA results... Useful if you have a directory full of chip files, and gmx files.
#
# Parameters:
#	dir: a GSA directory. The rpt file within this dir will be modified.
#	chip.dir: the directory containing your GSEA .chip files
#	gmx.dir: the directory containing your GSEA .gmx/.gmt files
# 
# Value:
#	none. the rpt file will be modified.
#
# Mark Cowley, 2009-12-11
#


##' Rewire the file locations referenced within a GSEA rpt file.
##' 
##' Why? If you run a GSEA on a GenePattern server, you can end up with files
##' in a rpt file that you can
##' no longer access on the machine where you're analysing the GSEA results...
##' Useful if you have a directory full of chip files, and gmx files.
##' 
##' @param dir a GSA directory. The rpt file within this dir will be modified.
##' @param chip.dir the directory containing your GSEA .chip files
##' @param gmx.dir the directory containing your GSEA .gmx/.gmt files
##' @return none. the rpt file will be modified.
##' @author Mark Cowley, 2009-12-11
##' @export
gsea.hardwire.rpt.paths <- function(dir, chip.dir="/pwbc/data/GSEA/chips", gmx.dir="/pwbc/data/GSEA/gmt/v2.5") {
	rpt <- import.gsea.rpt(dir)
	rpt$chip <- file.path(chip.dir, basename(rpt$chip))
	rpt$gmx <- file.path(gmx.dir, basename(rpt$gmx))
	rpt.name <- dir(dir, pattern=".*\\.rpt", full=TRUE)

	export.gsea.rpt(rpt, rpt.name)
}
