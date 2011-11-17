#
# script to convert RMA & DABG objects into gct & res files, for all non-control probes,
# and those expressed probes with Pcount>0
#
# Parameters:
#	rma.apt.file: path to an rma file created by APT. This should work with files created by other algorithms like PLIER or gcrma, but that hasn't been tested. must be present.
#	dabg.apt.file: path to a dabg file created by APT. must be present.
#	clm.file: optional clm file describing which arrays to include, in which order & how you'd like them to be renamed.
#	chip.file: optional chip file describing each probeset. recommended if you want to filter out the control probesets.
#	out.dir: where to output the res/gct/cls files? Defaults to the dir that the clm file lives in.
#	create.gct, create.res: do you want gct files, and/or res files to be created?
#	create.cls: using the supplied clm file, create a cls file?
#	keep.controls: keep, or exclude the control probes?
#	all.genes: generate res/gct files for all probes on the array (after the optional control-probe filtering)
#	expressed.genes: generate res/gct files for expressed probes on the array (with at least one DABG p-val of <= detection.threshold)
#	detection.threshold: p-value threshold for detecting expressed genes
#
# WARNING re detection.threshold:
#	yes 1e-05 seems a bit low, but DABG was designed for exon arrays with just 4 probes per
# probeset. the Gene ST arrays have ~25-30 probes which necessitates a more conservative
# threshold. In practice 1e-05 gives similar numbers of expressed genes as using the 'P'
# threshold on older 3' affy arrays.
#
# Mark Cowley, 2010-07-06


##' script to convert RMA & DABG objects into gct & res files, for all
##' non-control probes, and those expressed probes with Pcount>0
##' 
##' @param rma.apt.file path to an rma file created by APT. This should work
##'   with files created by other algorithms like PLIER or gcrma, but that
##'   hasn't been tested. must be present.
##' @param dabg.apt.file path to a dabg file created by APT. must be present.
##' @param clm.file optional clm file describing which arrays to include, in
##'   which order & how you'd like them to be renamed.
##' @param chip.file optional chip file describing each probeset. recommended
##'   if you want to filter out the control probesets.
##' @param out.dir where to output the res/gct/cls files? Defaults to the dir
##'   that the clm file lives in.
##' @param create.gct do you want gct files, and/or res files to be created?
##' @param create.res do you want gct files, and/or res files to be created?
##' @param create.cls using the supplied clm file, create a cls file?
##' @param keep.controls keep, or exclude the control probes?
##' @param all.genes generate res/gct files for all probes on the array (after
##'   the optional control-probe filtering)
##' @param expressed.genes generate res/gct files for expressed probes on the
##'   array (with at least one DABG p-val of <= detection.threshold)
##' @param detection.threshold p-value threshold for detecting expressed genes
##' @author Mark Cowley, 2010-07-06
##' @export
convert.apt2gct <- function(
	rma.apt.file="./normalised/rma.summary.txt",
	dabg.apt.file="./normalised/dabg.summary.txt",
	clm.file=NULL,
	chip.file=c(
		"/pwbc/private/data/GSEA/chips/HuGene_1_0_st_v1.na30.hg19.chip",
		"/pwbc/private/data/GSEA/chips/MoGene_1_0_st_v1.na30.1.mm9.chip",
		"/pwbc/private/data/GSEA/chips/RaGene_1_0_st_v1.na30.1.rn4.chip"
	)[1],
	out.dir=ifelse(file.exists(clm.file), dirname(clm.file), ""),
	create.gct=TRUE,
	create.res=TRUE,
	create.cls=TRUE,
	keep.controls=FALSE,
	all.genes=TRUE,
	expressed.genes=TRUE,
	detection.threshold=1e-05
	) {
	# chip.file <- "/pwbc/private/data/GSEA/chips/HuGene_1_0_st_v1.na30.hg19.chip"
	# chip.file <- "/pwbc/private/data/GSEA/chips/MoGene_1_0_st_v1.na30.1.mm9.chip"
	# chip.file <- "/pwbc/private/data/GSEA/chips/RaGene_1_0_st_v1.na30.1.rn4.chip"
	stopifnot( all(file.exists(rma.apt.file, dabg.apt.file)) )

	out.gct.all.file <- file.path(out.dir, sub(".clm", "_rma_all.gct", basename(clm.file)))
	out.res.all.file <- sub(".gct", ".res", out.gct.all.file)
	out.gct.expressed.file <- sub("_all", "_expressed", out.gct.all.file)
	out.res.expressed.file <- sub(".gct", ".res", out.gct.expressed.file)
	out.cls.file <- sub(".clm$", ".cls", clm.file)

	cat("Importing RMA & DABG data.\n")
	rma <- import.APT(rma.apt.file)
	dabg <- import.dabg(dabg.apt.file)

	if( file.exists(clm.file) ) {
		clm <- import.gsea.clm(clm.file)
		cat("Comparing RMA & DABG data to the clm file.\n")
		if( !all(colnames(rma) %in% clm$CEL) ) stop("Mismatch re array names in APT file and CLM file!!\n")
		if( !all(colnames(dabg) %in% clm$CEL) ) stop("Mismatch re array names in APT file and CLM file!!\n")
		rma <- rma[,clm$CEL]
		colnames(rma) <- clm$Sample
		dabg <- dabg[,clm$CEL]
		colnames(dabg) <- clm$Sample
		#
		# clm to cls
		#
		cat("Converting clm file to a cls file.\n")
		if( create.cls ) gsea.convert.clm2cls(clm.file, out.cls.file)
	}

	calls <- dabg2calls(dabg, Pthresh=detection.threshold)
	Pcount <- rowSums(calls=="P")

	if( file.exists(chip.file) ) {
		chip <- import.gsea.chip(chip.file)
		chip <- chip[match(rownames(rma), chip[,1]), ]; chip[,1] <- rownames(rma)
		desc <- chip2description(chip=chip)
	}
	else {
		desc <- rownames(rma)
		names(desc) <- rownames(rma)
	}

	#
	# remove the controls
	#
	if( !keep.controls ) {
		cat("Excluding control probes.\n")
		ids <- names(desc)[setdiff(1:length(desc), grep("->|^AFFX", desc))]
		rma <- rma[ids, ]
		dabg <- dabg[ids, ]
		calls <- calls[ids, ]
		Pcount <- Pcount[ids]
		desc <- desc[ids]
	}

	if( all.genes ) {
		if( create.gct ) export.gsea.gct(rma, description=desc, file=out.gct.all.file)
		if( create.res ) export.broad.res(rma, calls, description=desc, file=out.res.all.file, unlog=FALSE)
	}

	# png("GenePattern/barplot.Pcount.png", 1024, 640)
	# hist.int(Pcount, xlab="Pcount", ylab="Freq", main="Pcount histogram")
	# dev.off()

	#
	# remove Pcount==0
	#
	if( expressed.genes ) {
		ids <- names(Pcount)[Pcount>0]
		rma <- rma[ids, ]
		dabg <- dabg[ids, ]
		calls <- calls[ids, ]
		Pcount <- Pcount[ids]
		desc <- desc[ids]
		if( create.gct ) export.gsea.gct(rma, description=desc, file=out.gct.expressed.file)
		if( create.res ) export.broad.res(rma, calls, description=desc, file=out.res.expressed.file, unlog=FALSE)
	}

}
