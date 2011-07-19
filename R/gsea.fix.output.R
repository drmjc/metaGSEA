#' Fix the output from a GSEA run.
#' GSEA and GseaPreRanked both produce rich HTML reports with hyperlinks to 
#' annotation databaes, images and heatmaps. GseaPreRanked, misses out on some
#' of the features of the GSEA reports, and both GSEA and GseaPreRanked output
#' can be improved somewhat (IMHO).
#' 
#' Here's a list of the fixes, all to the HTML files:\cr
#' GSEA and GseaPreRanked:\cr
#' 	- Within the GeneSet column, in the geneset reports, remove the link to Affymetrix, and add hyperlinks to
#' 	      Entrez and GeneCards\cr
#' 	- Within the geneset reports, add a hyperlink to the genesets in the header section\cr
#' 	- Within the up/down-regulated summaries, add hyperlinks to MSigDB for all genesets, not
#' 	      just the top N that had individual reports made\cr
#' 	- Within the geneset reports, add a gene decscription into the appropriate column, by searching
#' 	      within the chip file listed in the rpt file\cr
#' 
#' 
#' GseaPreRanked only:\cr
#' 	- Within \sQuote{index.html}, instead of naming the phenotypes na_pos and na_neg, call 
#' them up- and down-regulated genesets\cr
#' 	- within the individual geneset reports, add a heatmap (if a GCT file is provided)\cr
#'
#' Control over which elements are fixed:\cr
#' \code{fix.posneg}: fix the 2 html files that summarise the stats for the
#'   up/down regulated genesets. This only makes sense for
#'   GseaPreRanked output, but no harm in leaving it \code{TRUE} for gsea output.\cr
#' \code{fix.index}: fix the index.html file? This only makes sense for
#'   GseaPreRanked output, but no harm in leaving it \code{TRUE} for gsea output.\cr
#'  \code{fix.reports}: fix the individual geneset reports? This controls whether the
#' geneset report fixes reported above, including adding heatmaps is done.\cr
#' 
#' You can import individual, or multiple GSEA runs via the \code{dir} parameter:\cr
#' 	- a single gsea directory, eg "./GSEA/pos6_versus_neg6/c1_all.Gsea.1252052322484"\cr
#' 	- the parent dir containing multiple GSEA runs, eg "./GSEA/pos6_versus_neg6" which will fix all GSEA dirs within it.\cr
#' 	- a single geneset report file, eg "./GSEA/pos6_versus_neg6/c1_all.Gsea.1252052322484/CHR1P21.html"\cr
#' 
#' If you ran GseaPreRanked, then you don't get any
#'   heatmaps. Specifying a gct file via the \code{genes.gct} parameter will
#'  make heatmaps for each of the individual geneset reports. 
#' NB: the gct file should have one row per
#'   gene symbol only, and the row name should be the gene symbol, NOT the probesetID
#' TODO: allow a probe-level gct file which matches the probe-level rnk file & which 
#' is generally easier to obtain a priori.
#' 
#' @param dir see details
#' @param genes.gct The path to a GCT file. see details.
#' @param fix.index logical: fix the index.html file? see details.
#' @param fix.posneg logical: Fix the pair of up-/down-regulated geneset summary tables? see details.
#' @param fix.reports logical: fix the individual geneset reports? see details
#' @param debug logical: print debugging messages?
#' @param overwrite.heatmap.images logical: overwrite already existing heatmap images?
#' @param verbose logical: print verbose messages?
#' @return none. it just fixes gsea outputs in situ
#' @author Mark Cowley, 2009-07-09
#' @export
gsea.fix.output <- function(dir, genes.gct=NULL, fix.index=TRUE, fix.posneg=TRUE, fix.reports=TRUE, debug=FALSE, overwrite.heatmap.images=FALSE, verbose=FALSE) {
	stopifnot( file.exists(dir) )
	if( debug ) verbose <- TRUE

	if( !is.gsea.dir(dir) && grepl("html$", dir) )
		gsea.fix.output.geneset.report(dir, genes.gct, debug=debug)
	else if( !is.gsea.dir(dir) && is.dir(dir) ) {
		subdirs <- dir(dir, full=TRUE)
		if( any(is.gsea.dir(subdirs)) ) {
			subdirs <- subdirs[is.gsea.dir(subdirs)]
			if( verbose ) cat(sprintf("Found %d GSEA subdirs.\n", length(subdirs)))
			for(subdir in subdirs) {
				if( verbose ) cat("Fixing", basename(subdir), "...")
				gsea.fix.output(subdir, genes.gct=genes.gct,
					debug=debug, fix.index=fix.index, fix.posneg=fix.posneg, fix.reports=fix.reports)
				if( verbose ) cat("\n")
			}
		}
		return(TRUE)
	}
	
	#
	# error checking
	#
	index.html <- file.path(dir, "index.html")
	if( !file.exists(index.html) )
		stop("Not a GSEA directory (it does not contain index.html).\n")


	if( fix.index )
		gsea.fix.output.index.html(index.html, debug=debug)
	
	if( fix.posneg ) {
		# if gsea preranked was run, these 2 files will exist:
		pos.html <- dir(dir, pattern="gsea_report_for_na_pos.*html", full=TRUE)
		neg.html <- dir(dir, pattern="gsea_report_for_na_neg.*html", full=TRUE)
		if( length(pos.html) == 0 || length(neg.html) == 0 ) {
			dirs <- dir(dir, pattern="gsea_report_for_.*html", full=TRUE)
			pos.html <- dirs[1]
			neg.html <- dirs[2]
		}
		gsea.fix.output.detailed.enrichment.results.html(pos.html, neg.html, debug=debug)
	}
	
	if( fix.reports ) {	
		addHeatMap <- FALSE
		#
		# make gct files for the top N genesets
		#
		if( !is.null(genes.gct) && file.exists(genes.gct) ) {
			gsea.genesets2gct(dir, genes.gct, leading.edge=FALSE, which.sets="best")
			if( verbose ) cat(" gct files made...")
			#
			# make heatmaps
			#
			gct.files <- dir(dir, pattern=".*gct$", full=TRUE)
			# gp.gct2heatmap(gct.files, fix.html=FALSE)
			for( gct.file in gct.files ) {
				genepattern.HeatMapImage(gct.file, method="local", format="png", overwrite=overwrite.heatmap.images)
			}
			if( verbose ) cat(" heatmaps made...")
			addHeatMap <- TRUE
		}


		#
		# import the chip file and determine the Description for each Gene on the chip.
		#
		sym2desc <- NULL
		try( {
			chip <- import.gsea.chip(dir)
			sym2desc <- chip[,c(2,3)]
			sym2desc <- sym2desc[!is.na(sym2desc[,1]), ]
			sym2desc <- collapse.rows(sym2desc, 1)
		}, silent=TRUE)

		#
		# fix the html files.
		#
		files <- dir(dir, pattern=".*\\.html$") # don't include the dirname path (yet)
		files <- setdiff(files, c("index.html", "neg_snapshot.html", "pos_snapshot.html", "heat_map_corr_plot.html"))
		files <- setdiff(files, grep("gsea_report", files, value=TRUE))
		files <- file.path(dir, files) # convert filenames to paths
	
		for(file in files) {
			if( debug ) cat(file)
			gsea.fix.output.geneset.report(file, sym2desc, addHeatMap=TRUE)
			if( debug ) cat("\n")
		}
		if( verbose ) cat(" individual reports fixed.\n")
	}
}



#' fix the index.html file
#' 
#' @param index.html the path to the index.html report
#' @param debug logical: debugging mode?
#' @author Mark Cowley
#' @@return none. fixes html files.
#' @seealso \code{\link{gsea.fix.output}}
#' @export
gsea.fix.output.index.html <- function(index.html, debug=FALSE) {
	stopifnot( length(index.html) == 1 && file.exists(index.html) && grepl("html$", index.html) )
	# patterns <- c(
	# "'s|<b>na_pos</b>|<b>AvsB</b>|g'",
	# "'s|<b>na_neg</b>|<b>BvsA</b>|g'",
	# "'s|<h4>Enrichment in phenotype: <b>na</b></h4>|<h4>Enrichment in phenotype: <b>AvsB</b></h4>|'",
	# "'s|<h4>Enrichment in phenotype: <b>na</b></h4>|<h4>Enrichment in phenotype: <b>BvsA</b></h4>|'",
	# "'s|<b>AvsB</b><i> versus </i><b>BvsA</b>|<b><i>A versus B</i></b>|'")
	patterns <- c(
	"'s|<h4>Enrichment in phenotype: <b>na</b></h4>|<h4>Up-regulated genesets</h4>|'",
	"'s|<h4>Enrichment in phenotype: <b>na</b></h4>|<h4>Down-regulated genesets</h4>|'",
	"'s|<h4>Enrichment in phenotype: <b>AvsB</b></h4>|<h4>Up-regulated genesets</h4>|'",
	"'s|<h4>Enrichment in phenotype: <b>BvsA</b></h4>|<h4>Down-regulated genesets</h4>|'",
	"'s|AvsB|na_pos|'",
	"'s|BvsA|na_neg|'",
	"'s|<b><i>A versus B</i></b>|<b>na_pos</b><i> versus </i><b>na_neg</b>|'")
	perl.oneliner(patterns, index.html, debug=debug)
}

#' fix the "Detailed enrichment results in html format" 2 pages.
#' 
#' @param pos.html the path to the positive enrichemnt html report
#' @param neg.html the path to the negative enrichemnt html report
#' @param debug logical: debugging mode?
#' @author Mark Cowley
#' @@return none. fixes html files.
#' @seealso \code{\link{gsea.fix.output}}
#' @export
gsea.fix.output.detailed.enrichment.results.html <- function(pos.html, neg.html, debug=FALSE) {

	stopifnot( length(pos.html) == 1 && file.exists(pos.html) && grepl("html$", pos.html) )
	stopifnot( length(neg.html) == 1 && file.exists(neg.html) && grepl("html$", neg.html) )

	# -- pos
	patterns <- c(
		"'s|<title>Report for na_pos|<title>Report for Up-regulated genesets|'",
		"'s|<title>Report for AvsB|<title>Report for Up-regulated genesets|'",
		"'s|Table: Gene sets enriched in phenotype <b>na<b>|Table: Up-regulated gene sets|'")
	perl.oneliner(patterns, pos.html, debug=debug)
	
	# -- neg
	patterns <- c(
		"'s|<title>Report for na_neg|<title>Report for Down-regulated genesets|'",
		"'s|<title>Report for [AB]vs[AB]|<title>Report for Down-regulated genesets|'",
		"'s|Table: Gene sets enriched in phenotype <b>na<b>|Table: Down-regulated gene sets|'")
	perl.oneliner(patterns, neg.html, debug=debug)
	
	#
	# The c5_bp|mf|cc webpages don't have hyperlinks to MSigDB...
	#
	.fixOne <- function(html) {
		toGeneCards <- function(x) {
			sprintf("<a href='http://www.broad.mit.edu/gsea/msigdb/cards/%s.html'>%s</a>", x, x)
		}
		toCell <- function(x) {
			paste("<td>",trim(x),"</td>", sep="")
		}
		xls <- sub("html", "xls", html)
		tmp <- read.delim(xls)
		
		# genesets <- tmp$NAME[nchar(tmp[,3]) > 1]
		genesets <- tmp$NAME
		from <- toCell(genesets)
		to <- toCell(toGeneCards(genesets))
		patterns <- sprintf('"s|%s|%s|"', from, to)
		perl.oneliner(patterns, html, debug=debug)
	}
	.fixOne(pos.html)
	.fixOne(neg.html)
}




#' fix an individual GSEA GeneSet html report.
#' 
#' @param htmlfile the path to the gsea report file
#' @param symbol2description a 2 column data.frame mapping from symbols to description
#' @param addHeatMap logical: generate & include a heatmap into the html page
#' @param debug logical: debugging mode?
#' @return none. fixes html files.
#' @author Mark Cowley
#' @seealso \code{\link{gsea.fix.output}}
#' @export
gsea.fix.output.geneset.report <- function(htmlfile, symbol2description, addHeatMap=TRUE, debug=FALSE) {
	# cmd <- paste("perl -pi~ -e 's|>([^<]+)</a></td><td></td><td></td>|>$1</a></td><td>$1<br><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=gene&term=${1}[sym]%20AND%209606[taxid]\">Entrez</a>, &nbsp <a href=\"http://genome-www5.stanford.edu/cgi-bin/SMD/source/sourceResult?option=Name&choice=Gene&organism=Hs&criteria=${1}\">Source</a></td><td>${1} Description</td>|g'",
	# shQuote(htmlfile))
	# system(cmd, intern=FALSE)

	# THIS IS FOR DEBUGGING - IT SHOULD BE DONE JUST ONCE IN THE CALLING FUNCTION, AND NOT HERE FOR EVERY SINGLE HTML REPORT.
	if( missing(symbol2description) || is.null(symbol2description) ) {
		dir <- dirname(htmlfile)
		sym2desc <- NULL
		try( {
			chip <- import.gsea.chip(dir)
			sym2desc <- chip[,c(2,3)]
			sym2desc <- sym2desc[!is.na(sym2desc[,1]), ]
			sym2desc <- collapse.rows(sym2desc, 1)
			}, silent=TRUE)
			symbol2description <- sym2desc
	}
	
	#
	# fix the tab delimited file
	# -- this fills in the description field
	# 
	txtfile <- sub("html$", "xls", htmlfile)
	# file.copy(txtfile, paste(txtfile, "~", sep=""))
	data <- read.delim(txtfile, na="null")
	if( all(is.na(data$GENE.SYMBOL)) ) {
		cols <- 1:grep("CORE.ENRICHMENT", colnames(data))
		data <- data[,cols]
	
		data$GENE.SYMBOL <- data$PROBE
		data$GENE_TITLE <- lookUp(data$PROBE, symbol2description, 1, 2, "first")
		# colnames(data) <- gsub("\\.", "", colnames(data))
		write.delim(data, txtfile, row.names=FALSE, col.names=TRUE)
	}
	# else it has already been fixed
	
	
	
	#
	# fix the HTML file
	# 
	
	# header fix 1:hyperlink the geneset name to MSigDB in the header.
	cmd <- "'s|<td>GeneSet</td><td>([^<]+)</td>|<td>GeneSet</td><td><a href=\"http://www.broad.mit.edu/gsea/msigdb/cards/$1.html\">$1</a></td>|'"
	perl.oneliner(cmd, htmlfile, debug=debug)
	
	# # header fix 2:
	# # na_pos -> AvsB
	# # na_neg -> BvsA
	# # <td>na_pos</td>
	# cmd <- paste("perl -pi -e",
	# 	"'s|<td>na_pos</td>|<td>AvsB</td>|g'",
	# 	shQuote(htmlfile))
	# # cat(cmd, "\n")
	# system(cmd, intern=FALSE)
	# cmd <- paste("perl -pi -e",
	# 	"'s|<td>na_neg</td>|<td>BvsA</td>|g'",
	# 	shQuote(htmlfile))
	# # cat(cmd, "\n")
	# system(cmd, intern=FALSE)
	
	
	#
	# some internal functions to create HTML table elements & hyperlinks to Entrez Gene and
	# Stanford Source, just like GSEA does.
	#
	fixGeneSymbol <- function(symbol) {
		toEntrez <- function(symbol) {
			sprintf("<a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=gene&term=%s[sym]%%20AND%%209606[taxid]'>Entrez</a>", symbol)
		}
		toSource <- function(symbol) {
			sprintf("<a href='http://genome-www5.stanford.edu/cgi-bin/SMD/source/sourceResult?option=Name&choice=Gene&organism=Hs&criteria=%s'>Source</a>", symbol)
		}
		paste(symbol,"<br>",toEntrez(symbol),",&nbsp",toSource(symbol), sep="")
	}
	toCell <- function(x) {
		paste("<td>",trim(x),"</td>", sep="")
	}
	toCoreEnrich <- function(x) {
		res <- ifelse(x=="<td>Yes</td>", "<td bgcolor='#CCFFCC'>Yes</td>", "<td>No</td>")
	}
	toItalics <- function(x) {
		sprintf("<i>%s</i>", x)
	}
	#
	# using the tab delimited file that we've read in as the template, start converting it into HTML.
	# important note re quotes:
	# normally, in a shell, you would go:   perl -pi -e 's|||' file
	# BUT in this instance, we need to search for a pattern with single quotes: <table cols='8' border='1'>
	# and you can't escape a single quote in the shell...
	# SO
	# I have to use:   perl -pi -e "s|||" file
	# which means all of the html elements need to use single quotes,
	# AND since some gene symbols & titles contain single quotes (but never double quotes?? i hope)
	# then these need to be replaced with their asvii values. (\047 for ' and \042 for "")
	#
	output <- data
	output$NAME <- NULL
	output$RANK.METRIC.SCORE <- round(output$RANK.METRIC.SCORE, 3)
	output$RUNNING.ES <- round(output$RUNNING.ES, 4)
	#
	# get rid of single quotes, double quotes & back quotes from gene symbols and title's. 
	# eg 6'-phosphodiesterase; SMARCAD1: SWI/SNF-rel...ox 1`
	#
	.bashSafeQuotes <- function(x) {
		x <- gsub("'", "\\047", x, fixed=TRUE)
		x <- gsub('"', "\\042", x, fixed=TRUE)
		x <- gsub('`', "\\140", x, fixed=TRUE)
		x <- gsub('@', "\\100", x, fixed=TRUE)
		x
	}
	output$PROBE <- .bashSafeQuotes(output$PROBE)
	output$GENE_TITLE <- .bashSafeQuotes(output$GENE_TITLE)
	output$GENE.SYMBOL <- .bashSafeQuotes(output$GENE.SYMBOL)

	output$PROBE <- toItalics(output$PROBE)
	output$GENE.SYMBOL <- fixGeneSymbol(output$GENE.SYMBOL)
	output <- cbind("i"=1:nrow(output), output)
	
	output <- colapply(output, toCell)
	output <- as.data.frame(output, stringsAsFactors=FALSE)
	
	# these 2 columns aren't wrapped with the <td> <\td> that the others do.
	output$CORE.ENRICHMENT <- toCoreEnrich( output$CORE.ENRICHMENT )
	output[,1] <- sprintf("<td class='lessen'>%d</td>", 1:nrow(output))

	colnames <- gsub("\\.", " ", colnames(output))
	colnames <- sprintf("<th class='richTable'>%s</th>", colnames)
	colnames[1] <- "<th class='richTable'>"
	
	output <- data.frame(lhs="<tr>", output, rhs="</tr>", stringsAsFactors=FALSE)
	output <- rowPaste(output, "")
	output <- c(colnames, output)

	#
	# set up the perl regex
	#
	leading.tag <- "<table cols='8' border='1'>"
	trailing.tag <- "<caption class='table'>"
	output <- c(leading.tag, output, trailing.tag)
	output <- paste(output, collapse="\n")

	cmd <- sprintf('"s|%s.*%s|%s|"', leading.tag, trailing.tag, output)
	if( debug ) cat(cmd, "\n")

	perl.oneliner(cmd, htmlfile, debug=debug)



	#
	# change the html to include the HeatMap (as a png)
	#
	png.file <- sub("html", "png", htmlfile)
	if( addHeatMap && file.exists(png.file) ) {
		png.file <- basename(sub("$", '\\$', png.file, fixed=TRUE))
		cmd <- sprintf('\'s|<br></body></html>|<br>\n<div class="image"><img name="heatmap" src="%s"><br><br><caption>Fig 3: Heatmap of leading edge genes</caption></div><br>\n</body></html>|\'', png.file)
		perl.oneliner(cmd, htmlfile, debug=debug)
	}

}


#####
# todo: convert this to pure R code - ie no perl oneliners in here!
#
# here's my first attempts for the simpler case of the individual reports:
# tmp.gsea.fix.output.geneset.report <- function(htmlfile, symbol2description, addHeatMap=TRUE, debug=FALSE) {
# 	#
# 	# fix the HTML file
# 	# 
# 	html <- readLines(htmlfile, warn=FALSE)
# 	split.tags <- c("<html>", "<body>", "<div", "</div>", "<caption>", "<tr>", "</body>", "</html>")
# 	for(tag in split.tags) {
# 		html[2] <- gsub(tag, paste("\n",tag,sep=""), html[2])
# 	}
# 	html <- c(html[1], strsplit(html[2], "\n")[[1]])
# 
# 	invisible(html)
# }	
# 	div.start.tags <- grep("<div", html)
# 	div.close.tags <- which(html == "</div>")
# 	
# 	richTable.start <- grep("<div class='richTable'><table cols='8'", fixed=TRUE, html)
# 	richTable.end <- grep("<div class='image'><img name='gset_rnd_es_dist'", fixed=TRUE, html) - 2
# 	
# 	writeLines(html[1:(richTable.start-1)], OUT)
# 	writeLines(richTable, OUT)
# 	writeLines(html[(richTable.end+1):length(html)], OUT)
# 	
# 	invisible(html)
# }
# # debug(tmp.gsea.fix.output.geneset.report)
# tmp <- tmp.gsea.fix.output.geneset.report(file2)
# nchar(tmp)
# tmp[1:10]
