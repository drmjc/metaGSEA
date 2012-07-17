#' From a GCT file, create a HeatMap Image, using the jar from Broad's GenePattern.
#' The Broad's HeatMapViewer tool generates nice looking, customisable heatmaps,
#' IMHO nicer than anything that \R can give, except maybe \code{\link[gplots]{heatmap.2}}.
#' If you have a bunch of GCT files that need converting to heatmap images fast,
#' then this is the code for you.
#' 
#' You can either connect to GenePattern server, and submit HeatMapImage jobs to the
#' server, or directly execute the jar file locally,
#'
#' Shout out to Joshua Gould, Broad Institute who wrote the jar in the first
#' place.
#' 
#' @param in.file a single gct/res or odf file.
#' @param out.file if NULL, the outfile will be in the same dir as in.file,
#'   with the extenstion changed dependant on format
#' @param method \dQuote{local} or \dQuote{server}
#' @param gp.connection if \code{method="server"} then you should specify a
#'   connection to the server. see examples.
#' @param format picture format, one of: \dQuote{png}, \dQuote{tiff}, \dQuote{jpeg},
#'  \dQuote{eps}, \dQuote{bmp}. \dQuote{png} is default & recommended. \dQuote{tiff} 
#' look good for publication.
#' @param row.size see HeatMapImage help at Broad's GenePattern server.
#' @param column.size see HeatMapImage help at Broad's GenePattern server.
#' @param show.grid see HeatMapImage help at Broad's GenePattern server.
#' @param grid.color see HeatMapImage help at Broad's GenePattern server.
#' @param show.row.descriptions see HeatMapImage help at Broad's GenePattern server.
#' @param show.row.names see HeatMapImage help at Broad's GenePattern server.
#' @param rows.to.highlight see HeatMapImage help at Broad's GenePattern server.
#' @param row.highlight.color see HeatMapImage help at Broad's GenePattern server.
#' @param color.scheme see HeatMapImage help at Broad's GenePattern server.
#' @param color.palette see HeatMapImage help at Broad's GenePattern server.
#' @param use.color.gradient see HeatMapImage help at Broad's GenePattern server.
#' @param memory see HeatMapImage help at Broad's GenePattern server.
#' @param verbose logical: verbose messages?
#' @param overwrite overwrite existing images?
#' @return none. creates heatmaps.
#' @author Mark Cowley, 2009-08-19
#' @examples
#' # not run
#' ## 1) running locally:
#' # genepattern.HeatMapImage(in.file=<myfile.gct>, method="local", format="png")
#' ## 2) 
#' # library(GenePattern)
#' # gp.con <- gp.login("http://pwbc.garvan.unsw.edu.au/gp", "<username>", "<password>")
#' # genepattern.HeatMapImage(in.file=<myfile.gct>, method="server", gp.connection=gp.con, format="png")
#' @export
genepattern.HeatMapImage <- function(in.file, out.file=NULL, 
	method=c("local", "server"), gp.connection=NULL, 
	format="png", 
	row.size=16, column.size=16, show.grid="yes", grid.color="0:0:0", 
	show.row.descriptions="yes", show.row.names="yes", 
	rows.to.highlight="", row.highlight.color="255:0:0", 
	color.scheme="row normalized", color.palette="", use.color.gradient="no",
	memory="-Xmx512m",
	verbose=FALSE, overwrite=FALSE) {
	method=="local" || (method=="server" && !is.null(gp.connection)) || stop("Either run the code locally, or provide a valid gp.connection object")
	
	ext <- sub(".*\\.", "", in.file)
	if( ! ext %in% c("res", "gct", "odf") )
		stop("Invalid in.file type. must be res/gct/odf.\n")

	if( length(in.file) > 1 )
		stop("use a loop if you've got multiple input files.\n")

	method <- method[1]

	if( is.null(out.file) )
		out.file <- sub("gct$|res$|odf$", format, in.file)

	if( file.exists(out.file) && !overwrite ) {
		cat(sprintf("'%s' already exists. skipping.\n", out.file))
		return(FALSE)
	}

	
	if( method == "local" ) {
		jar <- file.path(.path.package('metaGSEA'), 'bin', 'heatmapimage-o.jar')
		stopifnot( file.exists(jar) )

		# tmp.in <- mktemp()
		# file.copy(in.file, tmp.in, TRUE)
		# tmp.out <- file.path(dirname(tmp.in))
		
		# cmd <- sprintf('java -Djava.awt.headless=true -jar "%s" "%s" "%s" %s -c%d -r%d -g%s -l%s -a%s -s%s -f"%s" -h%s -n"%s" -m"%s" -u%s', jar, in.file, out.file, format, row.size, column.size, show.grid, grid.color, show.row.descriptions, show.row.names, rows.to.highlight, row.highlight.color, color.scheme, color.palette, use.color.gradient)
		cmd <- sprintf("java %s -Djava.awt.headless=true -jar '%s' '%s' '%s' %s -c%d -r%d -g%s -l%s -a%s -s%s -f'%s' -h%s -n'%s' -m'%s' -u%s", memory, jar, in.file, out.file, format, row.size, column.size, show.grid, grid.color, show.row.descriptions, show.row.names, rows.to.highlight, row.highlight.color, color.scheme, color.palette, use.color.gradient)
		if( verbose )
			cat(cmd, "\n")
		system(cmd, intern=TRUE)
	}
	else if( method == "server" ) {
		require(GenePattern) || stop("required package 'GenePattern' is not installed")
		tmp.dir <- tempdir()

		HeatMapImage.result <- run.analysis(gp.connection, "urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00032:6", input.dataset=in.file, output=out.file, output.format=format, column.size=as.character(column.size), row.size=as.character(row.size), show.grid=show.grid, grid.color=grid.color, show.row.descriptions=show.row.descriptions, show.row.names=show.row.names, rows.to.highlight=rows.to.highlight, row.highlight.color=row.highlight.color, color.scheme=color.scheme, color.palette=color.palette, use.color.gradient=use.color.gradient)
		preprocess.out.files <- job.result.download.files(HeatMapImage.result, tmp.dir)
		file.copy(preprocess.out.files[[1]], out.file, overwrite=TRUE)
	}
}
