#!/bin/bash
# Fix the html output from gsea, or gseaPreranked.
# This annotates the gene tables with proper gene symbols, descriptions and hyperlinks,
# and (in the case of gseaPreranked output which does not have heatmaps maps), if a gct
# file has been made with one row per genesymbol, then heatmaps of the leading edge genes 
# will be made too.
#
# Mark Cowley, 2009-08-19
#

# DEBUG=0
DEBUG=1

if [ $# -eq 0 ]; then
	echo "Usage: $0 <dir> [<gct.file>]
  where <dir> is either a gsea output directory (eg './TvsCon/c1_all.GseaPreranked.1250580478355'), or 
  a directory of GSEA dirs (eg './TvsCon'),
  and gct.file is an optional gct file, where rows are indexed by unique gene-symbols (used for making heatmaps)."
	exit 1
fi

if [ $DEBUG -eq 0 ]; then
	echo "program name: ${0}"
	echo "$# args specified"
fi


dir=$1
gct=${2-NULL}
if [ $DEBUG -eq 0 ]; then
	echo arguments $dir $gct
fi

script=`mktemp -t gseafix.XXXXX`.R

cat <<EOF > $script
require(metaGSEA, quietly=TRUE, warn.conflicts=FALSE)

stopifnot( file.exists("$dir") )
gct <- "$gct"
if( !file.exists(gct) )
	gct <- NULL

gsea.fix.output("$dir", gct, overwrite=TRUE)
# traceback()
EOF

R --vanilla --quiet < $script
