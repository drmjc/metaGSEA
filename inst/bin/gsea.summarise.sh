#!/bin/bash
# Summarise a GSEA Preranked run from 1 rnk vs N gmt files.
# Creates a genuine xls file.
# 
# Usage:
#	gsea.summarise.sh dir [out.xls]
#	gsea.summarise.sh <dir>, which runs from specified dir and writes to gsea.summary.xls
#	gsea.summarise.sh <dir> <out.xls>, which runs from specified dir and writes to specified xls file
#	
# Mark Cowley, 2009-03-23
# 2009-08-19: updated to work with the new import.gsea, and cat <<EOF
# 2009-12-09: must supply at least the first directory argument. scripts that run with no arguments are too dangerous!
#
usage() {
	echo "Usage: gsea.summarise.sh <directory> [<out.xls>]"
}
if [ $# -eq 0 ]; then 
	usage
	exit 1
fi

# DIR=${1-`pwd`}
DIR=${1}
OUT=${2-gsea.summary.xls}
echo $DIR $OUT

script=`mktemp -t gsea.XXXX`.R # works on OSX (..../gsea.XXXX.zQ7Mkbq0) and SunOS (/tmp//tmp/gsea.aG4q)

cat <<EOF > $script
require(metaGSEA, quietly=TRUE, warn.conflicts=FALSE)
require(excelIO, quietly=TRUE, warn.conflicts=FALSE)
res <- import.gsea("$DIR")
summary <- gsea.summarise(res)
if( is.vector(summary) ) {
	summary <- as.data.frame(summary)
}
write.xls(summary, "$OUT", row.names='')
EOF

R --vanilla --quiet < $script
