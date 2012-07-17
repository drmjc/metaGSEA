#!/bin/bash
#
# generate a PDF plot of all NES scores with FDR < 5%
#
# This uses the default settings provided in the gsea.utils.R file,
# specifically, the plot.gseacmp.barplot function
#
# Mark Cowley, 2009-03-23
#


DEBUG=0
# DEBUG=1

if [ $DEBUG -eq 1 ]; then
	echo "program name: ${0}"
	echo "$# args specified"
fi

if [ $# -ne 1 ]; then
	echo "Usage: gsea.compare.runs.1gmt.plot.sh <gsea.comp.xls>
	  where the only argument is an XLS produced by running
	  gsea.compare.runs.1gmt.sh"
	exit 1
fi

IN=${1}
PDF=${IN/xls/pdf}

if [ $DEBUG -eq 1 ]; then
	echo "inputfile: $IN"
	echo "outputfile: $PDF"
fi

script=`mktemp -t gseacmp.XXXXX`.R # works on solaris and MacOSX
cat <<EOF > $script
require(metaGSEA, quietly=TRUE, warn.conflicts=FALSE)
NROW <- 3
NCOL <- 4
gseacmp <- list(import.gsea.compare.runs("$IN", 1))
names(gseacmp) <- ''
try({
pdf.A4("$PDF")
par(mfrow=c(NROW,NCOL), mar=c(7,4,3,1)+0.1)
for( i in 1:length(gseacmp) ) {
	tmp <- gsea.compare.runs.filter(gseacmp[[i]], fdr=NULL, sort.by='avgabs')
	if( nrow(tmp) == 0 )
		next
	else
		# plot.gseacmp.barplot(tmp, rows=NULL, fdr.thresh=0.05, do.par=TRUE, sub=names(gseacmp)[i])
		plot.gseacmp.barplot(tmp, fdr.thresh=0.05, do.par=TRUE, sub=names(gseacmp)[i])
}
dev.off()
}) # end try
traceback()
EOF

if [ $DEBUG -eq 1 ]; then
	echo "Script file is: $script"
	cat $script
fi

if [ $DEBUG -eq 1 ]; then
	R --vanilla --quiet < "$script"
else
	R --vanilla --quiet < "$script" &> /dev/null
fi
