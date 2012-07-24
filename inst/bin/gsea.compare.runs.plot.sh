#!/bin/bash
#
# generate a PDF plot of all NES scores with FDR < 5%
#

# DEBUG=0
DEBUG=1

if [ $DEBUG -eq 1 ]; then
	echo "program name: ${0}"
	echo "$# args specified"
fi

if [ $# -ne 1 ]; then
	echo "Usage: gsea.compare.runs.plot.sh <gsea.comp.xls>
	  where the only argument is an XLS produced by running
	  gsea.compare.runs.sh"
	exit 1
fi

IN=${1}
PDF=${IN/xls/pdf}
# echo $IN $PDF

if [ $DEBUG -eq 1 ]; then
	echo "inputfile: $IN"
	echo "outputfile: $PDF"
fi

script=`mktemp`.R
cat <<EOF > $script
require(metaGSEA, quietly=TRUE, warn.conflicts=FALSE)

gseacmp <- import.gsea.compare.runs("$IN")
pdf.A4("$PDF")
par(mfrow=c(3,4))
for( i in 1:length(gseacmp) ) {
	tmp <- gsea.compare.runs.filter(gseacmp[[i]], 0.05, 'any')
	if( nrow(tmp) == 0 )
		next
	else
		# plot_gseacmp.barplot(tmp, NULL, 0.05, do.par=TRUE, sub=names(gseacmp)[i])
		plot_gseacmp.barplot(tmp, fdr.thresh=0.05, do.par=TRUE, sub=names(gseacmp)[i])
}
dev.off()
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
