#!/bin/bash
#
# Following a typical GSEA run comparing 1 rnk to N gmt files, this function combines
# the 'top tables' from each of the N gmt files into a single xls file, where
# each sheet in the workbook contains both the enriched GSEA terms for both
# up and down regulated genes.
#
# Run this function from the top-level dir where the next level are all of the
# N output dirs, usually starting with c1_all, c2_all, ... and ending in a timestamp
# Alternatively, provide the single path to the top-level dir as the 1st arg to
# this program.
#
# This creates gsea.reports.xls in the top-level dir, and also a 'gsea_report_combined.<TIMESTAMP>.xls'
# file in each of the individual output dirs (see gsea.combine.updown.sh)
#
# Usage:
#	gsea2xls.sh
#	gsea2xls.sh .
#	gsea2xls.sh ~/projects/username/GSEA/run1
#
# Value:
#	creates gsea.reports.xls which is a genuine, multi-paged XLS file (not tab delimited)
#
# Mark Cowley, 2009-03-19
#

#
# check the cmd line arguments
#
if [ $# -eq 1 ]; then
	if [ ! -d "$1" ]; then
		echo "Command line arg 1 is not a directory"
		exit 1
	else
		cd "$1"
	fi
fi

#
# generate all of the gsea_report_combined_... files in each dir
#
dirs=$(ls . | grep ^c[1-9]_*)
gsea.combine.updown.sh $dirs

#
# Combine all of the reports into a single XLS.
# -- since the file names become the worksheet names, and these are too long,
# I clean those file names up prior to making the new xls file
#

# clean up the file names into a tmp dir
tmpdir=`mktemp -d`
# echo $tmpdir
reports=$(find . -regex .*/gsea_report_combined.* -print)
for report in $reports; do
	gmt=${report/.Gsea*}
	gmt=${gmt/.\/}
	# echo $gmt
	cp $report $tmpdir/${gmt}.txt
done
# combine the files into an XLS
reports=$(find $tmpdir -regex .*txt -print)
# echo $reports
OUT=gsea.reports.xls
OFS=" "
tab2xls $reports "./$OUT"
exit $?
