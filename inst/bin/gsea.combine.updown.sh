#!/bin/bash
#
# GSEA report parser
#
# This program takes the output from running GSEA PreRanked (notably the pos and neg tables)
# and combines them into a single tab delim file.
# Why? Each combined file contains all of the genesets that passed the size filter criteria,
# regardless of if they are up or down-regulated. This faciltates cross-GSEA-result comaprisons.
# 
# Given a GSEA results dir this:
# - combines the pos and neg reports,
# - adds a column to rhs called DIRECTION which is either up or down
# - removes columns 2 and 3 (what was a hyperlink to MSigDB and GeneSet details)
# - sorts by NES, descending
# 
# Usage:
#    gsea.combine.updown.sh [dir1 [dir2 ...]]
# eg gsea.combine.updown.sh     <<< uses the current working dir
# eg gsea.combine.updown.sh .   <<< uses '.', also the cwd
# eg gsea.combine.updown.sh GSEA/c[1-9]+ <<< repeats for all specified directories
#
# Output:
#  Creates a file called gsea_report_combined_<timestamp>.xls in the dir where the other output
#  files are.
# 
# Mark Cowley, 2009-03-19
#
#


# if there are more than 1 command line arguments, then run this script on each argument,
# where each argument should be a directory name
#
if [ $# -gt 1 ]; then
	for arg in $@; do
		"$0" "$arg"
		echo -n "."
	done
	echo ""
	exit 0
elif [ $# -eq 1 ] && [ -d "$1" ]; then
	cd "$1"
fi

#
# at this point we should be inside a GSEA directory with the pos and neg output files
#
UP=`ls gsea_report_for_na_pos_?*.xls`
DOWN=`ls gsea_report_for_na_neg_?*.xls`
if [ ! -f "$UP" ]; then
	echo "No GSEA reports found"
	exit 1
fi
OUT=${UP/for_na_pos/combined}

# echo "Writing combined results file to $OUT"

tmp=`mktemp`
cp $UP $tmp
perl -pi -e 's/\t$/\tup/' $tmp
# fix the header row so it doesn't end in 'up'
perl -pi -e 's/EDGE\tup/EDGE\tDIRECTION/' $tmp
cp $tmp $OUT
rm $tmp

tmp=`mktemp`
cp $DOWN $tmp
perl -pi -e 's/\t$/\tdown/' $tmp
# skip the header row
tail +2 $tmp >> $OUT
# skip $tmp 1 >> $OUT
rm $tmp

#
# Sort by NES score from most up to most down.
#
header=`mktemp`
head -n1 $OUT > $header

tmp=`mktemp`
tail +2 $OUT > $tmp
# skip $OUT 1 > $tmp
# sort args: numerically, from hi-to-low, using tabs, by the 6th field, writing to tmp, reading from tmp
sort -n -r -t$'\t' -k 6 -o $tmp $tmp
cat $header > $OUT
cat $tmp >> $OUT
rm $tmp $header

# remove columns 2 and 3
tmp=`mktemp`
cut -d '	' -f 1,4-12 $OUT > $tmp
cp $tmp $OUT
# echo "Wrote file to $OUT"
