#!/bin/bash
#
# Compare the results from running GSEA preranked N times for N different rnk files vs
# the same gmt file. This combines the results from each run into a wide table, keeping
# just the NES and FDR scores.
#
# Usage:
#    gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* [run3/c1_all* ...] [output.xls]
#    which equates to these possibilities
#    gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all*
#    gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* output.xls
#    gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* [run3/c1_all* ...]
#    gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* [run3/c1_all* ...] output.xls
#    gsea.compare.runs.1gmt.sh $(find . -name c1_all* -type d) output.xls
# WARNING, please use your shell's tab complete to expand the *, since if you leave the *
# bash may expand this to all of the files underneath the directory, not the dir itself!!
# 
# Output:
# If a file argument is not specified at the end, then gsea.compare.1gmt.xls will be written,
# otherwise the resulting tab delim file will be named as you specified.
# The file has GeneSet and Size as first 2 columns, then 3 columns per run, which are
# NES, FDR and DIRECTION. The name of the rnk file used to generate each run will appear
# in the column titles
# 
# Mark Cowley, 2009-03-20
# 2009-12-16. major updates. engine has been rewritten in R, this script is now just a wrapper to the R script, using commandArgs to pass the dir names into R.
#

DEBUG=0
# DEBUG=1 

if [ $DEBUG -eq 1 ]; then
	echo "$# args specified"
fi

if [ $# -lt 3 ]; then
	echo "Only $# arguments were specified. You nead at least 4:
Usage:
   gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* [run3/c1_all* ...] [output.xls]
Which equates to these possibilities
   gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all*
   gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* output.xls
   gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* [run3/c1_all* ...]
   gsea.compare.runs.1gmt.sh run1/c1_all* run2/c1_all* [run3/c1_all* ...] output.xls
Easiest approach is to use find:
   gsea.compare.runs.1gmt.sh $(find . -name c1_all* -type d) output.xls

WARNING, use your shell's tab complete to expand the *, since if you leave the *
bash may expand this to all of the files underneath the directory, not the dir itself!!"
	exit -1
fi
#
# If the last arg is not a directory, then this will be used for the output file.
#
OUT=${!#} # << gets the last cmdline arg
if [ -d "$OUT" ]; then
	echo "Writing results to ./gsea.compare.1gmt.xls"
	OUT="gsea.compare.1gmt.xls"
fi
if [ -f "$OUT" ]; then
	rm "$OUT"
fi
touch "$OUT"

# echo $@

script=`mktemp -t gseacmp.XXXXX`.R # works on solaris and MacOSX
cat <<EOF > $script
require( pwbc, quietly=TRUE )
require( metaGSEA, quietly=TRUE )
dirs <- commandArgs(trailingOnly=TRUE)
dirs <- dirs[is.dir(dirs)]
gsea.compare.runs.1gmt(dirs, outfile="$OUT")
traceback()
EOF

if [ $DEBUG -eq 1 ]; then
	echo "Script file is: $script"
	cat $script
fi

if [ $DEBUG -eq 1 ]; then
	R --vanilla --quiet --args $@ < "$script" 
else
	R --vanilla --quiet --args $@ < "$script" &> /dev/null
fi
# 
# tmp=`mktemp`
# firstFile=1
# for dir in $@; do
# 	if [ ! -d "$dir" ]; then
# 		break
# 	fi
# 	if [ $DEBUG -eq 1 ]; then
# 		echo "Processing $dir"
# 	fi
# 
# 	#
# 	# 1. for each GSEA run, combine the results, and copy to tmp
# 	#
# 	gsea.combine.updown.sh "$dir"
# 	file=$(find "$dir" -type f -name gsea_report_combined* -depth 1)
# 	if [ $DEBUG -eq 1 ]; then
# 		echo "Combined file is: $file"
# 	fi
# 	cp $file $tmp
# 	
# 	#
# 	# 2. determine the name of the rnk file used (via parsing the rpt file)
# 	#
# 	rpt=$(find "$dir" -type f -name *rpt -depth 1)
# 	if [ $DEBUG -eq 1 ]; then
# 		echo "Rpt file: $rpt"
# 	fi
# 	# grab the line that refers to the rnk file, then get the 3rd word, then strip out the full path
# 	rpt=$(grep -w rnk $rpt | cut -f3 | xargs basename)
# 	rpt=${rpt/.rnk}
# 	
# 	#
# 	# 3. sort the results by column 1, keeping the header row intact
# 	#
# 	in=$tmp
# 	out=`mktemp`
# 	hdr=`mktemp`; head -1 $in > $hdr; body=`mktemp`; tail +2 $in > $body
# 	sort -t$'\t' -k 1 -o $body $body
# 	cat $hdr $body > $out; rm $hdr $body
# 	cp $out $tmp
# 	
# 	#
# 	# 4. transpose the txt file, so then they can be combined using cat
# 	#
# 	transpose $tmp $tmp
# 	
# 	#
# 	# 5. fix the column headers so they indicate which rnk file was used,
# 	# and eliminate 'useless' columns
# 	#
# 	perl -pi -e "s/^([0-9a-zA-Z -_]+)\t/\$1.${rpt}\t/" $tmp
# 	if [ $firstFile -eq 1 ]; then
# 		tmp2=`mktemp`
# 		egrep "^NAME|^SIZE|^NES|^FDR|^DIRECTION" $tmp > $tmp2
# 		mv $tmp2 $tmp
# 
# 		perl -pi -e "s/NAME.${rpt}/NAME/" $tmp
# 		perl -pi -e "s/SIZE.${rpt}/SIZE/" $tmp
# 		let "firstFile=0"
# 	else
# 		tmp2=`mktemp`
# 		egrep "^NES|^FDR|^DIRECTION" $tmp > $tmp2
# 		mv $tmp2 $tmp
# 	fi
# 	
# 	#
# 	# 6. grow the result file.
# 	#
# 	cat $tmp >> $OUT
# done
# 
# #
# # 7. transpose the result file. now it's the correct orientation
# #
# transpose $OUT $OUT
# 
# #
# # sort by the first files NES score (column 3)
# #
# in=$OUT
# out=`mktemp`
# hdr=`mktemp`; head -1 $in > $hdr; body=`mktemp`; tail +2 $in > $body
# sort -t$'\t' -k 3nr -o $body $body
# cat $hdr $body > $out; rm $hdr $body
# cp $out $OUT
# 
# tmp=`mktemp`
# cp "$OUT" $tmp
# tab2xls $tmp "$OUT" &> /dev/null
# # echo open $OUT
# 
# # if [ $DEBUG -eq 1 ]; then
# # 	echo "Making plot"
# # fi
# # gsea.compare.runs.1gmt.plot.sh "$OUT"
