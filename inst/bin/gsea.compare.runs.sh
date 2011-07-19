#!/bin/bash
# 
# For the scenario where 2 (X) rnk files are each run against 10 gmt files, there will be 2 top-level 
# dirs each containing (for eg) 10 (Y) dirs (c1_all, c2_all, c2_kegg, ....), this program combines the
# results from each gmt comparison to create 10 (Y) files, each with the 2 (X) sets of NES and FDR scores
# 
# Essentially, this program runs 'gsea.combine.runs.1gmt' for all gmt's in a set of dirs.
# Caveats: the list of gmt files used is determined from the 1st argument, and all other runs
#	must have used the same set of gmt files.
#
# Usage:
#   gsea.compare.runs.sh <run1> <run2> [<run3 ...] [output.xls]
#  eg gsea.compare.runs.sh male0vs8 male0vs4 female0vs8 male0vs8
#   ^^^^ where each argument is the top-level dir containing 15 folders for each of 15 gmt files.
#  if the last argument is not a directory, then it will be used to store the output file in,
# otherwise, ./gsea.combo.xls will be used.
# 
# Creates one xls file per gmt, and a master xls file (called ./gsea.combo.xls) with 1 sheet per gmt
#  
# Mark Cowley, 2009-03-20
# 2009-12-16: minor updates to find command on line 67. tested to work with the new R-based gsea.comare.runs.1gmt.sh
#


DEBUG=0
# DEBUG=1


if [ $# -lt 2 ]; then
	echo "Usage:
	$0 <run1> <run2> [<run3> ...] [output.xls]
	which equates to the following combinations:
	$0 <run1> <run2>
	$0 <run1> <run2> output.xls
	$0 <run1> <run2> [<run3> ...]
	$0 <run1> <run2> [<run3> ...] output.xls
	where run1 is a dir that contains a number of GSEA runs, where a run is 1rnk vs 1gmt
	"
	exit
fi


OUT=${!#}
if [ -d "$OUT" ]; then
	echo "Writing results to ./gsea.combo.xls" 
	OUT="gsea.combo.xls"
fi
if [ -f "$OUT" ]; then
	rm "$OUT"
fi
touch "$OUT"



dirs=$@
if [ $DEBUG -eq 1 ]; then echo "arguments: $# ${dirs[0]}"; fi

#
# If there's an index.html inside the 1st argument, then this program is being
# used from one level too low.
#
if [ -f $1/index.html ]; then
	echo "You probably want to run gsea.compare.runs.1gmt.sh"
	exit 3
fi


tmpdir=`mktemp -d -t gseacmp.XXXXX` # works on solaris and Mac OSX, even though -t works differently in both OS's
subdirs=$( find -s "$1" -type d -name "*GseaPreranked*" -maxdepth 1 -print | xargs basename | sed 's/\..*//' )


if [ $DEBUG -eq 1 ]; then 
	echo "tmpdir is: $tmpdir"
	echo "subdirs: $subdirs length ${#subdirs[@]}"
fi
	
i=0
for gmt in $subdirs; do 
	if [ $DEBUG -eq 1 ]; then
		echo -n "for $gmt, "
		echo -n "i=$i  --  "
	fi 
	let "i+=1"
	echo -n "."
	
	# I can't work out why find fails to find any files IF ${gmt}.xls exists
	rm -f ${gmt}.xls
	files=$( find -s $dirs -type d -name ${gmt}* -maxdepth 1 -print )
	if [ $DEBUG -eq 1 ]; then 
		echo "I found these files: ${files[@]}"
		echo "Command to be executed:  gsea.compare.runs.1gmt.sh ${files[@]} \"${tmpdir}/${gmt}.xls\""
	fi
	# gsea.compare.runs.1gmt.sh "${files[@]}" "${tmpdir}/${gmt}.xls"
	# gsea.compare.runs.1gmt.sh "$files" "${tmpdir}/${gmt}.xls"
	gsea.compare.runs.1gmt.sh ${files[@]} "${tmpdir}/${gmt}.xls" # WARNING: DON'T put the ${files[@]} in double quotes! bash only passes the first name on as an argument - the others are on 2nd++ lines and are forgotten
	retval=$?
	if [ $retval -ne 0 ]; then
		echo "gsea.compare.runs.1gmt.sh command failed for $gmt"
		exit
	fi
	if [ $DEBUG -eq 1 ]; then
		if [ -f "${tmpdir}/${gmt}.xls" ]; then
			echo "${gmt}.xls written to tmp"
		else
			echo "${gmt}.xls NOT written to tmp"
		fi
	fi
done
echo ""



if [ $DEBUG -eq 1 ]; then 
	echo "contents of tmpdir:"
	ls ${tmpdir}
	for file in $tmpdir/*; do
		echo $file
	done
fi

echo -n "Converting to tab delimited xls"
for f in ${tmpdir}/*.xls; do
	xls2tab $f &> /dev/null
	echo -n "."
done
echo ""

if [ $DEBUG -eq 1 ]; then 
	echo "contents of tmpdir:"
	ls ${tmpdir}
	for file in $tmpdir/*; do
		echo $file
	done
fi

echo -n "Saving xls"
tab2xls ${tmpdir}/*.tab "$OUT" &> /dev/null
echo ""

if [ $? -eq 0 ]; then
	echo ""

	echo -n "Generating PDF"
	gsea.compare.runs.plot.sh "$OUT"
	echo ""
fi
