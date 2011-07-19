#!/bin/bash
# from a root dir containing the results from 4 GSEA analyses,
# for curated, motif, positional and computational, add symlinks
# from the ./.........curated.........../index.html to ./curated.html
#

for f in *; do
	if [ -d "$f" ]; then
		if [ -f "$f/index.html" ]; then
			out="$f.index.html"
			ln -s "$f/index.html" "$out"
		fi
	fi
done
