\name{convert.apt2gct}
\alias{convert.apt2gct}
\title{script to convert RMA & DABG objects into gct & res files, for all...}
\usage{convert.apt2gct(rma.apt.file="./normalised/rma.summary.txt",
    dabg.apt.file="./normalised/dabg.summary.txt", clm.file,
    chip.file=c("/pwbc/private/data/GSEA/chips/HuGene_1_0_st_v1.na30.hg19.chip",
    "/pwbc/private/data/GSEA/chips/MoGene_1_0_st_v1.na30.1.mm9.chip",
    "/pwbc/private/data/GSEA/chips/RaGene_1_0_st_v1.na30.1.rn4.chip")[1],
    out.dir=ifelse(file.exists(clm.file), dirname(clm.file), ""),
    create.gct=TRUE, create.res=TRUE, create.cls=TRUE,
    keep.controls=FALSE, all.genes=TRUE, expressed.genes=TRUE,
    detection.threshold=1e-05)
}
\description{script to convert RMA & DABG objects into gct & res files, for all
non-control probes, and those expressed probes with Pcount>0}
\author{Mark Cowley, 2010-07-06}
\arguments{\item{rma.apt.file}{path to an rma file created by APT. This should work
with files created by other algorithms like PLIER or gcrma, but that
hasn't been tested. must be present.}
\item{dabg.apt.file}{path to a dabg file created by APT. must be present.}
\item{clm.file}{optional clm file describing which arrays to include, in
which order & how you'd like them to be renamed.}
\item{chip.file}{optional chip file describing each probeset. recommended
if you want to filter out the control probesets.}
\item{out.dir}{where to output the res/gct/cls files? Defaults to the dir
that the clm file lives in.}
\item{create.gct}{do you want gct files, and/or res files to be created?}
\item{create.res}{do you want gct files, and/or res files to be created?}
\item{create.cls}{using the supplied clm file, create a cls file?}
\item{keep.controls}{keep, or exclude the control probes?}
\item{all.genes}{generate res/gct files for all probes on the array (after
the optional control-probe filtering)}
\item{expressed.genes}{generate res/gct files for expressed probes on the
array (with at least one DABG p-val of <= detection.threshold)}
\item{detection.threshold}{p-value threshold for detecting expressed genes}
}

