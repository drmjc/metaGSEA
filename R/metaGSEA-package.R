#' Meta-analysis of GSEA analyses, including intra- and inter-experiment
#' comparisons.
#' 
#' \code{metaGSEA} is a collection of \R and \code{bash} code to simplify the down-stream analysis of
#' GSEA results.  Efficient methods for importing and storing GSEA outputs,
#' comparing genesets within a single GSEA run, and comparing genesets between
#' different GSEA runs.  Useful visualisation is the key strength of this
#' package, utilising GSEA bar plots, correspondence at the top (CAT) plots,
#' and hierarchical clustering plots of geneset similarities.  Also the
#' ability to clean up GSEA output on \code{unix-alike} systems, by improving the
#' GSEA preranked html output, and the gene tables in each GSEA report by adding
#' the gene symbol, description and hyperlinks to Entrez Gene.
#' 
#' Background\cr
#' GSEA is a popular program from the Broad Institute for performing functional
#' analysis of gene expression data. It generates rich, interactive outputs, however
#' each GSEA run is limited to one biological comparison, vs 1 collection of genesets.
#' In more complex experimental designs, its easy to end up with many GSEA result sets,
#' and comparing between them all to find consistent genesets, or inconsistent genesets
#' across all comparisons becomes very challenging. 
#' Enter \code{metaGSEA}.
#' 
#' Leading edge genes\cr
#' One of the unique properties of GSEA is the concept of leading edge genes within a geneset.
#' It's very rare for \emph{all} genes within a geneset to be dramatically up- or down-regulated; by
#' way of example, if we identify a 
#' strongly up-regulated geneset, there's usually ~60\% of the genes that are up-regulated, 
#' ~30\% that are unchanged, and the remaining ~10\% are down-regulated. The leading edge
#' genes are those that contribute towards the geneset getting its up-regulated score.
#' It's quite possible for 2 GSEA runs to find the same geneset strongly up-regulated,
#' however with almost opposite usage of genes within their leading edges. 
#' If we were studying hypoxia via 2 different experimental systems, we might be pleased
#' to see that the hypoxic geneset is ranked #1 in 2 different microarray studies, however
#' knowing whether the hypoxic signature is driven by the same genes in the 2 systems is
#' extremely important.
#' The bulk of the uniqueness of metaGSEA is that it can compare between multiple GSEA runs
#' using these leading edge genes, thereby revealing more molecular detail than just comparing
#' 2 GSEA runs at the level of the geneset names.
#' 
#' Expected input types\cr
#' metaGSEA supports a number of usage modes, including: 
#' analysis of just 1 GSEA result,
#' analysis of multiple GSEA results generated using different comparisons to the same GMT file
#' (eg c2_all), or
#' analysis of multiple GSEA runs on the same comparison vs multiple GMT files. metaGSEA supports
#' GSEA and GseaPreRanked results.
#' 
#' Usage\cr
#' You can \code{\link[=import.gsea]{import}}, \code{\link[=gsea.filter]{filter}} and 
#' \code{\link[=export.gsea]{export}} GSEA results, allowing better control over what
#' data you send to an external GSEA visualisation tool, such as the LeadingEdge viewer tool in the
#' GSEA GUI, or via GenePattern. Filtering can be done on all the columns that you see in the
#' up/down-regulated geneset summaries, namely geneset name, and all the statistics.
#' 
#' If you have just 1 GSEA result, besides filtering it, you can assess the similarity between
#' genesets using \code{\link{plot.gsea.leadingedge}}, which compares genesets using the leading
#' edge genes, and creates a \code{\link[=plot.gsea.leadingedge.HCL]{Heirarchical clustering dendrogram}}, 
#' \code{\link[=plot.gsea.leadingedge.heatmap]{heatmap}}, \code{\link[=plot.gsea.leadingedge.barplot]{barplot}}, 
#' and \code{\link[=plot.gsea.leadingedge.adjmat]{adjacency matrix}}.
#' 
#' If you have multiple GSEA results compared to the \emph{same} GMT file (ie the geneset names will
#' overlap), then you can compare results using \code{\link{gsea.compare.runs.1gmt}}, then filter
#' those results via \code{\link{gsea.compare.runs.filter}}, and plot the [dis]similarities via:
#' \code{\link{plot.gsea.venn}}, \code{\link{plot.gseacmp.barplot}}, \code{\link[=plot.CAT.GSEA]{CAT plots}},
#' and by a \code{\link[=plot.gsea.leadingedge.HCL]{combined HCL plot}}.
#' 
#' If you have multiple GSEA resuls vs \emph{different} GMT files (ie the geneset names will mostly
#' not match), then in addition to performing the analyses described above on each individual result,
#' we find a multi-panel \code{\link[=plot.CAT.GSEA]{CAT plot}} to be quite useful.
#' 
#' EnrichmentMap\cr
#' The \code{EnrichmentMap} plugin for \code{Cytoscape} lets you load GSEA results. it then generates
#' networks of geneset similarity, just like the HCL plots that metaGSEA creates. However, \code{EnrichmentMap}
#' compares \emph{all} genes in the geneset, so you end you looking at the structure of the genesets as
#' they are stored in MsigDB, not as they relate to your own dataset. You can create a custom
#' leading-edge only genes version of the GMT file, using \code{\link{import.gsea.leadingedge}} (which
#' is automatically run when you do an \code{\link{import.gsea}}), followed by \code{\link{export.gsea.gmt}}.
#' 
#' \tabular{ll}{ Package: \tab metaGSEA\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2009-09-24\cr License: \tab GPL\cr LazyLoad: \tab yes\cr
#' }
#' 
#' @name metaGSEA-package
#' @aliases metaGSEA-package metaGSEA
#' @docType package
#' @author Mark Cowley
#' 
#' Maintainer: Mark Cowley <m.cowley@@garvan.org.au>
#' @seealso \code{\link{import.gsea}}, \code{\link{gsea.filter}}, \code{\link{export.gsea}}
#' @references \url{http://www.broadinstitute.org/gsea}
#' @keywords package
#'
roxygen()
