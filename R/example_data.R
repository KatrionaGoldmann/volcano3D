#' PEAC synovial sample data
#'
#' A dataset containing sample data for 81 synovial biopsies from the PEAC 
#' cohort
#'
#' @format A data frame with 81 rows and 4 variables:
#' \describe{
#'   \item{ID}{IDs which match column names in `syn_example_rld``}
#'   \item{Pathotype}{The sample pathotype}
#' }
#' @source \url{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1
#' }
"syn_example_meta"

#' PEAC synovial gene expression data
#'
#' A dataset containing the gene expression data for 81 synovial biopsies from 
#' the PEAC cohort
#'
#' @format A data frame with 100 rows representing the most significant 
#' genes/probes and 81 columns representing samples. 
#' @source \url{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1
#' }
"syn_example_rld"

#' Synovial differential expression of genes across pathotypes in PEAC cohort 
#'
#' A dataset containing the differential expression parameters between different
#' pathotype groups for 81 synovial biopsies from the PEAC cohort.
#'
#' @format A data frame with 100 rows representing the most significant 
#' genes/probes and 13 columns for each statistical parameter:
#' \describe{
#'     \item{Module}{The gene name}
#'     \item{Fibroid-Lymphoid pvalue}{pvalue from fibroid vs lymphoid 
#'     comparison}
#'     \item{Fibroid-Lymphoid padj}{adjusted pvalue from fibroid vs lymphoid 
#'     comparison}
#'     \item{Fibroid-Lymphoid log2FoldChange}{logarithmic fold change from 
#'     fibroid vs lymphoid comparison}
#'     \item{Lymphoid-Myeloid pvalue}{pvalue from lymphoid vs myeloid 
#'     comparison}
#'     \item{Lymphoid-Myeloid padj}{adjusted pvalue from lymphoid vs myeloid 
#'     comparison}
#'     \item{Lymphoid-Myeloid log2FoldChange}{logarithmic fold change from 
#'     lymphoid vs myeloid comparison}
#'     \item{Myeloid-Fibroid pvalue}{pvalue from myeloid vs fibroid comparison}
#'     \item{Myeloid-Fibroid padj}{adjusted pvalue from myeloid vs fibroid 
#'     comparison}
#'     \item{Myeloid-Fibroid log2FoldChange}{logarithmic fold change from 
#'     myeloid vs fibroid comparison}
#'     \item{LRT pvalue}{pvalue from three-way likelihood ratio comparison}
#'     \item{LRT padj}{adjusted pvalue from three-way likelihood ratio 
#'     comparison}
#'     \item{LRT log2FoldChange}{logarithmic fold change from three-way 
#'     likelihood ratio comparison}
#' }
#' @source \url{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1
#' }
"syn_example_p"
