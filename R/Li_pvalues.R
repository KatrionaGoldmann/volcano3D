#' Li Module expression data in Synovium
#'
#' Modular expression for synovial biopsies from the PEAC cohort
#' Li module expression data for synovial 81 samples
#'
#' A dataset containing the module expression of Li modules among 81 synovial 
#' PEAC biopsies. Calculated using QuSAGE
#'
#' @format A data frame with 346 rows representing modules and 81 columns 
#' representing samples
#' @source \url{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}
"syn_mod"

#' Synovial differential expression of modules across pathotypes in PEAC cohort 
#'
#' A dataset containing the relative differential expression parameters 
#' comparing module expression between different pathotype groups for 81 
#' synovial biopsies from the PEAC cohort.
#'
#' @format A data frame with 346 rows representing modules and 10 columns for 
#' each statistical parameter:
#' \describe{
#'     \item{Module}{The gene name}
#'     \item{Lymphoid-Fibroid p.value}{pvalue from lymphoid vs fibroid 
#'     comparison}
#'     \item{Lymphoid-Fibroid q.value}{adjusted pvalue from lymphoid vs 
#'     fibroid comparison}
#'     \item{Lymphoid-Fibroid logFC}{logarithmic fold change from lymphoid vs 
#'     fibroid comparison}
#'     \item{Lymphoid-Myeloid p.value}{pvalue from lymphoid vs myeloid 
#'     comparison}
#'     \item{Lymphoid-Myeloid q.value}{adjusted pvalue from lymphoid vs 
#'     myeloid comparison}
#'     \item{Lymphoid-Myeloid logFC}{logarithmic fold change from lymphoid vs 
#'     myeloid comparison}
#'     \item{Myeloid-Fibroid p.value}{pvalue from myeloid vs fibroid 
#'     comparison}
#'     \item{Myeloid-Fibroid q.value}{adjusted pvalue from myeloid vs fibroid 
#'     comparison}
#'     \item{Myeloid-Fibroid logFC}{logarithmic fold change from myeloid vs 
#'     fibroid comparison}
#' }
#' @source \url{https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31007-1}
"syn_mod_pvalues"
