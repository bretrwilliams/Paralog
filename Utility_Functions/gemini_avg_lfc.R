#' log fold change
#'
#' @description retrieve the median of individual log fold changes (LFC) from a gemini.model
#'
#' @param Model an object of class gemini.model
#' @param LFC_center function for average calculation
#' @return A matrix of log fold changes
#' 
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#'
#' @export
#' 
#' @examples
#' data("Model", package = "gemini")
#' lfc <- gemini_avg_lfc(Model)
gemini_avg_lfc <- function(Model, LFC_center = 'median') {
  stopifnot("gemini.model" %in% class(Model))
  anno <- Model$Input$guide.pair.annot
  colnames(anno) <- c("key", "g1", "g2")
  Model$Input$LFC %>% 
    as.data.frame() %>%
    rownames_to_column("key") %>%
    gather(cellline, lfc, -key) %>%
    left_join(anno, by = "key") %>%
    mutate(dep_pair = ifelse(g1 < g2, 
                             paste0(gsub(" \\(.*", "", g1), Model$pattern_join, gsub(" \\(.*", "", g2)),
                             paste0(gsub(" \\(.*", "", g2), Model$pattern_join, gsub(" \\(.*", "", g1)))) %>%
    group_by(dep_pair, cellline) %>%
    summarize(lfc = match.fun(LFC_center)(lfc)) %>%
    ungroup() %>%
    spread(cellline, lfc) %>%
    column_to_rownames("dep_pair") %>%
    as.matrix()
}