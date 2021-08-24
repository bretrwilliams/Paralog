# ------------- reproduce the old ceres calculation from unscaled data ---------
calc_unscaled_to_old_ceres <- function(gene_effect_unscaled, essential.genes, nonessential.genes) {
  # calc median ceres_unscaled per sample for nonessetial genes
  nonessential_median <- gene_effect_unscaled %>% 
    filter(gene %in% nonessential.genes$gene) %>%
    group_by(Broad_ID) %>%
    summarise(ceres_nonessential_unscaled_median = median(ceres_unscaled, na.rm = TRUE))
  
  # shift ceres_unscaled by median of the nonessential genes, then calc median ceres_unscaled per sample for essetial genes
  common.essentials_median <- gene_effect_unscaled %>% 
    filter(gene %in% essential.genes$gene) %>%
    inner_join(nonessential_median, by = "Broad_ID") %>%
    mutate(ceres_shifted = ceres_unscaled - ceres_nonessential_unscaled_median) %>%
    group_by(Broad_ID) %>%
    summarise(ceres_common.essentials_shifted_median = median(ceres_shifted, na.rm = TRUE))
  
  gene_effect_unscaled %>%
    inner_join(nonessential_median, by = "Broad_ID") %>%
    inner_join(common.essentials_median, by = "Broad_ID") %>%
    mutate(ceres = (ceres_unscaled - ceres_nonessential_unscaled_median)/abs(ceres_common.essentials_shifted_median)) %>%
    select(-ceres_nonessential_unscaled_median, -ceres_common.essentials_shifted_median)
}


calc_ceres_summary <- function(gene_effect, essential.genes, nonessential.genes) {
  cat("nonessential genes:\n")
  gene_effect %>% filter(gene %in% nonessential.genes$gene) %>% select_if(is.numeric) %>% summary() %>% print()
  cat("essential genes:\n")
  gene_effect %>% filter(gene %in% essential.genes$gene) %>% select_if(is.numeric) %>% summary()
}