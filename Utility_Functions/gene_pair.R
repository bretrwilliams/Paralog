make_pairs_alphabetical <- function(symbol1, symbol2, all_permutations = FALSE) {
  if (all_permutations) {
    res <- tibble(symbol1, symbol2 = list(symbol2)) %>% 
      unnest() %>% 
      filter(symbol1 != symbol2)
  } else {
    res <- tibble(symbol1, symbol2) %>% 
      filter(symbol1 != symbol2)
  }
  res %>% 
    mutate(gene_pair = ifelse(symbol1 < symbol2, 
                              paste0(symbol1, ";", symbol2), 
                              paste0(symbol2, ";", symbol1))) %>% 
    .$gene_pair
}

get_separate_cols_from_pairs <- function(genepair, label = "symbol", sep = ";", remove = TRUE) {
  genepair %>%
    enframe(.) %>% 
    select(gene_pair = value) %>%
    separate(gene_pair, paste0(label, 1:2), sep = sep, remove = remove)
}

get_single_from_pairs <-  function(genepair, ...) {
  get_separate_cols_from_pairs(genepair, ...) %>% gather() %>% .$value %>% sort() %>% unique()
}

extract_symbol <- function(x) gsub(" \\(.*$", "", x)

sort_pair <- function(s1, s2) {
  if (s1 < s2) return(list(symbol1 = s1, symbol2 = s2))
  return(list(symbol1 = s2, symbol2 = s1))
}