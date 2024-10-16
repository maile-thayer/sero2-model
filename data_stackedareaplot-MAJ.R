require(tidyverse)

data_stackedareaplot2 <- function(mats = list(), nsims=NULL, tmax=NULL) {
  if (is.null(nsims)) nsims <- ncol(mats[[1]])
  if (is.null(tmax)) tmax <- nrow(mats[[1]])
  medians <- tibble()
  for (i in 1:length(mats)) {
    medians <- bind_rows(medians, 
      tibble(
        time = 1:nrow(mats[[1]]),
        group = names(mats)[i],
        median = apply(mats[[i]], 1, median)))
  }
  
  medians <- medians %>%
    group_by(time) %>%
    mutate(value_perc = median / sum(median)) %>%
    ungroup() %>%
    mutate(group = factor(group, levels = names(mats)))
  
  return(medians)
}


