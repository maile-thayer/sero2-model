require(tidyverse)

data_stackedareaplot2 <- function(mats = list(), nsims=NULL, tmax=NULL, times=NULL) {
  if (!is.matrix(mats[[1]])) {
    for (i in 1:length(mats)) {
      mats[[i]] <- as.matrix(mats[[i]], nrow=length(mats[[i]])) 
    }
  }
  if (is.null(nsims)) nsims <- ncol(mats[[1]])
  if (is.null(times)) {
    if (is.null(tmax)) times <- 1:nrow(mats[[1]])
    else times <- 1:tmax
  }
  
  medians <- tibble()
  for (i in 1:length(mats)) {
    medians <- bind_rows(medians, 
      tibble(
        time = times,
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


