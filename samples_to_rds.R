sample_module_size <- function(module_size, n_repeats, taxa_list){

  random_sample <- taxa_list[sample(length(taxa_list), module_size)]
  
  samples_list <- vector("list", n_repeats)
  
  for (i in 1:n_repeats) {
    samples_list[[i]] <- taxa_list[sample(length(taxa_list), module_size)]
  }
  
  samples <- sapply(samples_list, module_function_search2)
  metrics <- unlist(samples[4,])
  
  return(
    list(samples = samples
    ))
}

for (i in 1:7) {
  sample <- sample_module_size(i, 300)
  saveRDS(sample, paste0("sample_module_size_", i, ".rds"))
}

