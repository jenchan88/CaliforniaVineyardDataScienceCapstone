library(microeco)
library(stringr)

search_taxa <- function(taxa) {
  taxa_lower <- tolower(taxa)
  functions_with_taxa <- c()
  for (func in names(regex)) {
    taxa_list_lower <- tolower(regex[[func]])
    if (taxa_lower %in% taxa_list_lower) {
      functions_with_taxa <- c(functions_with_taxa, func)
    }
  }
  return(functions_with_taxa)
}

#Create Function Matrix 
create_function_matrix <- function(module, function_list){
  function_matrix <- matrix(0, nrow = length(module),  
                            ncol =  length(unique(function_list)),
                            dimnames = 
                              list(module,unique(function_list)))
  
  return(function_matrix)
}

# Shared Function Metric
shared_function_metric <- function(matrix) {
  cross <- tcrossprod(matrix)
  connected <- sum(cross[row(cross) > col(cross)])
  unconnected <- sum(diag(cross))
  
  return(connected/(connected + unconnected))
}

lower_level_search_taxa <- function(level, class) {
  taxa_list <- all_genus |>
    filter(Class == class) |>
    distinct({{level}}, .keep_all = TRUE) |>
    pull({{level}})
  
  functions <- c()
  for (t in taxa_list) {
    funcs <- search_taxa(t)
    functions <- c(functions, funcs)
  }
  return(unique(functions))
}

module_function_search <- function(module, t_level=Genus) {
  taxa_functions <- list()
  all_functions <- c()
  
  for (taxa in module) {
    func_list <- lower_level_search_taxa({{t_level}}, taxa)
    taxa_functions[[taxa]] <- func_list
    all_functions <- c(all_functions, func_list)
  }
  
  function_counts <- table(all_functions)
  
  # Compute score
  total_functions <- length(function_counts)
  shared_functions <- sum(function_counts > 1)
  score <- (shared_functions / total_functions) * 100
  
  # Compute similarity metric 
  function_matrix <- create_function_matrix(module, all_functions)
  
  for (i in seq_along(taxa_functions)) {
    function_matrix[i, taxa_functions[[i]]] <- 1
  }
  
  if (length(all_functions) == 0){ 
    metric <- 0 
  } else{
    metric <- shared_function_metric(function_matrix)
  }
  
  
  return(list(
    functions_by_class = taxa_functions,
    function_counts = function_counts,
    shared_percentage = score,
    metric = metric
  ))
}

sample_module <- function(module, n_repeats, taxa_list){
  
  module_size <- length(module)
  module_score <- module_function_search(module)$metric
  random_sample <- taxa_list[sample(length(taxa_list), module_size)]
  
  samples_list <- vector("list", n_repeats)
  
  for (i in 1:n_repeats) {
    samples_list[[i]] <- taxa_list[sample(length(taxa_list), module_size)]
  }
  
  samples <- sapply(samples_list, module_function_search2)
  metrics <- unlist(samples[4,])
  pval <- (sum(metrics >= module_score) / n_repeats)
  
  return(
    list(samples = samples,
         pval = pval
    ))
}

