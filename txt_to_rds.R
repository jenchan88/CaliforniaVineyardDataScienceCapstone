lines <- readLines("C:/Users/tyler/OneDrive/Desktop/vineyard-data-capstone/module_text_files/herbicideNH_H.txt")

parsed_list <- list()
current_key <- NULL

for (line in lines) {
  if (grepl("^\\$`", line)) {
    current_key <- gsub("[^0-9]", "", line)
    parsed_list[[current_key]] <- character()
  } else if (grepl('^\\[1\\]', line)) {
    # Clean and extract values
    values <- unlist(strsplit(line, '"'))[c(FALSE, TRUE)]
    parsed_list[[current_key]] <- c(parsed_list[[current_key]], values)
  }
}

saveRDS(parsed_list, "HerbicideNH_H.rds")
