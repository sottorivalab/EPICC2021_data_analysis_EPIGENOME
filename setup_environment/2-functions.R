null = # source all R files in function dir
  .function_source_file_dir %>%
    list.files("[.]R$",  full.names=TRUE) %>% 
    (function(x) x[!grepl("[.]not[.]R$", x)]) %>% # exclude ".not.R" files
    sapply(function(x) {source(x)})

