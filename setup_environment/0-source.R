library(dplyr)
library(cowplot)
library(ggplot2)
library(GenomicRanges)
theme_set(theme_cowplot())
options(ignore.interactive=TRUE) # always print progress bar.

.pipeline_ddir = "~/Dropbox (ICR)/EPICC_MENDELEY"
if (!file.exists(.pipeline_ddir)) {
  .pipeline_ddir = "../EPICC_data_mendeley"
}
stopifnot(file.exists(.pipeline_ddir))

if (!exists(".sourced_common")) {
  
  "./setup_environment" %>%
    list.files("[1-9]+.*[.]R$", full.names = TRUE) %>% 
    (function(x) x[!grepl("[.]not[.]R$", x)]) %>% # exclude ".not.R" files
    sapply(base::source)
  
  .sourced_common = TRUE
}
