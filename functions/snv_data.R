epicc_vcf_loading_wrapper = function(file) {
  
  checkmate::assertFileExists(file)
  
  # load and modify vcf data:
  THmisc::load_vcf_file(file) %>%
    (function(x) {colnames(x) = gsub("-GRCh38$", "", colnames(x)); x}) %>% 
    (function(x) x[,!colnames(x) %in% .excluded_samples]) %>% 
    THmisc::insert_cnas(.sequenza_cna_data) %>%
    THmisc::insert_ccf(get_purity(colnames(.)))
  
}

