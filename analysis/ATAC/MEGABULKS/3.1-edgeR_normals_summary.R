source("setup_environment/0-source.R")
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(magrittr)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_summary_normal_glands"
input_dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_results_plus_cn"
regex_result_files="-normal_gland"

cpm_cutoff = 0
fc_cutoff = 0

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
annodb = "org.Hs.eg.db"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load edger test results ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

test_results_all =
  input_dat_dir %>%
  list.files(regex_result_files, full.names = TRUE, recursive = TRUE) %>% 
  magrittr::set_names(gsub("[.]rds", "", basename(.))) %>% 
  lapply(readRDS) %T>%
  saveRDS(file.path(dat_dir, "test_results_list.rds"))

fix_cn = function(x) {
  for (el in c("p","p_no_cn","max_cpm","fc","sig","cpm_baseline")) {
    colnames(x[[el]]) = 
      gsub(
        "all_normals_vs_C[0-9]+-", 
        "normal_vs_", 
        colnames(x[[el]])
      )
  }
  
  return(x)
}

analysis_results = 
  test_results_all %>% 
  test_results_to_analysis_set() %>% 
  add_peak_annotation_to_analysis_set() %>% 
  add_tss_overlap_to_analysis_set() %>% 
  add_genhancer_overlap_to_analysis_set() %>%
  fix_cn() %T>%
  saveRDS(file.path(dat_dir, "analysis_results.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Recurrence dataset ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

rec_counts = 
  analysis_results %>% 
  get_final_reccurence_count(min_fc=fc_cutoff, min_cpm=cpm_cutoff) %T>%
  saveRDS(file.path(dat_dir, "final_reccurence_data.rds"))

rec_counts_no_cn = 
  analysis_results %>% 
  get_final_reccurence_count(p_col = "p_no_cn", min_fc=fc_cutoff, min_cpm=cpm_cutoff) %T>%
  saveRDS(file.path(dat_dir, "final_reccurence_data_no_cn.rds"))
