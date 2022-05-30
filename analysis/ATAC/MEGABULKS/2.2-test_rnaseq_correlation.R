source("setup_environment/0-source.R")
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(magrittr)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_summary"
#input_dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_results_plus_cn"
#dataset_file = "tmp_old_datafiles/megabulk_data.rds"
#atac_purity_data = "tmp_old_datafiles/genotyping_estimates_per_sample.rds"
min_rec = 6
n_cores = 40
n_bs = 10000

#txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
#annodb = "org.Hs.eg.db"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load datasets ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#megabulk_data = dataset_file %>% readRDS()
#atacseq_purity = atac_purity_data %>% readRDS()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Peak annotation ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#peak_annotation  = 
#  megabulk_data$data_matrix$no_nucleosome %>%
#  rownames() %>% as("GRanges") %>%
#  ChIPseeker::annotatePeak(TxDb=txdb, annoDb=annodb)

#tss_relationship = 
#  peak_annotation@detailGenomicAnnotation$Promoter %>% 
#  ifelse("proximal", "distal")

#wh_region = grepl("C.*region_[E]", colnames(megabulk_data$data_matrix$no_nucleosome))
#wh_not_sex_chr = !grepl("chr[XY]", rownames(megabulk_data$data_matrix$no_nucleosome))
#wh_promoter = tss_relationship == "proximal"
#cpm_data_normal = edgeR::cpm(megabulk_data$data_matrix$no_nucleosome[wh_not_sex_chr&wh_promoter,wh_region])
#mean_cpm_normal = apply(cpm_data_normal, 1, mean)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load edger test results ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

analysis_results = readRDS(file.path(dat_dir, "analysis_results.rds"))
rec_counts = readRDS(file.path(dat_dir, "final_reccurence_data.rds"))
rec_counts_no_cn = readRDS(file.path(dat_dir, "final_reccurence_data_no_cn.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# RNA-seq association ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# peaks to test:
rnaseq_test_results_rec_counts = 
  rec_counts$summary %>% 
  dplyr::filter(recurrence >= min_rec) %>% 
  arrange(-recurrence) %>%
  (function(x) split(x, (seq_len(NROW(x)) - 1) %% (n_cores * 20))) %>% 
  pbmcapply::pbmclapply(test_rnaseq_association_peaks, rec_counts, .dds_data, n_bs = n_bs, test_data=T, mc.cores=n_cores, verbose=FALSE) %>% 
  do.call(what=rbind)

out_file = file.path(dat_dir, "rnaseq_test_results_rec_counts.rds")
saveRDS(rnaseq_test_results_rec_counts, out_file)

out_file = file.path(dat_dir, "rnaseq_test_results_rec_counts.tsv.gz")
readr::write_tsv(rnaseq_test_results_rec_counts, out_file)
