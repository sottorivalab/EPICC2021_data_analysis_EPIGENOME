source("setup_environment/0-source.R")
source("functions/atac.R")
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(magrittr)
library(edgeR)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

fig_dir = "analysis/ATAC/MEGABULKS/plots/edger_results_plus_cn"
dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_results_plus_cn"
dis_source = "common" # source of dispersion estimate

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
annodb = "org.Hs.eg.db"

dataset_file = "analysis/ATAC/MEGABULKS/datasets/counts/data_matrix.rds"
atac_purity_data = "created_datasets/genotyping_estimates_per_sample.rds"
set = "nucleosome_free"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Loading datasets...\n")
data = dataset_file %>% readRDS()
atacseq_purity = atac_purity_data %>% readRDS()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Calculating CN data of each peak...\n")

cn_table = data %>% 
  rownames %>% 
  unlist %>% 
  as("GRanges") %>% "-"(250) %>%
  THmisc::get_cnas(.cna_data, names(.cna_data)) %>% 
  magrittr::set_rownames(rownames(data))

annot = THmisc::annotation_from_barcode_epicc(colnames(cn_table))
wh_tumour = annot$tissue_type == "cancer"

avg_cn_peaks = 
  apply(
    X = cn_table,
    MARGIN = 1,
    FUN = tapply,
    INDEX = paste0(annot$patient, annot$region),
    mean, na.rm = TRUE
  ) %>% t()

avg_cn_peaks_per_cancer = 
  apply(
    X = cn_table[, wh_tumour],
    MARGIN = 1,
    FUN = tapply,
    INDEX = annot$patient[wh_tumour],
    mean, na.rm = TRUE
  ) %>% t()


cn_and_purity_estimates = 
  list(
    cn_per_region = avg_cn_peaks,
    cn_per_cancer = avg_cn_peaks_per_cancer,
    purity = atacseq_purity
  )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Starting analysis...\n")

# groupings to plot:
pairs_to_plot = list(
  NvC = c("normal", "cancer"),
  NvCp = c("normal","pure")
)

test_group_labels = c(
  "cancer" = "Cancer",
  "impure" = "Impure cancer",
  "pure" = "Pure cancer",
  "NA" = "Unknown purity",
  "00" = "Purity < 20%",
  "20" = "Purity > 20%",
  "40" = "Purity > 40%",
  "60" = "Purity > 60%"
)


# annotate data:
peak_annotation  = 
  data %>%
  rownames() %>% 
  as("GRanges") %>%
  ChIPseeker::annotatePeak(TxDb=txdb, annoDb=annodb)

tss_relationship = 
  peak_annotation@detailGenomicAnnotation$Promoter %>% 
  ifelse("proximal", "distal")

for (peak_group in unique(tss_relationship)) {
  
  # create matrix of input count data
  peak_chr = as.character(seqnames(peak_annotation@anno))
  wh_use = tss_relationship == peak_group & peak_chr %in% paste0("chr", 1:22)
  matrix = data[wh_use,]
  colnames(matrix) = gsub("_nucleosome_free", "", colnames(matrix))
  
  # get input data
  input_data = get_input_data(matrix[grepl("EPICC_C5", colnames(matrix)),], cn_and_purity_estimates)
  dge_lists = with(input_data, get_dgelist_objects(data = counts, fc_data = fc))
  
  # dump the input data somewhere
  ofile = file.path(dat_dir, paste0("dge_lists_", peak_group, "_", set, ".rds"))
  saveRDS(dge_lists, ofile)
  
  ofile = file.path(dat_dir, paste0("input_data_", peak_group, "_", set, ".rds"))
  saveRDS(input_data, ofile)
  
  # do statistical tests for each group
  for (group in colnames(dge_lists$unadj)) {
    
    cat(group, "\n")
    if (group %in% c("all_normals")) next() 
    pair = c("all_normals",  group)
    case = gsub("-.*", "", group)
    
    # check if outfile exists
    fname = paste0(set, "-", case, "-", peak_group, "-", paste0(pair, collapse="_vs_"), ".rds")
    out_file = file.path(dat_dir, set, peak_group, fname)
    if (file.exists(out_file)) next()
    
    # testing
    tr = do_tests(dge_lists, input_data, group, dis_source)
    dir.create(dirname(out_file), FALSE, TRUE)
    saveRDS(tr, out_file)
    
    # plot expected vs observed fc vs coef:
    plot_exp_vs_obs_fc = plot_atac_cn_adjustment(rownames(tr), tr$coef_excl_cn, log(tr$exp_fc), tr$cn, group)
    fname = paste0(case, "-", peak_group, "-", paste0(pair, collapse="_vs_"),".png")
    out_file = file.path(fig_dir, "expected_coef_vs_observed_coef", set, fname)
    ggsave(out_file, plot_exp_vs_obs_fc, width=10, height=8)
    
    # plot expected vs observed fc vs coef:
    plot_exp_vs_obs_fc_adj = plot_atac_cn_adjustment(rownames(tr), tr$coef, 0, tr$cn, group)
    out_file = file.path(fig_dir, "expected_coef_vs_observed_coef_adjusted", set, fname)
    ggsave(out_file, plot_exp_vs_obs_fc_adj, width=10, height=8)
    
    
  }
}
